import numpy as np
import pickle
import pandas as pd
import re
import os
import json
import MDAnalysis as mda
from MDAnalysis.analysis import hole2
from MDAnalysis.analysis import align
import sys



temp = sys.argv[1]

traj_dir = '../trajs'

u = mda.Universe('ion_channel.pdb', f'{traj_dir}/rep_{temp}.xtc')
# need to align the whole dcd trajectory to ref of 7mij_renumber.pdb
# get center of pore coordinates between resi 675 by taking average of their CA coordinates
# specify to search along z axis since that pdb (in pymol) is oriented that way.
ref = mda.Universe('ion_channel.pdb')


# make a selection to center the pore
alignment = align.AlignTraj(u, ref, select = 'resnum 652:679',
                            in_memory=True)
alignment.run()


if os.path.exists(f'hole_profiles_{temp}'):
    os.makedirs(f'hole_profiles_{temp}')
else:
    pass

center_resis = u.select_atoms('resnum 675 and name CA')

center_coords = np.mean(center_resis.positions,axis=0)

hole_exe = '/PATH/TO/SOFTWARE/hole2/exe/hole'

with hole2.HoleAnalysis(u, executable=hole_exe, cpoint=center_coords,
                        cvect=[0,0,1]) as h2:
    h2.run()


hole_profiles = {}
for i in range(len(h2.results['profiles'])):
    hole_profiles[i] = {'rxn_coord':h2.results['profiles'][i]['rxn_coord'], 'radius': h2.results['profiles'][i]['radius'], 'cen_line_D':h2.results['profiles'][i]['cen_line_D']}

with open(f'trp_temp_{temp}_hole_profiles.pk', 'wb') as f:
    pickle.dump(hole_profiles, f, protocol=pickle.HIGHEST_PROTOCOL)

if os.path.exists('hole_traj.dcd'):
    os.remove('hole_traj.dcd')
else:
    print("The file does not exist")
