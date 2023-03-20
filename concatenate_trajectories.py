import MDAnalysis as mda
import os
import sys

construct = sys.argv[1]

run_out = 'run_output'

# might need a list comprehension with range() if you want to be sure trajectories are in original order
# can be a problem when the file indices get double digit if they don't start with (01,02...)
trajs = [f'{run_out}/{file}' for file in os.listdir(run_out) if file.endswith('.dcd')]

u = mda.Universe(f'{construct}_minimized.pdb', trajs)

# can get rid of solvent here
selection = u.select_atoms('all')

# change to xtc or smaller file type
with mda.Writer(f'gromacs/{construct}_concatenated.xtc', selection.n_atoms) as W:
    for ts in u.trajectory:
        W.write(selection)
