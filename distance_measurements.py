import MDAnalysis as mda
from MDAnalysis.analysis.distances import dist
import pandas as pd
import json
import pickle
import pmda
from pmda.custom import AnalysisFromFunction

def measure_distance(ag1,ag2):
    return dist(ag1,ag2)[2][0]

def rename_df_indices(df, append_string):
    # rename indices for a dataframe

    indices = df.index
    new_indices = {}
    for index in indices:
        new_indices[index]=(f'{index}_{append_string}')
    return df.rename(mapper=new_indices)

# provide a file of mdanalysis selections (keys can be simple name of what you're measuring, values are mdanalysis selections)
selections = json.load(open('extended_distance_selections.json','r'))
# do the distance calculation
structure = '../7mij_renumber.pdb'


for i,trajectory_file in enumerate(trajectory_files):
    traj = f'{file}.xtc'
    u = mda.Universe(structure, traj)
    distances = {key:[] for key in selections.keys()}
    for key in selections.keys():
        analysis = AnalysisFromFunction(measure_distance,u,u.select_atoms(selections[key][0]),u.select_atoms(selections[key][1]))
        analysis.run()
        distances[key] = analysis.results
    ddf = pd.DataFrame(distances)
    ddf = rename_df_indices(ddf, temp)
    ddf.to_pickle(f'trajectory_{i}.pd')
