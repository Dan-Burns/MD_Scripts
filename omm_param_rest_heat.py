#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from openmm.app import *
from openmm import *
from openff.toolkit.topology import *
from openmm.unit import *
import openmmtools as omt
from openmmtools.integrators import *
import os
import numpy as np
import MDAnalysis as mda
import pandas as pd
import matplotlib.pyplot as plt
import parmed 
import re
from sys import stdout


# In[ ]:


construct = 'eeic'
ligand = 'pep'
structures = '../../input_structures/'
pdb = PDBFile(f'{structures}{ligand.upper()}/{construct}.pdb')
pdb_mg = PDBFile(f'{structures}{ligand.upper()}/{construct}_mg.pdb')
pdb_ligand = PDBFile(f'{structures}{ligand.upper()}/{construct}_{ligand}_conect.pdb')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.add(pdb_mg.topology, pdb_mg.positions)
modeller.add(pdb_ligand.topology, pdb_ligand.positions)

# parameterize akg
molecule = Molecule()
ligand_molecule = molecule.from_file(f'{structures}{ligand.upper()}/{ligand.upper()}.mol')
ligand_molecule.compute_partial_charges_am1bcc()
from openmmforcefields.generators import GAFFTemplateGenerator
gaff = GAFFTemplateGenerator(molecules=ligand_molecule)
# Create an OpenMM ForceField object with AMBER ff14SB and TIP3P with compatible ions
from simtk.openmm.app import ForceField
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
# Register the GAFF template generator
forcefield.registerTemplateGenerator(gaff.generator)
# You can now parameterize an OpenMM Topology object that contains the specified molecule.
# forcefield will load the appropriate GAFF parameters when needed, and antechamber
# will be used to generate small molecule parameters on the fly.

# solvate
modeller.addSolvent(forcefield, model='tip3p', padding=1*nanometer, ionicStrength=0.1*molar )
# create system object
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,nonbondedCutoff=1*nanometer, constraints=HBonds)
# define temperature and pressure
temperature = 310 * kelvin
pressure = 1 * bar
# Add pressure control
system.addForce(MonteCarloBarostat(pressure, temperature))
# create integrator object
integrator = LangevinMiddleIntegrator(temperature, 1/picosecond, 2*femtoseconds)
# create simulation object
simulation = Simulation(modeller.topology, system, integrator)



##### Save a gromacs topology for future trjconv use - Use a no-constraints version of system to avoid parmed error
parmed_system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME,nonbondedCutoff=1*nanometer, rigidWater=False)
pmd_structure = parmed.openmm.load_topology(simulation.topology, system=parmed_system, xyz=modeller.positions)

if os.path.exists('gromacs'):
    pass
else:
    os.makedirs('gromacs')
    
pmd_structure.save(f"gromacs/{construct}_{ligand}_SYSTEM.top", overwrite=True)
pmd_structure.save(f"gromacs/{construct}_{ligand}_SYSTEM.gro", overwrite=True)
##### 


# In[ ]:


simulation.context.reinitialize()
simulation.context.setPositions(modeller.positions)


# In[ ]:


print(f'Before adding restraints there are {len(system.getForces())} in the system.', flush=True)


# In[ ]:


############### Create Dihedral Restraint for PEP or SEP ################
# Get IDs of pep torsion atoms
# dihedral is C1 C2 O2 P
# checked and works for eeic with atoms 9729 9731 9732 9724  AND 9741 9743 9744 9736
if ligand == 'pep':
    topology = simulation.topology
    peps = {}
    for res in topology.residues():
        if res.name == 'PEP':
            peps[str(res.index)] = {'C1':0,'C2':0,'O2':0,'P':0}
            for atom in res.atoms():
                if atom.name in ['C1', 'C2', 'O2', 'P']:
                    peps[str(res.index)][atom.name]=atom.index
if ligand == 'sep':
    topology = simulation.topology
    seps = {}
    for res in topology.residues():
        if res.name == 'SEP':
            seps[str(res.index)] = {'C1':0,'C2':0,'O01':0,'S01':0}
            for atom in res.atoms():
                if atom.name in ['C1', 'C2', 'O01', 'S01']:
                    seps[str(res.index)][atom.name]=atom.index




# Add the dihedral restraint to both PEPs
if ligand == 'pep':
    restraint = PeriodicTorsionForce()
    system.addForce(restraint)
    for key in peps.keys():
        restraint.addTorsion(peps[key]['C1'], peps[key]['C2'], peps[key]['O2'], peps[key]['P'], 1, 
                             (270)*degrees, 1000*kilojoules_per_mole)
if ligand == 'sep':
    restraint = PeriodicTorsionForce()
    system.addForce(restraint)
    for key in seps.keys():
        restraint.addTorsion(seps[key]['C1'], seps[key]['C2'], seps[key]['O01'], seps[key]['S01'], 1, 
                             (270)*degrees, 1000*kilojoules_per_mole)

################################ END Dihedral Restraint Section ######################################################


############################ Create Gly 452 to PEP/SEP C2 distance restraint ##################
# Check which gly 452 is closest to each PEP C2 atom (can't use pdb resnums...) by finding the gly within 7 angstroms
# and get the C2 id

#gly_indices = {'eeic':452, 'eteic':452, 'teeic':452, 'teic':452}
gly_index = 452

if ligand == 'pep':
    pep_gly_pairs = []
    for pep in peps.keys():
        for res in topology.residues():
            if res.id == str(gly_index):
                for atom in res.atoms():
                    if atom.name == 'CA':
                        if np.linalg.norm(modeller.getPositions()[peps[pep]['C2']] - 
                                          modeller.getPositions()[atom.index]) < 0.7*nanometer:
                            pep_gly_pairs.append([peps[pep]['C2'],atom.index])

    pep_mg_pairs = []
    for pep in peps.keys():
        for res in topology.residues():
            if res.name == 'MG':
                for atom in res.atoms():
                    if atom.name == 'MG':
                        if np.linalg.norm(modeller.getPositions()[peps[pep]['C2']] - 
                                      modeller.getPositions()[atom.index]) < 0.7*nanometer:
                            pep_mg_pairs.append([peps[pep]['P'],atom.index])

if ligand == 'sep':
    sep_gly_pairs = []
    for sep in seps.keys():
        for res in topology.residues():
            if res.id == str(gly_index):
                for atom in res.atoms():
                    if atom.name == 'CA':
                        if np.linalg.norm(modeller.getPositions()[seps[sep]['C2']] - 
                                          modeller.getPositions()[atom.index]) < 0.7*nanometer:
                            sep_gly_pairs.append([seps[sep]['C2'],atom.index])

    sep_mg_pairs = []
    for sep in seps.keys():
        for res in topology.residues():
            if res.name == 'MG':
                for atom in res.atoms():
                    if atom.name == 'MG':
                        if np.linalg.norm(modeller.getPositions()[seps[sep]['C2']] - 
                                      modeller.getPositions()[atom.index]) < 0.7*nanometer:
                            sep_mg_pairs.append([seps[sep]['S01'],atom.index])
                


# In[ ]:


if ligand == 'pep':
    pep_gly_rest1 = omt.forces.FlatBottomRestraintForce(1000*kilojoules_per_mole, 5*angstrom,
                                                            restrained_atom_indices1=[pep_gly_pairs[0][0]],
                                                            restrained_atom_indices2=[pep_gly_pairs[0][1]],
                                                           controlling_parameter_name='pep_gly_1')

    pep_gly_rest2 = omt.forces.FlatBottomRestraintForce(1000*kilojoules_per_mole, 5*angstrom,
                                                            restrained_atom_indices1=[pep_gly_pairs[1][0]],
                                                            restrained_atom_indices2=[pep_gly_pairs[1][1]],
                                                           controlling_parameter_name='pep_gly_2')
    system.addForce(pep_gly_rest1)
    system.addForce(pep_gly_rest2)
'''
    pep_mg_rest3 = omt.forces.FlatBottomRestraintForce(1000*kilojoules_per_mole, 3.5*angstrom,
                                                            restrained_atom_indices1=[pep_mg_pairs[0][0]],
                                                            restrained_atom_indices2=[pep_mg_pairs[0][1]],
                                                           controlling_parameter_name='pep_mg_3')
    pep_mg_rest4 = omt.forces.FlatBottomRestraintForce(1000*kilojoules_per_mole, 3.5*angstrom,
                                                            restrained_atom_indices1=[pep_mg_pairs[1][0]],
                                                            restrained_atom_indices2=[pep_mg_pairs[1][1]],
                                                           controlling_parameter_name='pep_mg_4')
    system.addForce(pep_mg_rest3)
    system.addForce(pep_mg_rest4)
''' 
if ligand == 'sep':
    sep_gly_rest1 = omt.forces.FlatBottomRestraintForce(1000*kilojoules_per_mole, 5*angstrom,
                                                            restrained_atom_indices1=[sep_gly_pairs[0][0]],
                                                            restrained_atom_indices2=[sep_gly_pairs[0][1]],
                                                           controlling_parameter_name='sep_gly_1')

    sep_gly_rest2 = omt.forces.FlatBottomRestraintForce(1000*kilojoules_per_mole, 5*angstrom,
                                                            restrained_atom_indices1=[sep_gly_pairs[1][0]],
                                                            restrained_atom_indices2=[sep_gly_pairs[1][1]],
                                                           controlling_parameter_name='sep_gly_2')
    system.addForce(sep_gly_rest1)
    system.addForce(sep_gly_rest2)

    sep_mg_rest3 = omt.forces.FlatBottomRestraintForce(1000*kilojoules_per_mole, 3.5*angstrom,
                                                            restrained_atom_indices1=[sep_mg_pairs[0][0]],
                                                            restrained_atom_indices2=[sep_mg_pairs[0][1]],
                                                           controlling_parameter_name='sep_mg_3')
    sep_mg_rest4 = omt.forces.FlatBottomRestraintForce(1000*kilojoules_per_mole, 3.5*angstrom,
                                                            restrained_atom_indices1=[sep_mg_pairs[1][0]],
                                                            restrained_atom_indices2=[sep_mg_pairs[1][1]],
                                                           controlling_parameter_name='sep_mg_4')
    system.addForce(sep_mg_rest3)
    system.addForce(sep_mg_rest4)
    
######################################## END Distance Restraint Setup #####################################
simulation.context.reinitialize()
simulation.context.setPositions(modeller.positions)
print(f'After adding the dihedral and distance restraints there are {len(system.getForces())} in the system.',flush=True)


# In[ ]:


####### Energy Minimize
simulation.context.reinitialize()
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()
with open(f'{construct}_{ligand}_minimized.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
#######


# In[ ]:


############################ Create Position Restraints for Water Equilibration #########################
# Get all the heavy protein, ligand, and MG atoms
atoms_to_restrain = []
for residue in topology.residues():
    if residue.name not in ['HOH','CL','NA']:
        for atom in residue.atoms():
            if atom.element.symbol != 'H':
                atoms_to_restrain.append(atom.index)


# In[ ]:
# Add position restraints to heavy atoms to allow water to relax around protein
# Create the restraint object, force, and add the particles to it
positions = simulation.context.getState(getPositions=True).getPositions()
reference_coordinates = positions.in_units_of(unit.nanometer)
restraint_weight = 2 * unit.kilocalories_per_mole / unit.angstrom ** 2
restraint_force = CustomExternalForce('K*periodicdistance(x, y, z, x0, y0, z0)^2')
# Add the restraint weight as a global parameter in kcal/mol/nm^2
restraint_force.addGlobalParameter("K", restraint_weight)
# Define the target xyz coords for the restraint as per-atom (per-particle) parameters
restraint_force.addPerParticleParameter("x0")
restraint_force.addPerParticleParameter("y0")
restraint_force.addPerParticleParameter("z0")

for index in range(0, len(positions)):
    if index in atoms_to_restrain:
        xyz = reference_coordinates[index].in_units_of(unit.nanometers) / unit.nanometers
        restraint_force.addParticle(index, xyz)
custom_forces ={}        
custom_forces['positional_restraints'] = system.addForce(restraint_force)


##################################### END Create Position Restraints #################
simulation.context.reinitialize(preserveState=True)
print(f'After adding the position restraints there are {len(system.getForces())} in the system.',flush=True)


# In[ ]:


# def file output directory
eq_out = 'equilibration_output'
if os.path.exists(eq_out):
    pass
else:
    os.makedirs(eq_out)

########## 1ns water equilibration with protein/lig/Mg restrained
simulation.context.reinitialize()
simulation.context.setPositions(positions)
simulation.reporters.append(DCDReporter(f'{eq_out}/{construct}_{ligand}_restrained.dcd', 2500))
simulation.step(500000)
##################

######################## Remove Position Restraints ####################
# After the equilibration, get the positions and then remove the restraints
# set the positions to the post-equilibration state but discard the velocities for slow heating of the protein
system.removeForce(custom_forces['positional_restraints'])
simulation.context.reinitialize(preserveState=True)
######################################################

print(f'After removing the position restraints there are {len(system.getForces())} in the system.',flush=True)

### Save Equilibration Output
positions = simulation.context.getState(getPositions=True).getPositions()
with open(f'{eq_out}/{construct}_{ligand}_water_eq.xml', 'w') as outfile:
    outfile.write(XmlSerializer.serialize(system))
with open(f'{eq_out}/{construct}_{ligand}_water_eq.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, positions, f)
# Don't preserve state info
state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True)
with open(f'{eq_out}/{construct}_{ligand}_water_eq.xml', 'w') as f:
    f.write(XmlSerializer.serialize(state))
###############################

########################### Heat to 310 K over 10 ns ###############################
simulation.reporters.append(DCDReporter(f'{eq_out}/{construct}_{ligand}_heat.dcd', 5000))
simulation.reporters.append(StateDataReporter(f'{eq_out}/heat.log', 5000, step=True,
        potentialEnergy=True, temperature=True, elapsedTime=True, speed=True))
simulation.reporters.append(CheckpointReporter(f'{eq_out}/heat_checkpnt.chk',5000))
# 10ns heating
for i in range(1,311):
    temperature = i*kelvin
    simulation.context.setParameter(MonteCarloBarostat.Temperature(), temperature)
    integrator.setTemperature(temperature)
    # Do you need to simulation.context.reinitialize(preserveState=True) when adjusting the integrator?
    # is integrator.step() the same thing as simulation.step()?
    simulation.context.reinitialize(preserveState=True)
    simulation.step(16000)

# simulation.context.getState(getForces=True) to check if dihedral and distance restraints are on
# Save Everything
positions = simulation.context.getState(getPositions=True).getPositions()
with open(f'{eq_out}/{construct}_{ligand}_equilibrated_system.xml', 'w') as outfile:
    outfile.write(XmlSerializer.serialize(system))
with open(f'{eq_out}/{construct}_{ligand}_md_start.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, positions, f)
# have to use this line to avoid including the global parameter section in the state xml file which will cause the
# restart to error out.
state = simulation.context.getState(getPositions=True, getVelocities=True, getForces=True)
with open(f'{eq_out}/{construct}_{ligand}_md_start.xml', 'w') as f:
    f.write(XmlSerializer.serialize(state))
simulation.saveState('equilibrated_state.xml')
simulation.saveCheckpoint('md_start.chk')
    
###################### End heating - system is ready for run ##########################


# In[ ]:



############## Process out waters and save an xtc of the heating to check the system
u = mda.Universe(f'{construct}_{ligand}_minimized.pdb',f'{eq_out}/{construct}_{ligand}_heat.dcd')
selection = u.select_atoms(f'protein or resname {ligand.upper()} or resname MG')                
with mda.Writer(f'{construct}_{ligand}_heat_check.xtc', selection.n_atoms) as W:
    for ts in u.trajectory:
        W.write(selection)

