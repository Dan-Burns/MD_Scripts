#!/usr/bin/env python
# coding: utf-8

from openmm.app import *
from openmm import *
from openff.toolkit.topology import *
from openmm.unit import *
import openmmtools as omt
from openmmtools.integrators import *

import parmed
import MDAnalysis as mda

import os
import numpy as np
import pandas as pd
import re
from sys import stdout



######## Relative location of input files ############
path = '..'
#######################################
############## System Setup ###########################
name = 'fto-m1a-in2'
base = 'D1A'

prmtop = AmberPrmtopFile(f'{path}/{name}-sol.prmtop')
inpcrd = AmberInpcrdFile(f'{path}/{name}-sol.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
        constraints=HBonds)
temperature = 303.15*kelvin
integrator = LangevinMiddleIntegrator(temperature, 2/picosecond, 0.002*picoseconds)
pressure = 1 * bar
system.addForce(MonteCarloBarostat(pressure, temperature))
simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
############### End Setup #############################

##### Save a gromacs topology for future trjconv use - Use a no-constraints version of system to avoid parmed error
parmed_system = prmtop.createSystem(rigidWater=False)
pmd_structure = parmed.openmm.load_topology(simulation.topology, system=parmed_system, xyz=inpcrd.positions)

if os.path.exists('gromacs'):
    pass
else:
    os.makedirs('gromacs')

pmd_structure.save(f"gromacs/{name}_SYSTEM.top", overwrite=True)
pmd_structure.save(f"gromacs/{name}_SYSTEM.gro", overwrite=True)
########### End gromacs ######################################

'''Please restrain the AKG, Fe(II), the methylated base, and the side chains
of residues His 231, Asp 233, and His 307 to the starting position'''

####################### Get atom indices for position restraints ####################
######## Edit if statement for systems that don't have base/ligand/iron ###############
# +37 will return res.INDEX matching Vincenzo's residue numbers.
topology = simulation.topology
restraint_indices = []
# confirm you got all the residues you need to restrain
res_count = 0 # akg, iron, base, and 3 residues
for res in topology.residues():
    if (res.name == 'HIS' and int(res.index)+37 == 231) \
    or (res.name == 'HIS' and int(res.index)+37 == 307)\
    or (res.name == 'ASP' and int(res.index)+37 == 233)\
    or res.name == base\
    or res.name == 'AKG'\
    or res.name == 'FE2':
        res_count += 1
        print(f'Getting atom indices for {res.name} {int(res.index)+37} position restraints',flush=True)
        restraint_indices.extend([atom.index for atom in res.atoms() 
                                  if atom.element.symbol != 'H'])

assert res_count == 6

# Add position restraints to heavy atoms to allow water to relax around protein
# Create the restraint object, force, and add the particles to it
positions = simulation.context.getState(getPositions=True).getPositions()
reference_coordinates = positions.in_units_of(unit.nanometer)
restraint_weight = 5 * unit.kilocalories_per_mole / unit.angstrom ** 2
restraint_force = CustomExternalForce('K*periodicdistance(x, y, z, x0, y0, z0)^2')
# Add the restraint weight as a global parameter in kcal/mol/nm^2
restraint_force.addGlobalParameter("K", restraint_weight)
# Define the target xyz coords for the restraint as per-atom (per-particle) parameters
restraint_force.addPerParticleParameter("x0")
restraint_force.addPerParticleParameter("y0")
restraint_force.addPerParticleParameter("z0")
for index in range(0, len(positions)):
    if index in restraint_indices:
        xyz = reference_coordinates[index].in_units_of(unit.nanometers) / unit.nanometers
        restraint_force.addParticle(index, xyz)
custom_forces ={}
custom_forces['positional_restraints'] = system.addForce(restraint_force)


##################################### END Create Position Restraints #################
########## Setup Iron to AKG distance restraints - get the atom indices ################
akg_indices = []
fe2_index = []
for res in topology.residues():
    if res.name == 'AKG':
        for atom in res.atoms():
            if atom.name == 'O2' or atom.name == 'O5':
                print(f'adding {res.name} atom {atom.name} to the distance restraint indices',flush=True)
                akg_indices.append(atom.index)
    if res.name == 'FE2':
        for atom in res.atoms():
            print(f'adding {res.name} atom to the distance restraint indices',flush=True)
            fe2_index.append(atom.index)
            ### AKG atom names for restraints are O2 and O5
# 100Kj/nm^2/mol restraint
iron_akg_restraint = omt.forces.HarmonicRestraintForce(100000, 
                                                    restrained_atom_indices1=akg_indices,
                                                    restrained_atom_indices2=fe2_index,
                                                   )

system.addForce(iron_akg_restraint)
################# End distance restraint setup ################################


####### Energy Minimize #####################
simulation.context.reinitialize(preserveState=True)
simulation.context.setPositions(inpcrd.positions)
simulation.minimizeEnergy()
with open(f'{name}_minimized.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
####### End Minimization ##################

############################ Begin Equilibration ###############################
# def file output directory
eq_out = 'equilibration_output'
if os.path.exists(eq_out):
    pass
else:
    os.makedirs(eq_out)

########################### Heat to 303 K over 10 ns ###############################

simulation.reporters.append(DCDReporter(f'{eq_out}/{name}_heat.dcd', 5000))
simulation.reporters.append(StateDataReporter(f'{eq_out}/{name}_heat.log', 5000, step=True,
        potentialEnergy=True, temperature=True, elapsedTime=True, speed=True))
simulation.reporters.append(CheckpointReporter(f'{eq_out}/heat_checkpnt.chk',50000))

# 10ns heating
steps_per_ps = 500 # 500 * .002ps = 1ps
ps_to_run = 10000 # 10ns
steps = steps_per_ps * ps_to_run
# raise by 1K every n steps/303
for i in range(1,304):
    temperature = i*kelvin
    simulation.context.setParameter(MonteCarloBarostat.Temperature(), temperature)
    integrator.setTemperature(temperature)
    # Do you need to simulation.context.reinitialize(preserveState=True) when adjusting the integrator?
    # is integrator.step() the same thing as simulation.step()?
    simulation.context.reinitialize(preserveState=True)
    simulation.step(int(steps/303))
# clear the reporters for next step
simulation.reporters.clear()
 ################# End Heating ########################   
# Equilibrate for 100 ns
simulation.reporters.append(DCDReporter(f'{eq_out}/{name}_eq.dcd', 5000))
simulation.reporters.append(StateDataReporter(f'{eq_out}/{name}_eq.log', 5000, step=True,
        potentialEnergy=True, temperature=True, elapsedTime=True, speed=True))
simulation.reporters.append(CheckpointReporter(f'{eq_out}/eq_checkpnt.chk',50000))
simulation.step(steps_per_ps*150000)
################## End Equilibration ###############################


######################## Remove Position Restraints ####################
# After the equilibration, get the positions and then remove the restraints
# set the positions to the post-equilibration state but discard the velocities for slow heating of the protein
system.removeForce(custom_forces['positional_restraints'])
simulation.context.reinitialize(preserveState=True)
######################################################


####################### Save Everything  ################################
positions = simulation.context.getState(getPositions=True).getPositions()
with open(f'{eq_out}/{name}_equilibrated_system.xml', 'w') as outfile:
    outfile.write(XmlSerializer.serialize(system))
with open(f'{eq_out}/{name}_md_start.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, positions, f)

simulation.saveState(f'{eq_out}/{name}_equilibrated_state.xml')
simulation.saveCheckpoint(f'{eq_out}/{name}_md_start.chk')

############## Process out waters and save an xtc of the heating to check the system
u = mda.Universe(f'{name}_minimized.pdb',f'{eq_out}/{name}_eq.dcd')
selection = u.select_atoms(f'not (resname HOH or resname Na+')
with mda.Writer(f'{eq_out}/{name}_eq_check.xtc', selection.n_atoms) as W:
    for ts in u.trajectory:
        W.write(selection)

