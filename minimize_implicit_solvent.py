from openmm.app import *
from openmm import *
from openff.toolkit.topology import *
from openmm.unit import *
import openmmtools as omt
from openmmtools.integrators import *
from openmm import app
import parmed
import MDAnalysis as mda

import os
import numpy as np
import pandas as pd
import re
from sys import stdout





# load the pdb, probably could have skipped Modeller() and just used pdb.topology and pdb.positions
pdb = PDBFile(f'{dir_70C}/{i}_ready.pdb')
modeller = Modeller(pdb.topology, pdb.positions)

# pick a temp
temperature = 300*kelvin
forcefield = ForceField('amber14-all.xml')
integrator = LangevinMiddleIntegrator(temperature, 2/picosecond, 0.002*picoseconds)
# using implicit solvent, no cutoff is used
system = forcefield.createSystem(modeller.topology,nonbondedMethod=app.NoCutoff,
                               constraints=app.HBonds, implicitSolvent=app.OBC2,
                               implicitSolventSaltConc=0.1*moles/liter,
)

## CREATE THE SIMULATION
temperature = 303.15*kelvin
integrator = LangevinMiddleIntegrator(temperature, 2/picosecond, 0.002*picoseconds)
# No pressue with implicit solvent
#pressure = 1 * bar
#system.addForce(MonteCarloBarostat(pressure, temperature))
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

### GET RESTRAINED ATOMS
topology = simulation.topology
# atoms to restrain
restraint_indices = []
# get the backbone CA, N, C
for atom in topology.atoms():
    if atom.name == 'CA' or \
       atom.name == 'N' or \
       atom.name == 'C':
        restraint_indices.append(atom.index)



# Add position restraints to heavy atoms to allow water to relax around protein
# Create the restraint object, force, and add the particles to it
positions = simulation.context.getState(getPositions=True).getPositions()
reference_coordinates = positions.in_units_of(unit.nanometer)
restraint_weight = 500 * unit.kilocalories_per_mole / unit.angstrom ** 2
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

####### Energy Minimize #####################
simulation.context.reinitialize()
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()
with open(f'{dir_70C}/{i}_ready_relaxed.pdb', 'w') as f:
    PDBFile.writeFile(simulation.topology, simulation.context.getState(getPositions=True).getPositions(), f)
####### End Minimization ##################
