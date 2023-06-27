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


run = 6

name = 'fto-m1a-in2'

# def file output directory
run_out = 'run_output'
if os.path.exists(run_out):
    pass
else:
    os.makedirs(run_out)
    
# Check that this run hasn't already been done.
if run == 1:
    if os.path.exists(f'{run_out}/{name}_run_{run}.dcd'):
        while os.path.exists(f'{run_out}/{name}_run_{run}.dcd'):
            run +=1
        print(f'Run {run} is the next run to conduct. \n Running simulation {run}.',flush=True)
    else:
        pass
# If the trajectory file corresponding to this run number exists, check which run hasn't been done yet and start
# that simulation. This is probably all you need and don't even have to add a new run number.
elif os.path.exists(f'{run_out}/{name}_run_{run}.dcd'):
    print(f'Run {run} already executed')
    run += 1
    while os.path.exists(f'{run_out}/{name}_run_{run}.dcd'):
        run +=1
    print(f'Run {run} is the next run to conduct. \n Running simulation {run}.',flush=True)
# If the preceding file number does not exist, see how far back you have to go and exit the script without running any
# simulation
else:
    if os.path.exists(f'{run_out}/{name}_run_{run - 1}.dcd'):
        pass
    else:
        run -= 1
        while os.path.exists(f'{run_out}/{name}_run_{run}.dcd') == False:
            run -= 1
            if run < 0:
                print(f'Cant find a completed run. \n Not running any simulation.',flush=True)
                quit()
        print(f'Run {run} was the last run conducted. \n Not running any simulation.',flush=True)
        quit()
temperature = 303        
integrator = LangevinMiddleIntegrator(temperature, 2/picosecond, 0.002*picoseconds)
steps_per_ps = 500 # 500 * .002ps = 1ps

if run == 1:
    # Load the system
    with open(f'equilibration_output/{name}_equilibrated_system.xml') as infile:
        system = XmlSerializer.deserialize(infile.read())
    system_pdb = PDBFile(f'equilibration_output/{name}_md_start.pdb')

    simulation = Simulation(system_pdb.topology,system, integrator)

    # heat
    simulation.context.setPositions(system_pdb.positions)
    simulation.context.reinitialize(preserveState=True)
    ############################ Begin heat ###############################
    # def file output directory
    heat_out = 'heat_output'
    if os.path.exists(heat_out):
        pass
    else:
        os.makedirs(heat_out)

    ########################### Heat to 303 K over 10 ns ###############################
    simulation.reporters.append(DCDReporter(f'{heat_out}/{name}_heat.dcd', 5000))
    simulation.reporters.append(StateDataReporter(f'{heat_out}/{name}_heat.log', 5000, step=True,
            potentialEnergy=True, temperature=True, elapsedTime=True, speed=True))
    simulation.reporters.append(CheckpointReporter(f'{heat_out}/heat_checkpnt.chk',50000))

    # 10ns heating
    steps_per_ps = 500 # 500 * .002ps = 1ps
    # changed because this system ended up with base popping out early on
    ps_to_run = 100000 # 50ns
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
    positions = simulation.context.getState(getPositions=True).getPositions()
    # clear the reporters
    simulation.reporters.clear()
    with open(f'{heat_out}/{name}_heat.pdb', 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)
     ################# End Heating ########################   
    ############ new reporters
    simulation.reporters.append(DCDReporter(f'{run_out}/{name}_run_{run}.dcd', 5000))
    simulation.reporters.append(StateDataReporter(f'{run_out}/{name}_run_{run}.log', 5000, step=True,
            potentialEnergy=True, temperature=True, elapsedTime=True, speed=True))
    simulation.reporters.append(CheckpointReporter(f'{run_out}/run_checkpoint_{run}.chk',50000))
    ########### run
    simulation.step(steps_per_ps*250000)
    ########## Save
    positions = simulation.context.getState(getPositions=True).getPositions()
    
    with open(f'{run_out}/{name}_system_{run}.xml', 'w') as outfile:
        outfile.write(XmlSerializer.serialize(system))
    with open(f'{run_out}/{name}_{run}.pdb', 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)
    simulation.saveState(f'{run_out}/{name}_state_{run}.xml')
    simulation.saveCheckpoint(f'{run_out}/{name}_{run}.chk')
    
    
    
# Otherwise load the latest state.
else:
    # Load the system
    with open(f'{run_out}/{name}_system_{run-1}.xml') as infile:
        system = XmlSerializer.deserialize(infile.read())
    # print the forces
    forces = system.getForces()
    print("The following forces are in the system:")
    for force in forces:
        print(force, flush=True)

    system_pdb = PDBFile(f'{run_out}/{name}_{run-1}.pdb')

    simulation = Simulation(system_pdb.topology,system, integrator)
    

    # Load the state
    simulation.loadCheckpoint(f'{run_out}/{name}_{run-1}.chk')

    # Set up the reporters
    simulation.reporters.append(DCDReporter(f'{run_out}/{name}_run_{run}.dcd', 5000))
    simulation.reporters.append(StateDataReporter(f'{run_out}/{name}_run_{run}.log', 5000, step=True,
            potentialEnergy=True, temperature=True, elapsedTime=True, speed=True))
    simulation.reporters.append(CheckpointReporter(f'{run_out}/run_checkpoint_{run}.chk',50000))
    ### Run steps
    simulation.step(steps_per_ps*500000)
    ### Save Everything 
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open(f'{run_out}/{name}_system_{run}.xml', 'w') as outfile:
        outfile.write(XmlSerializer.serialize(system))
    with open(f'{run_out}/{name}_{run}.pdb', 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)
    simulation.saveState(f'{run_out}/{name}_state_{run}.xml')
    simulation.saveCheckpoint(f'{run_out}/{name}_{run}.chk')


############## Process out waters and save an xtc of the heating to check the system
if run ==1:
    pdb = f'{heat_out}/{name}_heat.pdb'
else:
    pdb = f'{run_out}/{name}_{run-1}.pdb'

u = mda.Universe(pdb,f'{run_out}/{name}_run_{run}.dcd')
selection = u.select_atoms(f'not (resname HOH or resname Na+)')
with mda.Writer(f'{run_out}/{name}_run_{run}.xtc', selection.n_atoms) as W:
    for ts in u.trajectory:
        W.write(selection)

