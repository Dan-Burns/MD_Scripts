#!/bin/bash

#SBATCH --time=2:52:00   # walltime limit (HH:MM:SS)
#SBATCH --nodes=1  # number of nodes
#SBATCH --job-name="fto_aKG_off2_process"
#SBATCH --mail-user=dburns@iastate.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module use  /opt/rit/singularity/devel.2022/singularity/modules
module load gromacs/2022.4-cuda11.6.2/cpu-nompi

construct="fto-aKG-off2"

# These steps (often the first two) are sufficient to recenter broken up, multi-chain trajectories.
# The selections that are piped to gmx trjconv are non-solvent
printf "%s\n" "24" "24" | gmx trjconv -f ${construct}_concatenated.xtc -s *.tpr -pbc nojump -center -o step1.xtc -n index.ndx
printf "%s\n" "24" "24" "24" | gmx trjconv -f step1.xtc -s *.tpr -pbc cluster -center -o ${construct}_centered.xtc -n index.ndx
#printf "%s\n" "24" | gmx trjconv -f ${construct}_centered.xtc -o ${construct}_no_water.xtc -n index.ndx
