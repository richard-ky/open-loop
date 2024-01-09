#!/bin/bash
#
#SBATCH --job-name=open-loop
#SBATCH --partition=compute
#SBATCH --output=output.txt
#
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-user=richard.ky@sjsu.edu
#SBATCH --mail-type=BEGIN,END

module load matlab

matlab -nodisplay -nosplash -nodesktop -r "run('main_sim_stim.m')", exit
