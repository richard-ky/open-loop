#!/bin/bash
#
#SBATCH --job-name=open-loop
#SBATCH --partition=compute
#SBATCH --output=output.txt
#
#SBATCH --cpus-per-task=4
#SBATCH --time=2-0
#SBATCH --array=25-50:5
#SBATCH --mem-per-cpu=1000
#SBATCH --mail-user=richard.ky@sjsu.edu
#SBATCH --mail-type=BEGIN,END

module load matlab

matlab -nojvm -batch "main_sim_stim(1,15,$SLURM_ARRAY_TASK_ID)"
