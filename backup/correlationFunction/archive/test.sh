#!/bin/bash -l
#SBATCH --job-name=test
#SBATCH --qos=debug
#SBATCH --nodes=1
#SBATCH --constraint=knl
#SBATCH --time=2
#SBATCH --array=0-2

#SBATCH --mail-type=ALL
#SBATCH --mail-user=jesteves@umich.edu
#SBATCH --output=%a.out

echo $SLURM_ARRAY_TASK_ID
echo $(($SLURM_ARRAY_TASK_ID+9))