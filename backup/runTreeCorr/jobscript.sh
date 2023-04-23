#!/bin/bash -l
#SBATCH --job-name=corr_func_jk
#SBATCH --nodes=2
#SBATCH --qos=debug
#SBATCH --time=10:00
#SBATCH --constraint=haswell
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jesteves@umich.edu
#SBATCH --output=log.out
#SBATCH --error=log.err

srun -N 2 -n 128 python run.py 'mock_test' --nPatches 20 --nCores 128