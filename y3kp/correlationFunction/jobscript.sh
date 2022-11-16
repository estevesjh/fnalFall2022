#!/bin/bash -l
#SBATCH --job-name=corr_func_jk
#SBATCH --nodes=1
#SBATCH --qos=debug
#SBATCH --time=05:00
#SBATCH --constraint=haswell
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jesteves@umich.edu
#SBATCH --output=log.out
#SBATCH --error=log.err

python run.py 'rmy3_3d_smallScales' --is_3d 1 --nPatches 20 --nCores 31
python run.py 'rmy3_smallScales' --is_3d 0 --nPatches 20 --nCores 31

# python run.py 'rmy3_3d_smallScales' --nPatches 30 --nCores 31