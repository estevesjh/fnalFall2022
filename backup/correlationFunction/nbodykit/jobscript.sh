#!/bin/bash
#SBATCH -C haswell
#SBATCH --nodes=12
#SBATCH --time=30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jesteves@umich.edu
#SBATCH --output=log.out
#SBATCH --error=log.err
srun -N 4 -n 256 python run_nbodyKit.py 0 'xi' &
srun -N 4 -n 256 python run_nbodyKit.py 1 'xi' &
srun -N 4 -n 256 python run_nbodyKit.py 2 'xi' &
wait

# srun -N 2 -n 96 python run_nbodyKit.py 0 'power' &
# srun -N 2 -n 96 python run_nbodyKit.py 1 'power' &
# srun -N 2 -n 96 python run_nbodyKit.py 2 'power' &
# wait