#!/bin/bash
#SBATCH -C haswell
#SBATCH --nodes=10
#SBATCH --time=30:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jesteves@umich.edu
#SBATCH --output=log.out
#SBATCH --error=log.err

srun -N 10 -n 640 python run_nbodyKit.py 'power' &

# srun -N 5 -n 320 python run_nbodyKit.py 'xi' &
# srun -N 10 -n 320 python run_nbodyKit.py 'power' &
wait

# srun -N 1 -n 64 python run_nbodyKit.py 'xi'
# srun -N 1 -n 64 python run_nbodyKit.py 'power' &
# wait
