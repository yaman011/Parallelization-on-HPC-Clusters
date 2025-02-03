#!/bin/bash
#SBATCH --job-name="fdtd_weak_scaling_64proc"
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=09:00
#SBATCH --output=weak_scaling_64proc.out

module load OpenMPI

srun ./fdtd param_3d_64proc.txt
