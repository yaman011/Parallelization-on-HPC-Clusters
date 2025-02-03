#!/bin/bash
#SBATCH --job-name="fdtd_weak_scaling_16proc"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00
#SBATCH --output=weak_scaling_16proc.out

module load OpenMPI

srun ./fdtd param_3d_16proc.txt
