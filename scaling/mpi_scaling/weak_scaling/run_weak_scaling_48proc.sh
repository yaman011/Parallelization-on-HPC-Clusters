#!/bin/bash
#SBATCH --job-name="fdtd_weak_scaling_48proc"
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00
#SBATCH --output=weak_scaling_48proc.out

module load OpenMPI

srun ./fdtd param_3d_48proc.txt
