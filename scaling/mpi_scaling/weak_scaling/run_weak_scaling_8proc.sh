#!/bin/bash
#SBATCH --job-name="fdtd_weak_scaling_8proc"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00
#SBATCH --output=weak_scaling_8proc.out

module load OpenMPI


srun ./fdtd param_3d_8proc.txt
