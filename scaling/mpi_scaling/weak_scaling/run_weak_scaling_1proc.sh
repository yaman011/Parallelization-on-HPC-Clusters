#!/bin/bash
#SBATCH --job-name="fdtd_weak_scaling_1proc"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00
#SBATCH --output=weak_scaling_1proc.out

module load OpenMPI


srun ./fdtd param_3d_1proc.txt
