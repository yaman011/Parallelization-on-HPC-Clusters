#!/bin/bash
#SBATCH --job-name="fdtd_mpi_strong_scaling_1"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00
#SBATCH --output=mpi_strong_1.out

module load OpenMPI

cd simple3d
srun ../fdtd param_3d.txt
