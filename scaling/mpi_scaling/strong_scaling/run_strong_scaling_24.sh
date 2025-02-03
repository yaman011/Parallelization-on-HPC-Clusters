#!/bin/bash
#SBATCH --job-name="fdtd_mpi_strong_scaling_24"
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00
#SBATCH --output=mpi_strong_24.out

module load OpenMPI

cd simple3d
srun ../fdtd param_3d.txt
