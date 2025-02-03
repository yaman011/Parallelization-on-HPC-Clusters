#!/bin/bash
#SBATCH --job-name="fdtd_mpi_strong_scaling_48"
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=05:00
#SBATCH --output=mpi_strong_48.out

module load OpenMPI

cd simple3d
srun ../fdtd param_3d.txt
