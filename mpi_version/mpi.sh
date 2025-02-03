#!/bin/bash
#SBATCH --job-name="fdtd mpi"
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00
#SBATCH --output=mpi.out

module load OpenMPI

cd simple3d
srun ../fdtd param_3d.txt