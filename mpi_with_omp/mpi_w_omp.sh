#!/bin/bash
#SBATCH --job-name="mpi_w_omp"
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --time=00:05:00
#SBATCH --output=mpi_w_omp.out

module load OpenMPI
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd simple3d
srun ../fdtd param_3d.txt