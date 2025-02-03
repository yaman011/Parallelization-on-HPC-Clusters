#!/bin/bash
#SBATCH --job-name="mpi_weak_1n"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --time=00:05:00
#SBATCH --output=mpi_weak_1n.out

module load OpenMPI
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}


srun ./fdtd param_3d_baseline.txt
