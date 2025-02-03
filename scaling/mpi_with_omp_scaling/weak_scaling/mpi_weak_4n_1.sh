#!/bin/bash
#SBATCH --job-name="mpi_weak_4n_1"
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --time=00:10:00
#SBATCH --output=mpi_weak_4n_1.out

module load OpenMPI
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}


srun ./fdtd param_3d_4nodes_1.txt
