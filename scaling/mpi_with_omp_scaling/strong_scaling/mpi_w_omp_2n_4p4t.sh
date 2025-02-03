#!/bin/bash
#SBATCH --job-name="mpi_w_omp_2n_4p4t"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=4
#SBATCH --time=00:05:00
#SBATCH --output=mpi_w_omp_2n_4p4t.out

module load OpenMPI
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd simple3d
srun ../fdtd param_3d.txt
