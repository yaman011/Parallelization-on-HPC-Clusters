#!/bin/bash
#SBATCH --job-name="mpi_w_omp_2n_8p2t"
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=2
#SBATCH --time=00:05:00
#SBATCH --output=mpi_w_omp_2n_8p2t.out

module load OpenMPI
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

cd simple3d
srun ../fdtd param_3d.txt
