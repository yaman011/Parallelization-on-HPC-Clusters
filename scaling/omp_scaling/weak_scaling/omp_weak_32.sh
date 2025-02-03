#!/bin/bash
#SBATCH --job-name="fdtd omp weak 32"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=1024
#SBATCH --time=00:20:00
#SBATCH --output=omp_weak_32.out

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load releases/2021b
module load GCC/11.2.0


./fdtd param_3d_32cpu.txt
