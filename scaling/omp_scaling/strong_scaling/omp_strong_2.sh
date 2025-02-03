#!/bin/bash
#SBATCH --job-name="fdtd omp"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=1024
#SBATCH --time=00:20:00
#SBATCH --output=omp_strong_2.out

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load releases/2021b
module load GCC/11.2.0

cd simple3d
../fdtd param_3d.txt
