#!/bin/bash -l
#
#SBATCH --job-name="LIKWID OpenMP"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=15:00
#SBATCH --output=likwid_omp.out

module load likwid

export OMP_PROC_BIND=true
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

likwid-perfctr -g CACHE ./fdtd param_3d.txt

