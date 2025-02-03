# High-Performance FDTD Solver with MPI and OpenMP

## 📌 Project Overview
This project implements a **parallelized 3D Finite-Difference Time-Domain (FDTD) solver** using **MPI and OpenMP** for distributed and shared memory architectures. It efficiently models wave propagation by partitioning the computational grid across multiple processes and threads.

## 🚀 Key Features
- **MPI-Based Domain Decomposition**: Uses a 3D Cartesian grid with **MPI_Dims_create** for parallel subdomain partitioning.
- **Hybrid Parallelism (MPI + OpenMP)**: MPI for inter-node communication, OpenMP for intra-node threading.
- **Strong & Weak Scaling Analysis**: Performance benchmarking using varying numbers of processors.
- **Optimized Communication**: Implements **non-blocking MPI_Isend/MPI_Irecv** for efficient boundary data exchange.
- **Score-P Profiling**: Performance profiling to analyze execution time, computational load, and communication overhead.

## 🛠️ Installation & Compilation
### **Prerequisites**
- MPI (e.g., OpenMPI, MPICH)
- OpenMP support
- C compiler (e.g., GCC)
- Python (for visualization and analysis)
- LIKWID (optional, for hardware performance monitoring)

### **Compilation**
```sh
mpicc -fopenmp -o fdtd_solver fdtd.c -lm
