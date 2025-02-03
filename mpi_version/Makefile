# Makefile for compiling and running FDTD with MPI support using OpenMPI

# Compiler
MPICC = mpicc

# Compiler flags
CFLAGS = -g

# Executable name
TARGET = fdtd

# Source files
SRCS = fdtd.c

# Base directory (where this Makefile is located)
BASE_DIR = $(shell pwd)

# Directory for the parameter file
PARAM_DIR = simple3d

# Number of processes to run with mpirun
NP = 27

# Load the necessary module
LOAD_MODULE = module load OpenMPI

# Run command
RUN_CMD = mpirun -np $(NP) $(BASE_DIR)/$(TARGET) param_3d.txt

# Default target to build and run the application
all: $(TARGET) run

# Compile the source code
$(TARGET): $(SRCS)
	$(LOAD_MODULE) && $(MPICC) $(CFLAGS) -o $(TARGET) $(SRCS) -lm

# Run the application
run:
	cd $(PARAM_DIR) && $(RUN_CMD)

# Rule to clean up the build
clean:
	rm -f $(TARGET)

.PHONY: all clean
