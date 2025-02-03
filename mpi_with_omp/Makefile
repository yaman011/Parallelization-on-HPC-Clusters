# Compiler
MPICC = mpicc

# Compiler flags for OpenMP support
CFLAGS = -g -fopenmp

# Target executable
TARGET = fdtd

# Source file
SRC = fdtd.c

# Libraries
LIBS = -lm

# Number of processes to run with mpirun
NP = 27

# Base directory (where this Makefile is located)
BASE_DIR = $(shell pwd)

# Directory for the parameter file
PARAM_DIR = simple3d

# Load the necessary module
LOAD_MODULE = module load OpenMPI

# Run command
RUN_CMD = mpirun -np $(NP) $(BASE_DIR)/$(TARGET) param_3d.txt

# Default target to build and run the program
all: $(TARGET) run

# Rule to compile the program with OpenMP support
$(TARGET): $(SRC)
	$(LOAD_MODULE) && $(MPICC) $(CFLAGS) -o $(TARGET) $(SRC) $(LIBS)

# Rule to run the program
run:
	$(LOAD_MODULE) && cd $(PARAM_DIR) && $(RUN_CMD)

# Rule to clean up the build
clean:
	rm -f $(TARGET)

.PHONY: all clean
