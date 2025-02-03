# Makefile for compiling and running FDTD with OpenMP

# Compiler
CC = gcc

# Compiler flags
CFLAGS = -fopenmp -g

# Executable name
TARGET = fdtd

# Source files
SRCS = fdtd.c

# Directory for the parameter file
PARAM_DIR = simple3d

# Number of threads
NUM_THREADS = 27

# Base directory (where this Makefile is located)
BASE_DIR = $(shell pwd)

# Run command
RUN_CMD = export OMP_NUM_THREADS=$(NUM_THREADS) && ./$(TARGET) param_3d.txt

# Default target to build and run the application
all: $(TARGET) run

# Compile the source code
$(TARGET): $(SRCS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) -lm

# Run the application
run:
	cd $(PARAM_DIR) && $(RUN_CMD)

# Clean up the compiled files
clean:
	rm -f $(TARGET)
