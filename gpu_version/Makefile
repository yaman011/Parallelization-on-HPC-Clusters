# Makefile for compiling and running FDTD with GPU support using Clang

# Compiler
CC = clang

# Compiler flags for GPU
CFLAGS = -O3 -fopenmp -fopenmp-targets=nvptx64-nvidia-cuda

# Executable name
TARGET = fdtd

# Source files
SRCS = fdtd.c

# Base directory (where this Makefile is located)
BASE_DIR = $(shell pwd)

# Directory for the parameter file
PARAM_DIR = simple3d

# Load the necessary module
LOAD_MODULE = module load Clang/16.0.6-GCCcore-11.3.0-CUDA-11.7.0

# Library path (modify this to the correct path if needed)
LIBRARY_PATH = /path/to/your/libs

# Run command
RUN_CMD = LD_LIBRARY_PATH=$(LIBRARY_PATH):$$LD_LIBRARY_PATH $(BASE_DIR)/$(TARGET) $(BASE_DIR)/$(PARAM_DIR)/param_3d.txt

# Default target to build and run the application
all: $(TARGET) run

# Compile the source code
$(TARGET): $(SRCS)
	$(LOAD_MODULE) && $(CC) $(CFLAGS) -o $(TARGET) $(SRCS) -lm

# Run the application
run:
	$(LOAD_MODULE) && cd $(PARAM_DIR) && $(RUN_CMD)

# Rule to clean up the build
clean:
	rm -f $(TARGET)

.PHONY: all clean
