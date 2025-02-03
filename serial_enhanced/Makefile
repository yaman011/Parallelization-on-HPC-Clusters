# Compiler
CC = gcc

# Compiler flags
CFLAGS = -O3

# Target executable
TARGET = fdtd

# Source file
SRC = fdtd.c

# Libraries
LIBS = -lm

# Base directory (where this Makefile is located)
BASE_DIR = $(shell pwd)

# Input directory
INPUT_DIR = example_inputs/simple3d

# Run command
RUN_CMD = $(BASE_DIR)/$(TARGET) $(BASE_DIR)/$(INPUT_DIR)/param_3d.txt

# Default target to build and run the program
all: $(TARGET) run

# Rule to compile the program
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LIBS)

# Rule to run the program
run:
	cd $(INPUT_DIR) && $(RUN_CMD)

# Rule to clean up the build
clean:
	rm -f $(TARGET)

.PHONY: all clean
