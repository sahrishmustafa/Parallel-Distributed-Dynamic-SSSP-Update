# Compiler settings
MPICXX = mpicxx

# Add this to detect METIS configuration and enable OpenMP
CXXFLAGS = -std=c++17 -fopenmp -I/usr/include -Iinclude
LDFLAGS = -L/usr/lib/x86_64-linux-gnu -fopenmp
LIBS = -lmetis

# Source files
SRC_DIR = src
SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ = $(SRC:.cpp=.o)

# Target
TARGET = partition

# Default target
all: $(TARGET)

# Main build rule
$(TARGET): $(OBJ)
	$(MPICXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) $(LIBS)

# Compile source files
$(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(MPICXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(OBJ) $(TARGET) partition_output.part

# Run with mpiexec - calculates processes as partitions + 1
run: $(TARGET)
	@if [ "$(NUM_PARTS)" = "" ]; then \
		echo "Usage: make run NUM_PARTS=<number_of_partitions>"; \
		echo "Example: make run NUM_PARTS=4"; \
		exit 1; \
	fi; \
	echo "Running with $$(($(NUM_PARTS)+1)) MPI processes (${NUM_PARTS} partitions + 1 master)"; \
	mpiexec -n $$(($(NUM_PARTS)+1)) ./$(TARGET) ../dataset/email-enron/Email-Enron.graph $(NUM_PARTS) || \
	(echo "MPI execution failed with status $$?"; exit 1)

.PHONY: all clean run