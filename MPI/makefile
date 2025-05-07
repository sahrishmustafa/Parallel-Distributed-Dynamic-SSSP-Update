# Compiler settings
MPICXX = mpicxx
CXX = g++

# Directories
SRC_DIR = src
INC_DIR = include
BIN_DIR = bin

# Paths for METIS
METIS_INCLUDE = /home/kali/local/include
METIS_LIB = /home/kali/local/lib

# Flags
CXXFLAGS = -std=c++17 -I$(INC_DIR)
METIS_FLAGS = -I$(METIS_INCLUDE)
METIS_LDFLAGS = -L$(METIS_LIB) -Wl,-rpath,$(METIS_LIB)
METIS_LIBS = -lmetis -lGKlib

MPI_CXXFLAGS = -I/usr/include/x86_64-linux-gnu/mpich
MPI_LDFLAGS = -L/usr/lib/x86_64-linux-gnu -lmpichcxx -lmpich


# Output targets
PREPROCESS_EXE = $(BIN_DIR)/preprocess
SSSP_EXE = $(BIN_DIR)/sssp

# Object files
PREPROCESS_OBJS = $(SRC_DIR)/preprocess.o $(SRC_DIR)/preprocess_graph.o
SSSP_OBJS = $(SRC_DIR)/distributed_sssp.o $(SRC_DIR)/updation.o

# Default target
all: $(PREPROCESS_EXE) $(SSSP_EXE)

# Create bin directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Preprocessing (non-MPI)
$(PREPROCESS_EXE):$(PREPROCESS_OBJS)
	$(CXX) $(CXXFLAGS) $(METIS_FLAGS) -o $@ $^ $(METIS_LDFLAGS) $(METIS_LIBS)

# Distributed MPI computation
$(SSSP_EXE): $(SSSP_OBJS)
	$(MPICXX) $(CXXFLAGS) $(MPI_CXXFLAGS) -o $@ $^ $(MPI_LDFLAGS)

# Compilation rules
$(SRC_DIR)/preprocess_graph.o: $(SRC_DIR)/preprocess_graph.cpp
	$(CXX) $(CXXFLAGS) $(METIS_FLAGS) -c $< -o $@

$(SRC_DIR)/preprocess.o: $(SRC_DIR)/preprocess.cpp
	$(CXX) $(CXXFLAGS) $(METIS_FLAGS) -c $< -o $@

# Compile MPI source files with mpicxx
$(SRC_DIR)/distributed_sssp.o: $(SRC_DIR)/distributed_sssp.cpp
	$(MPICXX) $(CXXFLAGS) $(MPI_CXXFLAGS) -c $< -o $@

$(SRC_DIR)/updation.o: $(SRC_DIR)/updation.cpp
	$(MPICXX) $(CXXFLAGS) $(MPI_CXXFLAGS) -c $< -o $@


# $(SRC_DIR)/%.o: $(SRC_DIR)/%.cpp
# 	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -f $(SRC_DIR)/*.o $(PREPROCESS_EXE) $(SSSP_EXE)
	rm -f pre-parts/*


# Run Preprocessing
run-preprocess: $(PREPROCESS_EXE)
	@if [ "$(NUM_PARTS)" = "" ]; then \
		echo "Usage: make run-preprocess NUM_PARTS=<number_of_partitions>"; \
		exit 1; \
	fi; \
	./$(PREPROCESS_EXE) dataset/twitch.graph $(NUM_PARTS)


# Run on Cluster with MPI
cluster-run: $(SSSP_EXE)
	@if [ "$(NUM_PARTS)" = "" ]; then \
		echo "Usage: make cluster-run NUM_PARTS=<number_of_partitions>"; \
		exit 1; \
	fi; \
	echo "Running on cluster with $$(($(NUM_PARTS)+1)) MPI processes"; \
	mpiexec -machinefile machinefile -n $$(($(NUM_PARTS)+1)) ./$(SSSP_EXE)

# Sync code to slave, excluding the preprocess binary
sync-slave:
	@if [ "$(SLAVE)" = "" ]; then \
		echo "Usage: make sync-slave SLAVE=user@ip"; \
		exit 1; \
	fi; \
	rsync -avz --exclude=bin/sssp ./ $$SLAVE:~/Desktop/project



#.PHONY: all clean run-preprocess cluster-run


####################################################
####################################################
####################################################
# Profiling targets
.PHONY: profile-time profile-perf profile-strace

# GNU time profiling (install if missing)
profile-time: $(SSSP_EXE)
	@if [ "$(NUM_PARTS)" = "" ]; then \
		echo "Usage: make profile-time NUM_PARTS=<number_of_partitions>"; \
		exit 1; \
	fi; \
	if ! command -v gtime >/dev/null 2>&1; then \
		echo "Installing GNU time..."; \
		sudo apt install time; \
	fi; \
	/usr/bin/time -v mpiexec -machinefile machinefile -n $$(($(NUM_PARTS)+1)) ./$(SSSP_EXE)

# perf profiling
profile-perf: $(SSSP_EXE)
	@if [ "$(NUM_PARTS)" = "" ]; then \
		echo "Usage: make profile-perf NUM_PARTS=<number_of_partitions>"; \
		exit 1; \
	fi; \
	sudo perf stat -d mpiexec -machinefile machinefile -n $$(($(NUM_PARTS)+1)) ./$(SSSP_EXE)

# strace profiling
profile-strace: $(SSSP_EXE)
	@if [ "$(NUM_PARTS)" = "" ]; then \
		echo "Usage: make profile-strace NUM_PARTS=<number_of_partitions>"; \
		exit 1; \
	fi; \
	strace -c -f mpiexec -machinefile machinefile -n $$(($(NUM_PARTS)+1)) ./$(SSSP_EXE)