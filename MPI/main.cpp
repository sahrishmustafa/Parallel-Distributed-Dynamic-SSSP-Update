#include "graph_partition.h"
#include <iostream>

// Updated main.cpp
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    GraphPartitionData partition_data;
    std::vector<idx_t> xadj, adjncy, partitions;
    
    if (world_rank == 0) {
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " <graph_file.graph> <num_partitions>" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        int num_parts = std::stoi(argv[2]);
        if (num_parts != world_size - 1) {
            std::cerr << "Error: Number of partitions must match number of workers" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        if (!GraphPartition::loadGraph(argv[1], xadj, adjncy)) {
            return 1;
        }
        GraphPartition::partitionGraph(xadj, adjncy, num_parts, partitions);
    }

    // Broadcast partition mapping to all processes
    if (world_rank == 0) {
        idx_t size = partitions.size();
        MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(partitions.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
    } else {
        idx_t size;
        MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
        partitions.resize(size);
        MPI_Bcast(partitions.data(), size, MPI_INT, 0, MPI_COMM_WORLD);
    }

    GraphPartition::distributePartitions(argc, argv, xadj, adjncy, partitions, partition_data);

    MPI_Barrier(MPI_COMM_WORLD);

    // if (world_rank == 0) {
    //     int source_global = 1; // Set source node
    //     computeDistributedSSSP(source_global, partition_data);
    // } else {
    //     computeDistributedSSSP(-1, partition_data);
    // }
    
    MPI_Finalize();
    return 0;
}