#include "graph_partition.h"
#include <algorithm>
#include <iostream>
#include <fstream>

// Updated main.cpp
int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    GraphPartitionData partition_data;
    std::vector<idx_t> xadj, adjncy, partitions;
    int total_edges=0;
    
    if (world_rank == 0) {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " <graph_file.graph> <num_partitions> <updates_file>" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        int num_parts = std::stoi(argv[2]);
        if (num_parts != world_size - 1) {
            std::cerr << "Error: Number of partitions must match number of workers" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        if (!GraphPartition::loadGraph(argv[1], xadj, adjncy,total_edges)) {
            return 1;
        }
        GraphPartition::partitionGraph(xadj, adjncy, num_parts, partitions, total_edges);
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
    int source_global = 1;
    
    GraphPartition::computeDistributedSSSP(source_global, partition_data);

    MPI_Barrier(MPI_COMM_WORLD);
    // int flag;
    // MPI_Status status;
    // do {
    //     MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
    //     if (flag) {
    //         int dummy[2];
    //         MPI_Recv(dummy, 2, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //         std::cout << "[Rank " << world_rank << "] Drained late ghost message from rank " << status.MPI_SOURCE << "\n";
    //     }
    // } while (flag);

    MPI_Barrier(MPI_COMM_WORLD); // safe to finalize now

    // NEW: Perform update propagation via distributeEdgeUpdates  
    std::string update_file = "deletions.txt";
    std::vector<EdgeUpdate> local_updates;
    GraphPartition::distributeEdgeUpdates(update_file, local_updates, partition_data);

    MPI_Barrier(MPI_COMM_WORLD);

    update_file = "insertions.txt";
    std::vector<EdgeUpdate> local_updates_in;
    GraphPartition::distributeEdgeUpdates_Insertions(update_file, local_updates_in, partition_data);

    // Final SSSP output
    MPI_Barrier(MPI_COMM_WORLD);

    GraphPartition::gatherAndWriteSSSPToFile(partition_data);

    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Finalize();
    return 0;
}