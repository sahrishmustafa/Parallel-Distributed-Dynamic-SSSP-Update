#include "graph_partition.h"
#include <iostream>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
        
    std::vector<idx_t> xadj, adjncy, partitions;
    
    if (world_rank == 0){
        if (argc < 2){
            std::cerr << "Usage: " << argv[0] << " <graph_file.graph> <num_partitions>" << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        int num_parts = std::stoi(argv[2]);
        std::string filename = argv[1];
        filename = argv[1];

        if (num_parts != world_size - 1) {
            std::cerr << "Error: Number of partitions (" << num_parts 
                    << ") must match number of workers (" << world_size - 1 << ")\n";
            MPI_Abort(MPI_COMM_WORLD, 1);
            return 1;
        }

        if (!GraphPartition::loadGraph(filename, xadj, adjncy)) {
            return 1;

        }
        GraphPartition::partitionGraph(xadj, adjncy, num_parts, partitions);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    GraphPartition::distributePartitions(argc, argv, xadj, adjncy, partitions);

    MPI_Barrier(MPI_COMM_WORLD);

    // if (world_rank == 0) {
    //     // Master only initiates computation
    //     int source_global = 1; // Set source node
    //     computeDistributedSSSP(source_global);
    // } else {
    //     // Workers participate in computation
    //     computeDistributedSSSP(-1); // Source irrelevant for workers
    // }
    
    return 0;
}