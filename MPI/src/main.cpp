#include <fstream>
#include <chrono>
#include "graph_partition.h" // adjust as necessary

int countLinesInFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return -1; // or throw an exception if preferred
    }

    int count = 0;
    std::string line;
    while (std::getline(file, line)) {
        ++count;
    }

    file.close();
    return count;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    
    using clock = std::chrono::high_resolution_clock;
    std::chrono::duration<double> total_time{0};

    GraphPartitionData partition_data;
    std::vector<idx_t> xadj, adjncy, partitions;
    int total_edges = 0;

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


        if (!GraphPartition::loadGraph(argv[1], xadj, adjncy, total_edges)) {
            return 1;
        }

        auto start = clock::now();
        GraphPartition::partitionGraph(xadj, adjncy, num_parts, partitions, total_edges);
        auto end = clock::now();
        total_time += end - start;
    }

    // Broadcast partition mapping
    {
        auto start = clock::now();
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
        auto end = clock::now();
        total_time += end - start;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    auto start = clock::now();
    GraphPartition::distributePartitions(argc, argv, xadj, adjncy, partitions, partition_data);
    auto end = clock::now();
    total_time += end - start;

    MPI_Barrier(MPI_COMM_WORLD);

    int source_global = 1;
    start = clock::now();
    GraphPartition::computeDistributedSSSP(source_global, partition_data);
    end = clock::now();
    total_time += end - start;

    MPI_Barrier(MPI_COMM_WORLD);

    std::string update_file = "dataset/deletions.txt";
    int deletions_size = countLinesInFile(update_file);
    if (world_rank == 0) std::cout << "Total Deletions: " << deletions_size << std::endl;
    
    std::vector<EdgeUpdate> local_updates;
    start = clock::now();
    GraphPartition::distributeEdgeUpdates(update_file, local_updates, partition_data);
    end = clock::now();
    total_time += end - start;

    MPI_Barrier(MPI_COMM_WORLD);

    update_file = "dataset/insertions.txt";
    int insertions_size = countLinesInFile(update_file);
    if (world_rank == 0) std::cout << "Total Insertions: " << insertions_size << std::endl;

    std::vector<EdgeUpdate> local_updates_in;
    start = clock::now();
    GraphPartition::distributeEdgeUpdates_Insertions(update_file, local_updates_in, partition_data);
    end = clock::now();
    total_time += end - start;

    MPI_Barrier(MPI_COMM_WORLD);

    start = clock::now();
    GraphPartition::gatherAndWriteSSSPToFile(partition_data);
    end = clock::now();
    total_time += end - start;

    MPI_Barrier(MPI_COMM_WORLD);

    double local_time = total_time.count(), max_time;
    MPI_Reduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        std::cout << "Total execution time: " << max_time << " seconds" << std::endl;
    }

    MPI_Finalize();
    return 0;
}
