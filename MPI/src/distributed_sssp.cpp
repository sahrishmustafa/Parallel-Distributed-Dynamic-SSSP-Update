#include "updation.h"
#include <mpi.h>
#include <iostream>
#include <vector>
#include <numeric>

using namespace std;

#ifndef IDX_T_DEFINED
#define IDX_T_DEFINED
typedef int idx_t;
#endif

// Structure to hold timing data
struct TimingData {
    double read_metadata = 0.0;
    double distribute_partitions = 0.0;
    double initial_sssp = 0.0;
    double deletions = 0.0;
    double insertions = 0.0;
    double gather_results = 0.0;
};

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);

    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    TimingData local_timing, max_timing, min_timing, avg_timing;
    double start_time;
    GraphPartitionData partition_data;
    std::vector<idx_t> xadj, adjncy, partitions;

    // Section 1: Read metadata
    start_time = MPI_Wtime();
    if (world_rank == 0) {
        Updation::readPartitionMetadata("graph_partition_info.txt", xadj, adjncy, partitions);
    }
    local_timing.read_metadata = MPI_Wtime() - start_time;

    // Section 2: Broadcast partitions
    start_time = MPI_Wtime();
    idx_t part_size;
    if (world_rank == 0) {
        part_size = partitions.size();
    }
    MPI_Bcast(&part_size, 1, GetMPI_IDX_T(), 0, MPI_COMM_WORLD);
    partitions.resize(part_size);
    MPI_Bcast(partitions.data(), part_size, GetMPI_IDX_T(), 0, MPI_COMM_WORLD);
    
    if (part_size < 0 || part_size > 1000000) {
        std::cerr << "Invalid part_size on rank " << world_rank << ": " << part_size << std::endl;
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Section 3: Distribute partitions
    start_time = MPI_Wtime();
    Updation::distributePartitions(xadj, adjncy, partitions, partition_data);
    local_timing.distribute_partitions = MPI_Wtime() - start_time;

    MPI_Barrier(MPI_COMM_WORLD);

    // Section 4: Initial SSSP
    start_time = MPI_Wtime();
    int source_global = 1;
    Updation::computeDistributedSSSP(source_global, partition_data);
    local_timing.initial_sssp = MPI_Wtime() - start_time;

    MPI_Barrier(MPI_COMM_WORLD);

    // Section 5: Process deletions
    start_time = MPI_Wtime();
    std::vector<EdgeUpdate> local_updates;
    Updation::distributeEdgeUpdates("deletions.txt", local_updates, partition_data);
    local_timing.deletions = MPI_Wtime() - start_time;
    cout<<"del:"<<local_timing.deletions<<endl;

    MPI_Barrier(MPI_COMM_WORLD);

    // Section 6: Gather results after deletions
    start_time = MPI_Wtime();
    Updation::gatherAndWriteSSSPToFile(partition_data);
    local_timing.gather_results += MPI_Wtime() - start_time;

    MPI_Barrier(MPI_COMM_WORLD);

    // Section 7: Process insertions
    start_time = MPI_Wtime();
    std::vector<EdgeUpdate> local_updates_in;
    Updation::distributeEdgeUpdates("insertions.txt", local_updates_in, partition_data);
    local_timing.insertions = MPI_Wtime() - start_time;
    cout<<"in:"<<local_timing.insertions<<endl;


    MPI_Barrier(MPI_COMM_WORLD);

    // Section 8: Final gather results
    start_time = MPI_Wtime();
    Updation::gatherAndWriteSSSPToFile(partition_data);
    local_timing.gather_results += MPI_Wtime() - start_time;

    // Collect and analyze timing data across all processes
    MPI_Reduce(&local_timing.insertions, &max_timing.insertions, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_timing.deletions, &min_timing.deletions, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    // Repeat for all timing fields...

    if (world_rank == 0) {
        avg_timing.insertions /= world_size;
        avg_timing.deletions /= world_size;
        // Calculate averages for other fields...



        // Print other sections similarly...
        
        double total_max = max_timing.deletions + max_timing.insertions;
        cout << "Total Maximum Time: " << total_max << "s\n";
    }

    MPI_Finalize();
    return 0;
}