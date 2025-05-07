#ifndef UPDATION_H
#define UPDATION_H

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <mpi.h>

#ifndef IDX_T_DEFINED
#define IDX_T_DEFINED
typedef int idx_t;
#endif


// METIS/MPI type compatibility
static inline MPI_Datatype GetMPI_IDX_T() {
    if (sizeof(idx_t) == sizeof(int)) {
        return MPI_INT;
    } else if (sizeof(idx_t) == sizeof(long)) {
        return MPI_LONG;
    } else if (sizeof(idx_t) == sizeof(long long)) {
        return MPI_LONG_LONG;
    } else {
        std::cerr << "Error: Unsupported idx_t size" << std::endl;
        return MPI_DATATYPE_NULL;
    }
}


struct NodeSSSPInfo {
    int distance;
    int parent;
};

struct EdgeUpdate {
    idx_t u;
    idx_t v;
    int weight;
    bool is_insertion;  // true for insertion, false for deletion
};

class GraphPartitionData {
public:
    std::vector<idx_t> local_nodes;          // Global IDs of local nodes
    std::vector<idx_t> xadj;                 // CSR index array for local nodes
    std::vector<idx_t> adjncy;               // CSR adjacency array
    std::unordered_map<idx_t, idx_t> global_to_local;  // Global ID -> Local index
    std::unordered_map<idx_t, int> ghost_nodes;       // Global ID -> Owner rank
    std::unordered_map<int, std::vector<idx_t>> send_map; // Rank -> nodes to send
    
    // For quick access to partition info
    std::vector<idx_t> partitions;           // Global partition mapping
    int my_rank;

    std::unordered_map<int, NodeSSSPInfo> sssp_tree;

    
    void initializeFromMPI();
    void identifyGhostNodes(const std::vector<idx_t>& global_partitions);
};


class Updation {
public:
    // Distribute graph parts to MPI processes
    static void distributePartitions(const std::vector<idx_t>& xadj,
                                const std::vector<idx_t>& adjncy,
                                const std::vector<idx_t>& partitions,
                                GraphPartitionData& partition_data);

    static void computeDistributedSSSP(int source_global, GraphPartitionData& partition_data);

    static void distributeEdgeUpdates(const std::string& filename, std::vector<EdgeUpdate>& local_updates, GraphPartitionData& partition_data);

    static void gatherAndWriteSSSPToFile(const GraphPartitionData& partition_data);

    static void readPartitionMetadata(const std::string& filename,
                            std::vector<idx_t>& xadj,
                            std::vector<idx_t>& adjncy,
                            std::vector<idx_t>& partitions);
};

#endif // UPDATION_H
