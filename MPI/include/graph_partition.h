#ifndef GRAPH_PARTITION_H
#define GRAPH_PARTITION_H

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <metis.h>
#include <mpi.h>

using namespace std;

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


class GraphPartition {
public:
    // Load graph from METIS format file
    static bool loadGraph(const std::string& filename, 
                         std::vector<idx_t>& xadj, 
                         std::vector<idx_t>& adjncy,
                          int& total_edges);

    // Partition graph using METIS
    static void partitionGraph(const std::vector<idx_t>& xadj,
                             const std::vector<idx_t>& adjncy,
                             int num_parts,
                             std::vector<idx_t>& partitions,
                             int& total_edges);

    // Distribute graph parts to MPI processes
    static void distributePartitions(int argc, char* argv[],
                                const std::vector<idx_t>& xadj,
                                const std::vector<idx_t>& adjncy,
                                const std::vector<idx_t>& partitions,
                                GraphPartitionData& partition_data);

    static void computeDistributedSSSP(int source_global, GraphPartitionData& partition_data);

    static void distributeEdgeUpdates(const std::string& filename, std::vector<EdgeUpdate>& local_updates, GraphPartitionData& partition_data);

    static void gatherAndWriteSSSPToFile(const GraphPartitionData& partition_data);

};

// Updated in graph_partition.h
struct DistributedSSSP {
    std::vector<int> local_distances;
    std::unordered_map<idx_t, int> ghost_nodes;  // <global_node_id, owner_rank>
    std::unordered_map<int, std::vector<idx_t>> ghost_edges;  // <owner_rank, [nodes]>
    std::vector<idx_t> local_nodes;  // Global IDs of local nodes
};

#endif