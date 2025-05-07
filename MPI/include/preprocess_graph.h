#ifndef PREPROCESS_GRAPH_H
#define PREPROCESS_GRAPH_H

#include <vector>
#include <string>
#include <unordered_map>
#include <iostream>
#include <metis.h>

using namespace std;

class PreProcessGraph {
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

    static void writePartitionMetadata(const std::string& filename,
                                const std::vector<idx_t>& xadj,
                                const std::vector<idx_t>& adjncy,
                                const std::vector<idx_t>& partitions);

};

#endif // PREPROCESS_GRAPH_H
