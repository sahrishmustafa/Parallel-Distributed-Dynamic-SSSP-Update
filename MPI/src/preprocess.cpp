#include "preprocess_graph.h"
#include <iostream>
#include <vector>
#include <string>

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <graph_file.graph> <num_partitions>" << std::endl;
        return 1;
    }

    std::string graph_file = argv[1];
    int num_parts = std::stoi(argv[2]);

    std::vector<idx_t> xadj, adjncy, partitions;
    int total_edges = 0;

    if (!PreProcessGraph::loadGraph(graph_file, xadj, adjncy, total_edges)) {
        return 1;
    }

    PreProcessGraph::partitionGraph(xadj, adjncy, num_parts, partitions, total_edges);

    PreProcessGraph::writePartitionMetadata("graph_partition_info.txt", xadj, adjncy, partitions);

    std::cout << "Partitioning completed and written to files." << std::endl;
    return 0;
}
