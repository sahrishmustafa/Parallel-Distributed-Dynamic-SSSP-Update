#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <metis.h>
#include "preprocess_graph.h"
#include <unordered_set>
#include <map>
#include <functional>
#include <set>
#include <queue>
#include <climits>
#include <sstream> 
#include <sys/stat.h>
#include <unistd.h>   
using namespace std;

const int INF = INT_MAX;


bool PreProcessGraph::loadGraph(const std::string& filename, 
                             std::vector<idx_t>& xadj, 
                             std::vector<idx_t>& adjncy,
                             int& total_edges) {
    std::ifstream in(filename);
    if (!in.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    // Read header line
    idx_t nodes, edges, fmt, ncon;
    in >> nodes >> edges >> fmt >> ncon;

    // Initialize CSR structures
    xadj.resize(nodes + 1);
    xadj[0] = 0;
    adjncy.reserve(fmt == 0 ? edges * 2 : edges);

    std::string line;
    std::getline(in, line); // Read the rest of the header line

    for (idx_t i = 0; i < nodes; i++) {
        if (!std::getline(in, line)) {
            std::cerr << "Unexpected end of file at node " << i << std::endl;
            return false;
        }

        std::istringstream iss(line);
        idx_t node_number;
        iss >> node_number;

        // Verify node numbering
        if (node_number != i + 1) {
            std::cerr << "Node numbering error at position " << i 
                     << ": expected " << i << ", got " << node_number << std::endl;
            return false;
        }

        // Read all neighbors for this node
        idx_t neighbor;
        while (iss >> neighbor) {
            adjncy.push_back(neighbor - 1); // Convert to 0-based
        }

        xadj[i + 1] = adjncy.size();
    }

    // Verify edge count
    size_t expected_edges = (fmt == 0) ? edges * 2 : edges;
    if (adjncy.size() != expected_edges) {
        std::cerr << "Edge count mismatch! Expected " << expected_edges
                 << ", got " << adjncy.size() << std::endl;
        return false;
    }
    total_edges = edges;

    return true;
}


void PreProcessGraph::partitionGraph(const std::vector<idx_t>& xadj,
                                  const std::vector<idx_t>& adjncy,
                                  int num_parts,
                                  std::vector<idx_t>& partitions,
                                  int& total_edges) {
    idx_t nvtxs = xadj.size() - 1;
    // cout << "# of nodes:"<<nvtxs<<endl;
    idx_t size_adj = adjncy.size();
    // cout<<"size of adjcny: "<<size_adj<<endl;

    idx_t ncon = 1;
    partitions.resize(nvtxs);
    idx_t objval;

    // Create pre-parts directory (portable version)
    mkdir("pre-parts", 0777); // Creates directory with full permissions

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NUMBERING] = 0;

    if (adjncy.size() != 2 * total_edges) {
        std::cerr << "ERROR: adjncy size (" << adjncy.size() 
                << ") != 2 × edges (" << 2 * total_edges << ")\n";
        return;
    }

    int result = METIS_PartGraphKway(
        &nvtxs, &ncon, 
        const_cast<idx_t*>(xadj.data()),
        const_cast<idx_t*>(adjncy.data()),
        nullptr, nullptr, nullptr,
        &num_parts,
        nullptr, nullptr, options,
        &objval, partitions.data()
    );

    if (result != METIS_OK) {
        std::cerr << "METIS partitioning failed" << std::endl;
        return;
    }

    // Store partition mapping
    std::ofstream out("partition_output.part");
    for (auto p : partitions) out << p << "\n";
    out.close();

    // Create adjacency lists for each partition
    std::vector<std::ofstream> part_files(num_parts);
    for (int p = 0; p < num_parts; p++) {
        std::string filename = "pre-parts/part_" + std::to_string(p) + ".txt";
        part_files[p].open(filename);
        if (!part_files[p]) {
            std::cerr << "Error creating file: " << filename << std::endl;
            return;
        }
    }

    // For each node, write its edges to its partition file
    for (idx_t node = 0; node < nvtxs; node++) {
        int part = partitions[node];
        idx_t start = xadj[node];
        idx_t end = xadj[node + 1];

        part_files[part] << "Node " << node << " (Original ID): ";
        
        // Write all edges
        for (idx_t j = start; j < end; j++) {
            idx_t neighbor = adjncy[j] ;//+1; // Convert to 0-based for output
            //cout<<" neighbor:"<<adjncy[j];

            bool is_external = (partitions[adjncy[j]] != part);
            
            part_files[part] << neighbor;
            if (is_external) {
                part_files[part] << "(X)"; // Mark cross-partition edges
            }
            part_files[part] << " ";
        }
        // cout<<endl;
        part_files[part] << "\n";
    }

    // Close all files
    for (auto& f : part_files) {
        f.close();
    }

    std::cout << "Stored pre-partitioned adjacency lists in 'pre-parts' directory\n";
}


void PreProcessGraph::writePartitionMetadata(const std::string& filename,
                            const std::vector<idx_t>& xadj,
                            const std::vector<idx_t>& adjncy,
                            const std::vector<idx_t>& partitions) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Failed to write partition metadata file." << std::endl;
        return;
    }

    out << xadj.size() << "\n";
    for (auto val : xadj) out << val << " ";
    out << "\n";

    out << adjncy.size() << "\n";
    for (auto val : adjncy) out << val << " ";
    out << "\n";

    out << partitions.size() << "\n";
    for (auto val : partitions) out << val << " ";
    out << "\n";

    out.close();
}