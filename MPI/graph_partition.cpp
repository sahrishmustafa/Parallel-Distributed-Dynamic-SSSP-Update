#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <metis.h>
#include "graph_partition.h"
#include <unordered_set>
#include <mpi.h>
#include <functional>
#include <climits>
#include <sstream> 
#include <sys/stat.h>
#include <unistd.h>   
using namespace std;


bool GraphPartition::loadGraph(const std::string& filename, 
                             std::vector<idx_t>& xadj, 
                             std::vector<idx_t>& adjncy) {
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
            std::cerr << "Unexpected end of file at node " << i+1 << std::endl;
            return false;
        }

        std::istringstream iss(line);
        idx_t node_number;
        iss >> node_number;

        // Verify node numbering
        if (node_number != i + 1) {
            std::cerr << "Node numbering error at position " << i 
                     << ": expected " << i+1 << ", got " << node_number << std::endl;
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

    return true;
}


void GraphPartition::partitionGraph(const std::vector<idx_t>& xadj,
                                  const std::vector<idx_t>& adjncy,
                                  int num_parts,
                                  std::vector<idx_t>& partitions) {
    idx_t nvtxs = xadj.size() - 1;
    cout << "# of nodes:"<<nvtxs<<endl;
    idx_t size_adj = adjncy.size();
    cout<<"size of adjcny: "<<size_adj<<endl;

    idx_t ncon = 1;
    partitions.resize(nvtxs);
    idx_t objval;

    // Create pre-parts directory (portable version)
    mkdir("pre-parts", 0777); // Creates directory with full permissions

    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
    options[METIS_OPTION_NUMBERING] = 0;

    if (adjncy.size() != 2 * 150) {
        std::cerr << "ERROR: adjncy size (" << adjncy.size() 
                << ") != 2 × edges (" << 2 * 150 << ")\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
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
        MPI_Abort(MPI_COMM_WORLD, 1);
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
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    // For each node, write its edges to its partition file
    for (idx_t node = 0; node < nvtxs; node++) {
        int part = partitions[node];
        idx_t start = xadj[node];
        idx_t end = xadj[node + 1];
        cout<<"node:"<<node<<" in part:"<<part<<" start:"<<start<<" end:"<<end<<endl;

        part_files[part] << "Node " << node + 1 << " (Original ID): ";
        
        // Write all edges
        for (idx_t j = start; j < end; j++) {
            idx_t neighbor = adjncy[j] + 1; // Convert to 1-based for output
            //cout<<" neighbor:"<<adjncy[j];

            bool is_external = (partitions[adjncy[j]] != part);
            
            part_files[part] << neighbor;
            if (is_external) {
                part_files[part] << "(X)"; // Mark cross-partition edges
            }
            part_files[part] << " ";
        }
        cout<<endl;
        part_files[part] << "\n";
    }

    // Close all files
    for (auto& f : part_files) {
        f.close();
    }

    std::cout << "Stored pre-partitioned adjacency lists in 'pre-parts' directory\n";
}


// Updated distributePartitions in graph_partition.cpp
void GraphPartition::distributePartitions(int argc, char* argv[],
                                        const std::vector<idx_t>& xadj,
                                        const std::vector<idx_t>& adjncy,
                                        const std::vector<idx_t>& partitions,
                                        GraphPartitionData& partition_data) {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    partition_data.my_rank = world_rank;
    partition_data.partitions = partitions;

    if (world_rank == 0) {
        // MASTER PROCESS - ONLY DISTRIBUTES DATA
        std::cout << "\nMASTER (Rank 0) - Distributing partitions using pre-generated files\n";
        std::cout << "Total partitions to distribute: " << world_size-1 << "\n";

        // Verify partition files exist
        for (int dest = 1; dest < world_size; dest++) {
            std::string part_file = "pre-parts/part_" + std::to_string(dest-1) + ".txt";
            if (!std::ifstream(part_file)) {
                std::cerr << "Error: Partition file missing: " << part_file << std::endl;
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        // Distribute partitions from files
        for (int dest = 1; dest < world_size; dest++) {
            std::string part_file = "pre-parts/part_" + std::to_string(dest-1) + ".txt";
            std::ifstream in(part_file);
            
            // Count nodes in this partition file
            idx_t node_count = 0;
            std::string line;
            while (std::getline(in, line)) {
                if (line.find("Node ") == 0) node_count++;
            }
            in.clear();
            in.seekg(0);

            std::cout << "Sending partition " << dest-1 << " to Rank " << dest 
                     << " (" << node_count << " nodes)\n";

            // Send node count
            MPI_Send(&node_count, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);

            // Send node data
            while (std::getline(in, line)) {
                if (line.find("Node ") != 0) continue;
                
                // Parse node info: "Node X (Original ID): Y1 Y2(X) Y3 ..."
                size_t colon_pos = line.find(":");
                idx_t global_node = std::stoi(line.substr(5, colon_pos - 5 - 13)) - 1;
                cout<<"global node:"<<global_node<<" from file"<<dest-1;
                
                // Get edges from original graph data
                idx_t start = xadj[global_node];
                idx_t end = xadj[global_node+1];
                idx_t edge_count = end - start;
                
                // Send node data
                MPI_Send(&global_node, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
                MPI_Send(&edge_count, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
                MPI_Send(&adjncy[start], edge_count, MPI_INT, dest, 0, MPI_COMM_WORLD);

                std::cout << "  Sent node " << global_node << " with " << edge_count << " edges\n";
            }
            in.close();
        }
        std::cout << "MASTER: Distribution complete. Rank 0 holds no graph data.\n";
    } else {
        // WORKER PROCESSES - RECEIVE AND STORE THEIR PARTITION
        std::ofstream outfile("verified" + std::to_string(world_rank) + ".txt");
        if (!outfile.is_open()) {
            std::cerr << "Error creating verification file for rank " << world_rank << std::endl;
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        outfile << "WORKER (Rank " << world_rank << ") - Partition Verification Data\n";
        
        // Receive node count
        idx_t node_count;
        MPI_Recv(&node_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        outfile << "Total nodes received: " << node_count << "\n\n";

        // Initialize data structures
        partition_data.local_nodes.resize(node_count);
        partition_data.xadj.resize(node_count + 1);
        partition_data.xadj[0] = 0;
        
        // Receive nodes and write to file
        for (idx_t i = 0; i < node_count; i++) {
            idx_t global_node;
            MPI_Recv(&global_node, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            partition_data.local_nodes[i] = global_node;
            partition_data.global_to_local[global_node] = i;
            
            idx_t edge_count;
            MPI_Recv(&edge_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            std::vector<idx_t> edges(edge_count);
            MPI_Recv(edges.data(), edge_count, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
            // Store in CSR format
            partition_data.adjncy.insert(partition_data.adjncy.end(), edges.begin(), edges.end());
            partition_data.xadj[i+1] = partition_data.xadj[i] + edge_count;

            // Write node details to file
            outfile << "Node " << global_node << " (local index " << i << ")\n";
            outfile << "  Edge count: " << edge_count << "\n";
            outfile << "  Edges: ";
            for (size_t j = 0; j < edges.size(); j++) {
                outfile << edges[j];
                if (j < edges.size()-1) outfile << ", ";
            }
            outfile << "\n\n";
        }

        // Identify ghost nodes
        partition_data.identifyGhostNodes(partitions);
        
        // Write ghost node information
        outfile << "\nGhost Nodes (" << partition_data.ghost_nodes.size() << " total):\n";
        for (const auto& entry : partition_data.ghost_nodes) {
            outfile << "  Node " << entry.first << " (owned by Rank " << entry.second << ")\n";
        }

        // Write CSR structure summary
        outfile << "\nCSR Structure Summary:\n";
        outfile << "xadj: [";
        for (size_t i = 0; i < partition_data.xadj.size(); i++) {
            outfile << partition_data.xadj[i];
            if (i < partition_data.xadj.size()-1) outfile << ", ";
        }
        outfile << "]\n";
        
        outfile << "adjncy: [";
        for (size_t i = 0; i < partition_data.adjncy.size(); i++) {
            outfile << partition_data.adjncy[i];
            if (i < partition_data.adjncy.size()-1) outfile << ", ";
        }
        outfile << "]\n";

        // Final summary
        outfile << "\nPartition Summary:\n";
        outfile << "  Local nodes: " << node_count << "\n";
        outfile << "  Ghost nodes: " << partition_data.ghost_nodes.size() << "\n";
        outfile << "  Total edges: " << partition_data.adjncy.size() << "\n";
        outfile << "  xadj size: " << partition_data.xadj.size() << "\n";
        outfile << "  adjncy size: " << partition_data.adjncy.size() << "\n";

        outfile.close();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (world_rank == 0) {
        std::cout << "\nAll partitions successfully distributed\n";
    }
}

void GraphPartitionData::identifyGhostNodes(const std::vector<idx_t>& global_partitions) {
    // First pass: identify all ghost nodes
    for (size_t i = 0; i < local_nodes.size(); i++) {
        idx_t global_node = local_nodes[i];
        idx_t start = xadj[i];
        idx_t end = xadj[i+1];
        
        for (idx_t j = start; j < end; j++) {
            idx_t neighbor = adjncy[j];
            int neighbor_rank = global_partitions[neighbor] + 1; // partitions are 0-based
            
            if (neighbor_rank != my_rank) {
                ghost_nodes[neighbor] = neighbor_rank;
            }
        }
    }
    
    // Second pass: build send maps
    for (const auto& entry : ghost_nodes) {
        idx_t node = entry.first;
        int owner_rank = entry.second;
        send_map[owner_rank].push_back(node);
    }
}

// // Updated in graph_partition.cpp
// void GraphPartition::computeDistributedSSSP(int source_global) {
//     int world_size, world_rank;
//     MPI_Comm_size(MPI_COMM_WORLD, &world_size);
//     MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

//     // Only workers participate in computation
//     if (world_rank == 0) {
//         MPI_Finalize();
//         return;
//     }

//     DistributedSSSP sssp;
//     std::unordered_map<idx_t, idx_t> global_to_local;

//     // Phase 1: Build local data structures
//     idx_t node_count;
//     MPI_Recv(&node_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//     sssp.local_nodes.resize(node_count);
//     sssp.local_distances.resize(node_count, INT_MAX);
    
//     // Receive node data and identify ghost nodes
//     for (idx_t i = 0; i < node_count; i++) {
//         idx_t global_node;
//         MPI_Recv(&global_node, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         global_to_local[global_node] = i;
//         sssp.local_nodes[i] = global_node;

//         idx_t edge_count;
//         MPI_Recv(&edge_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//         std::vector<idx_t> edges(edge_count);
//         MPI_Recv(edges.data(), edge_count, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//         // Identify ghost nodes (cross-partition edges)
//         for (idx_t neighbor : edges) {
//             int neighbor_rank = partitions[neighbor] + 1; // partitions are 0-based
//             if (neighbor_rank != world_rank) {
//                 sssp.ghost_nodes[neighbor] = neighbor_rank;
//                 sssp.ghost_edges[neighbor_rank].push_back(neighbor);
//             }
//         }
//     }

//     // Phase 2: Initialize SSSP
//     if (partitions[source_global-1] == world_rank-1) { // Check if source is local
//         auto it = global_to_local.find(source_global-1);
//         if (it != global_to_local.end()) {
//             sssp.local_distances[it->second] = 0;
//         }
//     }

//     // Phase 3: Distributed Bellman-Ford
//     bool updated = true;
//     for (int iter = 0; iter < partitions.size() && updated; iter++) {
//         updated = false;
        
//         // Process local edges
//         for (size_t i = 0; i < sssp.local_nodes.size(); i++) {
//             idx_t global_node = sssp.local_nodes[i];
//             idx_t start = xadj[global_node];
//             idx_t end = xadj[global_node+1];
            
//             for (idx_t j = start; j < end; j++) {
//                 idx_t neighbor = adjncy[j];
//                 int weight = 1; // Assuming unweighted

//                 if (partitions[neighbor] == world_rank-1) { // Local neighbor
//                     auto local_neighbor = global_to_local[neighbor];
//                     if (sssp.local_distances[i] + weight < sssp.local_distances[local_neighbor]) {
//                         sssp.local_distances[local_neighbor] = sssp.local_distances[i] + weight;
//                         updated = true;
//                     }
//                 } else { // Ghost neighbor
//                     int owner_rank = partitions[neighbor] + 1;
//                     int update_data[2] = {neighbor, sssp.local_distances[i] + weight};
//                     MPI_Send(update_data, 2, MPI_INT, owner_rank, 0, MPI_COMM_WORLD);
//                 }
//             }
//         }

//         // Process incoming updates
//         MPI_Status status;
//         int flag;
//         do {
//             MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
//             if (flag) {
//                 int update_data[2];
//                 MPI_Recv(update_data, 2, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, &status);
//                 idx_t global_node = update_data[0];
//                 int new_dist = update_data[1];
                
//                 if (global_to_local.count(global_node)) {
//                     idx_t local_node = global_to_local[global_node];
//                     if (new_dist < sssp.local_distances[local_node]) {
//                         sssp.local_distances[local_node] = new_dist;
//                         updated = true;
//                     }
//                 }
//             }
//         } while (flag);

//         // Global convergence check (workers only)
//         int global_updated;
//         MPI_Allreduce(&updated, &global_updated, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
//         updated = global_updated;
//     }

//     // Phase 4: Result collection (optional)
//     if (world_rank == 1) { // Designate rank 1 as gather node
//         std::vector<int> global_distances(partitions.size(), INT_MAX);
        
//         // Include own results first
//         for (size_t i = 0; i < sssp.local_nodes.size(); i++) {
//             global_distances[sssp.local_nodes[i]] = sssp.local_distances[i];
//         }

//         // Receive from other workers
//         for (int rank = 2; rank < world_size; rank++) {
//             int node_count;
//             MPI_Recv(&node_count, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
//             std::vector<idx_t> nodes(node_count);
//             std::vector<int> dists(node_count);
//             MPI_Recv(nodes.data(), node_count, MPI_INT, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             MPI_Recv(dists.data(), node_count, MPI_INT, rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            
//             for (int i = 0; i < node_count; i++) {
//                 if (dists[i] < global_distances[nodes[i]]) {
//                     global_distances[nodes[i]] = dists[i];
//                 }
//             }
//         }

//         // Print results
//         std::cout << "Final SSSP distances from node " << source_global << ":\n";
//         for (size_t i = 0; i < global_distances.size(); i++) {
//             if (global_distances[i] != INT_MAX) {
//                 std::cout << "Node " << i+1 << ": " << global_distances[i] << "\n";
//             }
//         }
//     } else if (world_rank > 1) {
//         // Send results to rank 1
//         MPI_Send(&node_count, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);
//         MPI_Send(sssp.local_nodes.data(), node_count, MPI_INT, 1, 0, MPI_COMM_WORLD);
//         MPI_Send(sssp.local_distances.data(), node_count, MPI_INT, 1, 0, MPI_COMM_WORLD);
//     }
//     MPI_Finalize();
// }

