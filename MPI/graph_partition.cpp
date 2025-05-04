#include <iostream>
#include <algorithm>
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


void GraphPartition::partitionGraph(const std::vector<idx_t>& xadj,
                                  const std::vector<idx_t>& adjncy,
                                  int num_parts,
                                  std::vector<idx_t>& partitions,
                                  int& total_edges) {
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

    if (adjncy.size() != 2 * total_edges) {
        std::cerr << "ERROR: adjncy size (" << adjncy.size() 
                << ") != 2 × edges (" << 2 * total_edges << ")\n";
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

          
            // Send node count
            MPI_Send(&node_count, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);

            // Send node data
            while (std::getline(in, line)) {
                if (line.find("Node ") != 0) continue;
                
                // Parse node info: "Node X (Original ID): Y1 Y2(X) Y3 ..."
                size_t colon_pos = line.find(":");
                idx_t global_node = std::stoi(line.substr(5, colon_pos - 5 - 13)) ;//- 1;
               
                // Get edges from original graph data
                idx_t start = xadj[global_node];
                idx_t end = xadj[global_node+1];
                idx_t edge_count = end - start;
                
                // Send node data
                MPI_Send(&global_node, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
                MPI_Send(&edge_count, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);
                MPI_Send(&adjncy[start], edge_count, MPI_INT, dest, 0, MPI_COMM_WORLD);

             
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


void GraphPartition::computeDistributedSSSP(int source_global, GraphPartitionData& data) {
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    ////////////////////////////////
    ////////////////////////////////

    std::ofstream printFile("prints.txt", std::ios::app); // Append mode to combine outputs
    std::streambuf* originalCoutBuffer = std::cout.rdbuf(); // Save original buffer
    std::cout.rdbuf(printFile.rdbuf()); // Redirect cout to file

    ////////////////////////////////
    ////////////////////////////////

    const int INF = INT_MAX;
    const int MAX_ITER = 500;
    int iter_count = 0;

    std::unordered_map<int, int> last_sent_ghost_distance;
    std::vector<int> local_dist(data.local_nodes.size(), INF);
    std::vector<int> parent(data.local_nodes.size(), -1);

    int local_updated = 1;
    int global_updated = 1;

    bool is_active = (world_rank != 0);
    if (!is_active) {
        local_updated = 0;
        std::cout << "[Root] Participating passively in SSSP to support MPI_Allreduce.\n";
    }

    auto it = data.global_to_local.find(source_global);
    if (it != data.global_to_local.end()) {
        local_dist[it->second] = 0;
        parent[it->second] = source_global;  // root points to itself
        std::cout << "[Rank " << world_rank << "] Initializing source node " << source_global << " with dist 0.\n";
    }

    while (global_updated) {
        if (is_active) {
            iter_count++;
            if (iter_count > MAX_ITER) {
                std::cerr << "[Rank " << world_rank << "] Exceeded max iterations! Breaking.\n";
                break;
            }

            local_updated = 0;
            std::cout << "[Rank " << world_rank << "] Iteration " << iter_count << " started.\n";

            // Step 1: Local relaxation
            for (size_t u = 0; u < data.local_nodes.size(); ++u) {
                int u_dist = local_dist[u];
                if (u_dist == INF) continue;

                int u_global = data.local_nodes[u];
                int start = data.xadj[u];
                int end = data.xadj[u + 1];

                for (int i = start; i < end; ++i) {
                    int v_global = data.adjncy[i];
                    int new_dist = u_dist + 1;

                    if (data.global_to_local.count(v_global)) {
                        int v_local = data.global_to_local[v_global];
                        if (new_dist < local_dist[v_local]) {
                            local_dist[v_local] = new_dist;
                            parent[v_local] = u_global;
                            local_updated = 1;

                            std::cout << "[Rank " << world_rank << "] Relaxed local edge "
                                      << u_global << " → " << v_global << " dist=" << new_dist << "\n";
                        }
                    } else {
                        int target_rank = data.ghost_nodes[v_global];
                        if (!last_sent_ghost_distance.count(v_global) || last_sent_ghost_distance[v_global] > new_dist) {
                            int message[3] = {v_global, new_dist, u_global}; // include parent
                            MPI_Send(message, 3, MPI_INT, target_rank, 0, MPI_COMM_WORLD);
                            last_sent_ghost_distance[v_global] = new_dist;

                            // std::cout << "[Rank " << world_rank << "] ✅ Sent ghost update for " << v_global
                            //           << " to rank " << target_rank << " with dist=" << new_dist
                            //           << " parent=" << u_global << "\n";
                        }
                    }
                }
            }

            // Step 2: Receive ghost updates
            int recv_count = 0;
            MPI_Status status;
            while (true) {
                int flag = 0;
                MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
                if (!flag) break;

                int recv_data[3];
                MPI_Recv(recv_data, 3, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                int v_global = recv_data[0];
                int new_dist = recv_data[1];
                int parent_global = recv_data[2];

                if (data.global_to_local.count(v_global)) {
                    int v_local = data.global_to_local[v_global];
                    if (new_dist < local_dist[v_local]) {
                        local_dist[v_local] = new_dist;
                        parent[v_local] = parent_global;
                        local_updated = 1;
                        recv_count++;

                        // std::cout << "[Rank " << world_rank << "] 📥 Ghost update for " << v_global
                        //           << " dist=" << new_dist << " parent=" << parent_global << "\n";
                    }
                }
            }

        }

        MPI_Allreduce(&local_updated, &global_updated, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);

    }

    std::cout << "\n[Rank " << world_rank << "] Final SSSP distances and parents:\n";
    for (size_t i = 0; i < data.local_nodes.size(); ++i) {
        int g = data.local_nodes[i];
        int d = (local_dist[i] == INF) ? -1 : local_dist[i];
        int p = parent[i];
        
        data.sssp_tree[g] = NodeSSSPInfo{d, p};

        std::cout << "  Node " << g << ": dist=" << d << ", parent=" << p << "\n";
    }

    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
    std::vector<std::tuple<int, int, int>> local_results;

    for (size_t i = 0; i < data.local_nodes.size(); ++i) {
        int g = data.local_nodes[i]; // global node ID
        int d = (local_dist[i] == INF) ? -1 : local_dist[i]; // Replace INF with -1
        int p = parent[i];

        
        local_results.emplace_back(g, d, p);
    }

    // Flatten triples to int vector
    std::vector<int> flat_data;
    for (const auto& [u, d, p] : local_results) {
        flat_data.push_back(u);
        flat_data.push_back(d);
        flat_data.push_back(p);
    }

    if (world_rank == 0) {
        std::vector<std::vector<int>> all_data(world_size);
        all_data[0] = flat_data;

        for (int r = 1; r < world_size; ++r) {
            int count;
            MPI_Recv(&count, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            std::vector<int> buffer(count);
            MPI_Recv(buffer.data(), count, MPI_INT, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            all_data[r] = buffer;
        }

        // Combine all partitions
        std::vector<std::tuple<int, int, int>> combined_results;
        for (int r = 0; r < world_size; ++r) {
            const std::vector<int>& buf = all_data[r];
            for (size_t i = 0; i < buf.size(); i += 3) {
                combined_results.emplace_back(buf[i], buf[i+1], buf[i+2]);
            }
        }

        // Sort by global node ID
        std::sort(combined_results.begin(), combined_results.end(),
            [](const auto& a, const auto& b) {
                return std::get<0>(a) < std::get<0>(b);
            });

        // Write to file
        std::ofstream out("sssp_result.txt");
        if (!out.is_open()) {
            std::cerr << "Failed to open sssp_result.txt for writing\n";
        } else {
            for (const auto& [u, d, p] : combined_results) {
                out << u << "\t" << d << "\t" << p << "\n";
            }
            out.close();
            std::cout << "[Rank 0] Wrote sorted SSSP results to sssp_result.txt\n";
        }

    } else {
        int count = flat_data.size();
        MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(flat_data.data(), count, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }

}



  //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////////////////

// std::vector<std::tuple<int, int, int>> local_results;

//     for (size_t i = 0; i < data.local_nodes.size(); ++i) {
//         int g = data.local_nodes[i]; // global ID
//         int d = (local_dist[i] == INF) ? -1 : local_dist[i];
//         int p = -1; // we aren't storing parent info yet — optional
//         local_results.emplace_back(g, d, p);
//     }

//     int local_count = local_results.size();
//     std::vector<int> flat_data; // u, dist, parent triples
//     for (const auto& [u, d, p] : local_results) {
//         flat_data.push_back(u);
//         flat_data.push_back(d);
//         flat_data.push_back(p);
//     }

//     if (world_rank == 0) {
//         std::vector<std::vector<int>> all_data(world_size);
//         all_data[0] = flat_data;

//         for (int r = 1; r < world_size; ++r) {
//             int count;
//             MPI_Recv(&count, 1, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

//             std::vector<int> buffer(count);
//             MPI_Recv(buffer.data(), count, MPI_INT, r, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//             all_data[r] = buffer;
//         }

//         // Combine and write to file
//         std::vector<std::tuple<int, int, int>> combined_results;
//         for (int r = 0; r < world_size; ++r) {
//             const std::vector<int>& buf = all_data[r];
//             for (size_t i = 0; i < buf.size(); i += 3) {
//                 combined_results.emplace_back(buf[i], buf[i+1], buf[i+2]); // u, dist, parent
//             }
//         }

//         // Sort by node ID (u)
//         std::sort(combined_results.begin(), combined_results.end(),
//                 [](const auto& a, const auto& b) {
//                     return std::get<0>(a) < std::get<0>(b);
//                 });

//         // Write to file
//         std::ofstream out("sssp_result.txt");
//         if (!out.is_open()) {
//             std::cerr << "Failed to open output file!\n";
//             return;
//         }

//         //out << "u\tdist\tparent\n";
//         for (const auto& [u, d, p] : combined_results) {
//             out << u << "\t" << d << "\t" << p << "\n";
//         }
//         out.close();
//         std::cout << "[Rank 0] Wrote sorted SSSP output to sssp_result.txt\n";

//     } else {
//         int count = flat_data.size();
//         MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
//         MPI_Send(flat_data.data(), count, MPI_INT, 0, 1, MPI_COMM_WORLD);
//     }
