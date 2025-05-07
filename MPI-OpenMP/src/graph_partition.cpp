#include <omp.h>
#include <mpi.h>
#include <metis.h>

#include <fstream>
#include <sstream> 
#include <iostream>
#include <unistd.h>   

#include <set>
#include <queue>
#include <vector>
#include <string>
#include <unordered_set>

#include <climits>
#include <algorithm>
#include <sys/stat.h>
#include <functional>

#include "graph_partition.h"

using namespace std;

const int INF = INT_MAX;
const int VERBOSE = 0;

void rebuildCSRWithUpdates(GraphPartitionData& partition_data, const std::vector<std::pair<idx_t, idx_t>>& pending_insertions,
                           const std::vector<std::pair<idx_t, idx_t>>& pending_deletions, int world_rank) {

    int local_vertex_count = partition_data.global_to_local.size();
    std::vector<std::vector<idx_t>> temp_adj(local_vertex_count);

    // Convert old CSR to temp_adj
    for (idx_t u = 0; u < local_vertex_count; ++u) {
        for (idx_t i = partition_data.xadj[u]; i < partition_data.xadj[u + 1]; ++i) {
            idx_t v = partition_data.adjncy[i];
            temp_adj[u].push_back(v);
        }
    }

    // Apply deletions
    for (const auto& [global_u, global_v] : pending_deletions) {
        if (partition_data.global_to_local.count(global_u)) {
            idx_t local_u = partition_data.global_to_local.at(global_u);
            auto& neighbors = temp_adj[local_u];
            neighbors.erase(std::remove(neighbors.begin(), neighbors.end(), global_v), neighbors.end());
        }
    }

    // Apply insertions (avoid duplicates)
    for (const auto& [global_u, global_v] : pending_insertions) {
        if (partition_data.global_to_local.count(global_u)) {
            idx_t local_u = partition_data.global_to_local.at(global_u);
            auto& neighbors = temp_adj[local_u];
            if (std::find(neighbors.begin(), neighbors.end(), global_v) == neighbors.end()) {
                neighbors.push_back(global_v);
            }
        }
    }

    // Rebuild CSR
    partition_data.xadj.clear();
    partition_data.adjncy.clear();
    partition_data.xadj.push_back(0);
    for (const auto& neighbors : temp_adj) {
        partition_data.adjncy.insert(partition_data.adjncy.end(), neighbors.begin(), neighbors.end());
        partition_data.xadj.push_back(partition_data.adjncy.size());
    }

    // std::cout << "[Rank " << world_rank << "] CSR rebuilt with insertions and deletions. xadj size: " << partition_data.xadj.size() << ", adjncy size: " << partition_data.adjncy.size() << std::endl;
}

bool GraphPartition::loadGraph(const std::string& filename, std::vector<idx_t>& xadj, std::vector<idx_t>& adjncy, int& total_edges) {
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

    std::cout << "Total Nodes: " << nodes << std::endl;
    std::cout << "Total Edges: " << total_edges << std::endl;

    return true;
}


void GraphPartition::partitionGraph(const std::vector<idx_t>& xadj, const std::vector<idx_t>& adjncy, int num_parts, std::vector<idx_t>& partitions, int& total_edges) {
    idx_t nvtxs = xadj.size() - 1;
    // cout << "# of nodes:"<<nvtxs<<endl;
    idx_t size_adj = adjncy.size();
    // cout<<"size of adjcny: "<<size_adj<<endl;

    idx_t ncon = 1;
    partitions.resize(nvtxs);
    idx_t objval;

    // Creates directory that has the graph partitons.
    mkdir("pre-parts", 0777); 

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
        // cout<<endl;
        part_files[part] << "\n";
    }

    // Close all files
    for (auto& f : part_files) {
        f.close();
    }

    // std::cout << "Stored pre-partitioned adjacency lists in 'pre-parts' directory\n";
}


// Updated distributePartitions in graph_partition.cpp
void GraphPartition::distributePartitions(int argc, char* argv[], const std::vector<idx_t>& xadj, const std::vector<idx_t>& adjncy, const std::vector<idx_t>& partitions, GraphPartitionData& partition_data) {
                                                
    int world_size, world_rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    partition_data.my_rank = world_rank;
    partition_data.partitions = partitions;

    if (world_rank == 0) {
        // MASTER PROCESS - ONLY DISTRIBUTES DATA
        // std::cout << "\nMASTER (Rank 0) - Distributing partitions using pre-generated files\n";
        // std::cout << "Total partitions to distribute: " << world_size-1 << "\n";

        // Verify partition files exist.
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

        // std::cout << "MASTER: Distribution complete. Rank 0 holds no graph data.\n";
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
        // std::cout << "[Root] Participating passively in SSSP to support MPI_Allreduce.\n";
    }

    auto it = data.global_to_local.find(source_global);
    if (it != data.global_to_local.end()) {
        local_dist[it->second] = 0;
        parent[it->second] = -1;  // root points to itself
        // std::cout << "[Rank " << world_rank << "] Initializing source node " << source_global << " with dist 0.\n";
    }

    while (global_updated) {
        if (is_active) {
            iter_count++;
            if (iter_count > MAX_ITER) {
                std::cerr << "[Rank " << world_rank << "] Exceeded max iterations! Breaking.\n";
                break;
            }

            local_updated = 0;
            // std::cout << "[Rank " << world_rank << "] Iteration " << iter_count << " started.\n";

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

                            // std::cout << "[Rank " << world_rank << "] Relaxed local edge " << u_global << " → " << v_global << " dist=" << new_dist << "\n";
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

    // std::cout << "\n[Rank " << world_rank << "] Final SSSP distances and parents:\n";
    for (size_t i = 0; i < data.local_nodes.size(); ++i) {
        int g = data.local_nodes[i];
        int d = (local_dist[i] == INF) ? -1 : local_dist[i];
        int p = parent[i];
        
        data.sssp_tree[g] = NodeSSSPInfo{d, p};

        // std::cout << "  Node " << g << ": dist=" << d << ", parent=" << p << "\n";
    }

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
        for (int r = 1; r < world_size; ++r) {
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
            // std::cout << "[Rank 0] Wrote sorted SSSP results to sssp_result.txt\n";
        }

    } else {
        int count = flat_data.size();
        MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        MPI_Send(flat_data.data(), count, MPI_INT, 0, 1, MPI_COMM_WORLD);
    }

}

void disconnectSubtree(GraphPartitionData& partition_data, idx_t root_global) {
    // Get local index of the root
    auto it = std::find(partition_data.local_nodes.begin(), partition_data.local_nodes.end(), root_global);
    if (it == partition_data.local_nodes.end()) return;

    idx_t root_local = std::distance(partition_data.local_nodes.begin(), it);
    // std::cout << "[Rank " << partition_data.my_rank << "] root_global: " << root_global << " with distance: " << partition_data.sssp_tree[root_global].distance << std::endl;

    std::unordered_set<idx_t> visited;
    std::queue<idx_t> q;
    q.push(root_local);
    visited.insert(root_local);

    while (!q.empty()) {
        idx_t u_local = q.front(); q.pop();
        idx_t u_global = partition_data.local_nodes[u_local];

        // std::cout << "[Rank " << partition_data.my_rank << "] Disconnecting u_local: " << u_local << " (global: " << u_global << "), original distance: " << partition_data.sssp_tree[u_global].distance << std::endl;

        // Invalidate distance to mark it disconnected
        partition_data.sssp_tree[u_global].distance = INT_MAX;

        // Traverse only children in the SSSP tree
        for (idx_t i = partition_data.xadj[u_local]; i < partition_data.xadj[u_local + 1]; ++i) {
            idx_t neighbor_global = partition_data.adjncy[i];

            // Check if neighbor is local
            auto it_n = std::find(partition_data.local_nodes.begin(), partition_data.local_nodes.end(), neighbor_global);
            if (it_n == partition_data.local_nodes.end()) continue;

            idx_t neighbor_local = std::distance(partition_data.local_nodes.begin(), it_n);

            // Ensure neighbor is child of current node in SSSP tree
            if (partition_data.sssp_tree.count(neighbor_global) &&
                partition_data.sssp_tree[neighbor_global].parent == u_global &&
                !visited.count(neighbor_local)) {

                // std::cout << "[Rank " << partition_data.my_rank << "] --> Adding child: " << neighbor_global << " (local: " << neighbor_local << ") to queue" << std::endl;

                q.push(neighbor_local);
                visited.insert(neighbor_local);
            }
        }

        // std::cout << "[Rank " << partition_data.my_rank << "] Completed processing u_global: " << u_global << "\n" << std::endl;
    }

    // std::cout << "[Rank " << partition_data.my_rank << "] ✅ Disconnected subtree rooted at global node " << root_global << std::endl;
}


// Reads updates from file and distributes to relevant ranks
void GraphPartition::distributeEdgeUpdates(const std::string& filename, std::vector<EdgeUpdate>& local_updates, GraphPartitionData& partition_data) {
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_rank == 0) {
        std::ifstream infile(filename);
        std::string line;
        std::unordered_map<int, std::vector<EdgeUpdate>> updates_by_rank;
        bool is_insert = false;

        if (!infile.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return;
        }

        while (std::getline(infile, line)) {
            is_insert = false;
            std::istringstream iss(line);
            char type;
            idx_t u, v;
            int weight = 1; // default
            iss >> type >> u >> v;

            if (type == 'I') {
                is_insert = true;
            }

            EdgeUpdate upd = {u, v, weight, is_insert};
            int part_u = partition_data.partitions[u] + 1;
            int part_v = partition_data.partitions[v] + 1;

            std::unordered_set<int> target_ranks = {part_u, part_v};
            for (int r : target_ranks) {
                updates_by_rank[r].push_back(upd);
            }
        }

        for (int r = 1; r < world_size; ++r) {
            int count = updates_by_rank[r].size();
            MPI_Send(&count, 1, MPI_INT, r, 100, MPI_COMM_WORLD);
            if (count > 0) {
                MPI_Send(updates_by_rank[r].data(), count * sizeof(EdgeUpdate), MPI_BYTE, r, 101, MPI_COMM_WORLD);
            }
        }

    } else {
        int num_updates;
        MPI_Recv(&num_updates, 1, MPI_INT, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        local_updates.resize(num_updates);
        if (num_updates > 0) {
            MPI_Recv(local_updates.data(), num_updates * sizeof(EdgeUpdate), MPI_BYTE, 0, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        std::set<idx_t> affected_nodes;
        std::vector<std::pair<idx_t, idx_t>> pending_insertions;
        std::vector<std::pair<idx_t, idx_t>> pending_deletions;
        bool reconnecting = false;

        for (const auto& upd : local_updates) {
            idx_t u = upd.u;
            idx_t v = upd.v;

            bool owns_u = partition_data.global_to_local.count(u);
            bool owns_v = partition_data.global_to_local.count(v);

            if (upd.is_insertion) {
                if (owns_u) pending_insertions.emplace_back(u, v);
                if (owns_v) pending_insertions.emplace_back(v, u);
            } else {
                if (owns_u) pending_deletions.emplace_back(u, v);
                if (owns_v) pending_deletions.emplace_back(v, u);

                if (owns_u && partition_data.sssp_tree[upd.u].parent == upd.v) {
                    reconnecting = true;
                    disconnectSubtree(partition_data, upd.u);
                    pending_deletions.emplace_back(upd.u, upd.v);
                }
                else if (owns_v && partition_data.sssp_tree[upd.v].parent == upd.u) {
                    reconnecting = true;
                    disconnectSubtree(partition_data, upd.v);
                    pending_deletions.emplace_back(upd.v, upd.u);
                }
            }
            if (owns_u) affected_nodes.insert(u);
            if (owns_v) affected_nodes.insert(v);
        }

        for (const auto& [node, info] : partition_data.sssp_tree) {
            if (info.distance == INT_MAX) affected_nodes.insert(node);
        }

        if (!pending_insertions.empty() || !pending_deletions.empty()) {
            rebuildCSRWithUpdates(partition_data, pending_insertions, pending_deletions, world_rank);
        }

        std::priority_queue<std::pair<int, idx_t>, std::vector<std::pair<int, idx_t>>, std::greater<>> pq;
        std::unordered_map<idx_t, int> new_distances;
        std::unordered_set<idx_t> visited;

        for (auto node : affected_nodes) {
            int initial_dist = INT_MAX;
            if (partition_data.sssp_tree.count(node)) {
                initial_dist = partition_data.sssp_tree[node].distance;
            } else {
                partition_data.sssp_tree[node] = {initial_dist, -1};
            }
            pq.push({initial_dist, node});
        }

        while (!pq.empty()) {
            auto [dist, u] = pq.top(); pq.pop();
            if (dist > partition_data.sssp_tree[u].distance) continue;

            auto it = partition_data.global_to_local.find(u);
            if (it == partition_data.global_to_local.end()) continue;

            idx_t local_idx = it->second;

            #pragma omp parallel for
            for (int i = partition_data.xadj[local_idx]; i < partition_data.xadj[local_idx + 1]; ++i) {
                idx_t v = partition_data.adjncy[i];

                if (!partition_data.sssp_tree.count(v) || partition_data.sssp_tree[v].distance == INT_MAX) continue;
                int alt = partition_data.sssp_tree[v].distance + 1;

                #pragma omp critical
                if (alt < partition_data.sssp_tree[u].distance) {
                    partition_data.sssp_tree[u] = {alt, (int)v};
                    new_distances[u] = alt;
                    pq.push({alt, u});
                }
            }

            int new_dist = partition_data.sssp_tree[u].distance;

            #pragma omp parallel for
            for (int i = partition_data.xadj[local_idx]; i < partition_data.xadj[local_idx + 1]; ++i) {
                idx_t v = partition_data.adjncy[i];
                int alt = (new_dist == INT_MAX ? INT_MAX : new_dist + 1);

                #pragma omp critical
                if (!partition_data.sssp_tree.count(v) || alt < partition_data.sssp_tree[v].distance) {
                    partition_data.sssp_tree[v] = {alt, (int)u};
                    new_distances[v] = alt;
                    pq.push({alt, v});
                }
            }
        }

        for (const auto& [v, new_dist] : new_distances) {
            if (partition_data.ghost_nodes.count(v)) {
                int owner = partition_data.ghost_nodes.at(v);
                MPI_Send(&v, 1, GetMPI_IDX_T(), owner, 200, MPI_COMM_WORLD);
                MPI_Send(&new_dist, 1, MPI_INT, owner, 201, MPI_COMM_WORLD);
            }
        }

        MPI_Status status;
        int flag;
        bool any_update = true;

        while (any_update) {
            any_update = false;
            while (true) {
                MPI_Iprobe(MPI_ANY_SOURCE, 200, MPI_COMM_WORLD, &flag, &status);
                if (!flag) break;

                idx_t node;
                int dist;
                MPI_Recv(&node, 1, GetMPI_IDX_T(), status.MPI_SOURCE, 200, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                MPI_Recv(&dist, 1, MPI_INT, status.MPI_SOURCE, 201, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                if (!partition_data.sssp_tree.count(node) || dist < partition_data.sssp_tree[node].distance) {
                    partition_data.sssp_tree[node] = {dist, -1};
                    pq.push({dist, node});
                    any_update = true;
                }
            }

            while (!pq.empty()) {
                auto [dist, u] = pq.top(); pq.pop();
                if (visited.count(u)) continue;
                visited.insert(u);

                auto it = partition_data.global_to_local.find(u);
                if (it == partition_data.global_to_local.end()) continue;
                idx_t local_idx = it->second;

                #pragma omp parallel for
                for (int i = partition_data.xadj[local_idx]; i < partition_data.xadj[local_idx + 1]; ++i) {
                    idx_t v = partition_data.adjncy[i];
                    int alt = dist + 1;

                    #pragma omp critical
                    if (!partition_data.sssp_tree.count(v) || alt < partition_data.sssp_tree[v].distance) {
                        partition_data.sssp_tree[v] = {alt, (int)u};
                        new_distances[v] = alt;
                        pq.push({alt, v});
                    }
                }
            }
        }
    }
}

void GraphPartition::distributeEdgeUpdates_Insertions(const std::string& filename, 
                                         std::vector<EdgeUpdate>& local_updates,
                                         GraphPartitionData& partition_data) {
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    if (world_rank == 0) {
        // Read updates from file and separate insertions/deletions
        std::ifstream file(filename);
        std::unordered_map<int, std::vector<EdgeUpdate>> updates_by_rank;
        EdgeUpdate upd;
        char type;
        
        while (file >> type >> upd.u >> upd.v) {
            upd.is_insertion = (type == 'I');
            upd.weight = 1;
            
            // Determine which ranks need this update
            for (int r = 1; r < world_size; r++) {
                if (partition_data.partitions[upd.u] == r-1 || 
                    partition_data.partitions[upd.v] == r-1) {
                    updates_by_rank[r].push_back(upd);
                }
            }
        }

        // Send updates to workers
        for (int r = 1; r < world_size; r++) {
            int count = updates_by_rank[r].size();
            MPI_Send(&count, 1, MPI_INT, r, 100, MPI_COMM_WORLD);
            if (count > 0) {
                MPI_Send(updates_by_rank[r].data(), 
                        count * sizeof(EdgeUpdate),
                        MPI_BYTE, r, 101, MPI_COMM_WORLD);
            }
        }
    } else {
        // Receive updates
        int update_count;
        MPI_Recv(&update_count, 1, MPI_INT, 0, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        local_updates.resize(update_count);
        if (update_count > 0) {
            MPI_Recv(local_updates.data(),
                    update_count * sizeof(EdgeUpdate),
                    MPI_BYTE, 0, 101,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

        // Separate insertions and deletions
        std::vector<EdgeUpdate> insertions, deletions;
        for (const auto& upd : local_updates) {
            if (upd.is_insertion) insertions.push_back(upd);
            else deletions.push_back(upd);
        }


        if (!insertions.empty()) {
            processDistributedInsertions(partition_data, insertions);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

#include <omp.h>  // Don't forget this

void GraphPartition::processDistributedInsertions(GraphPartitionData& data, const std::vector<EdgeUpdate>& insertions) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // 1. Sequential: Update CSR structure
    for (const auto& upd : insertions) {
        if (!upd.is_insertion) continue;

        bool u_local = data.global_to_local.count(upd.u);
        bool v_local = data.global_to_local.count(upd.v);

        if (u_local || v_local) {
            if (u_local) {
                idx_t u_idx = data.global_to_local[upd.u];
                data.adjncy.insert(data.adjncy.begin() + data.xadj[u_idx + 1], upd.v);
                for (size_t i = u_idx + 1; i < data.xadj.size(); i++) {
                    data.xadj[i]++;
                }
            }

            if (v_local) {
                idx_t v_idx = data.global_to_local[upd.v];
                data.adjncy.insert(data.adjncy.begin() + data.xadj[v_idx + 1], upd.u);
                for (size_t i = v_idx + 1; i < data.xadj.size(); i++) {
                    data.xadj[i]++;
                }
            }

            if (data.global_to_local.count(upd.u) && !data.sssp_tree.count(upd.u)) {
                data.sssp_tree[upd.u] = {INF, -1};
            }

            if (data.global_to_local.count(upd.v) && !data.sssp_tree.count(upd.v)) {
                data.sssp_tree[upd.v] = {INF, -1};
            }

            if (u_local && !v_local) {
                data.ghost_nodes[upd.v] = data.partitions[upd.v];
                data.send_map[data.ghost_nodes[upd.v]].push_back(upd.v);
            }
            if (v_local && !u_local) {
                data.ghost_nodes[upd.u] = data.partitions[upd.u];
                data.send_map[data.ghost_nodes[upd.u]].push_back(upd.u);
            }
        }
    }

    // 2. Parallel: Collect all nodes affected by insertions
    std::vector<std::pair<int, idx_t>> local_pq_entries;

    #pragma omp parallel
    {
        std::vector<std::pair<int, idx_t>> thread_entries;

        #pragma omp for nowait
        for (size_t i = 0; i < insertions.size(); ++i) {
            const auto& upd = insertions[i];
            if (!upd.is_insertion) continue;

            auto enqueue_if_valid = [&](idx_t node) {
                if (data.global_to_local.count(node) && data.sssp_tree.count(node)) {
                    thread_entries.emplace_back(data.sssp_tree[node].distance, node);
                }
            };

            enqueue_if_valid(upd.u);
            enqueue_if_valid(upd.v);

            if (data.global_to_local.count(upd.u)) {
                idx_t u_local = data.global_to_local[upd.u];
                for (idx_t j = data.xadj[u_local]; j < data.xadj[u_local + 1]; j++) {
                    enqueue_if_valid(data.adjncy[j]);
                }
            }

            if (data.global_to_local.count(upd.v)) {
                idx_t v_local = data.global_to_local[upd.v];
                for (idx_t j = data.xadj[v_local]; j < data.xadj[v_local + 1]; j++) {
                    enqueue_if_valid(data.adjncy[j]);
                }
            }
        }

        #pragma omp critical
        local_pq_entries.insert(local_pq_entries.end(), thread_entries.begin(), thread_entries.end());
    }

    std::priority_queue<std::pair<int, idx_t>, std::vector<std::pair<int, idx_t>>, std::greater<>> pq;
    for (const auto& entry : local_pq_entries) {
        pq.push(entry);
    }

    // 3. Process updates
    std::unordered_map<idx_t, int> last_processed;

    while (!pq.empty()) {
        auto [dist_u, u] = pq.top(); pq.pop();
        if (last_processed.count(u) && last_processed[u] < dist_u) continue;
        last_processed[u] = dist_u;

        if (!data.global_to_local.count(u)) continue;

        idx_t u_local = data.global_to_local[u];

        for (idx_t i = data.xadj[u_local]; i < data.xadj[u_local + 1]; i++) {
            idx_t v = data.adjncy[i];
            int new_dist = dist_u + 1;

            if (data.global_to_local.count(v)) {
                if (!data.sssp_tree.count(v) || new_dist < data.sssp_tree[v].distance) {
                    data.sssp_tree[v] = {new_dist, (int)u};
                    pq.push({new_dist, v});
                }
            } else if (data.ghost_nodes.count(v)) {
                int owner = data.ghost_nodes[v];
                int msg[3] = {v, new_dist, u};
                MPI_Send(msg, 3, MPI_INT, owner, 200, MPI_COMM_WORLD);
            }
        }
    }

    // 4. Handle ghost node updates
    MPI_Status status;
    int flag, attempts = 0;
    const int max_attempts = 30, max_iterations = insertions.size() * 5;

    for (int iteration = 0; iteration < max_iterations && attempts < max_attempts; iteration++) {
        MPI_Iprobe(MPI_ANY_SOURCE, 200, MPI_COMM_WORLD, &flag, &status);
        if (!flag) {
            attempts++;
            continue;
        }

        attempts = 0;
        int msg[3];
        MPI_Recv(msg, 3, MPI_INT, status.MPI_SOURCE, 200, MPI_COMM_WORLD, &status);
        idx_t node = msg[0];
        int new_dist = msg[1], parent = msg[2];

        if (!data.sssp_tree.count(node) || new_dist <= data.sssp_tree[node].distance) {
            data.sssp_tree[node] = {new_dist, parent};
            if (data.global_to_local.count(node)) {
                pq.push({new_dist, node});
            }
        }
    }
}


// Gathers all the partitions together and gives the final result.
void GraphPartition::gatherAndWriteSSSPToFile(const GraphPartitionData& partition_data) {
    int world_rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    std::vector<std::tuple<idx_t, int, int>> local_entries;
    for (const auto& [u, info] : partition_data.sssp_tree) {
        local_entries.emplace_back(u, info.distance, info.parent);
    }

    // Serialize local entries
    int local_count = local_entries.size();
    std::vector<int> sendbuf(local_count * 3);
    for (int i = 0; i < local_count; ++i) {
        sendbuf[i * 3 + 0] = std::get<0>(local_entries[i]);
        sendbuf[i * 3 + 1] = std::get<1>(local_entries[i]);
        sendbuf[i * 3 + 2] = std::get<2>(local_entries[i]);
    }

    std::vector<int> recv_counts(world_size);
    MPI_Gather(&local_count, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    std::vector<int> displs;
    std::vector<int> recvbuf;
    if (world_rank == 0) {
        displs.resize(world_size);
        int total = 0;
        for (int i = 0; i < world_size; ++i) {
            displs[i] = total;
            recv_counts[i] *= 3;  
            total += recv_counts[i];
        }
        recvbuf.resize(total * 3);
    }

    MPI_Gatherv(sendbuf.data(), local_count * 3, MPI_INT,
                recvbuf.data(), recv_counts.data(), displs.data(), MPI_INT,
                0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        std::vector<std::tuple<idx_t, int, int>> all_entries;
        int total_entries = 0;
        for (int i = 0; i < world_size; ++i) total_entries += recv_counts[i];

        for (int i = 0; i < total_entries; ++i) {
            idx_t u = recvbuf[i * 3 + 0];
            int dist = recvbuf[i * 3 + 1];
            int parent = recvbuf[i * 3 + 2];
            all_entries.emplace_back(u, dist, parent);
        }

        std::sort(all_entries.begin(), all_entries.end());

        std::ofstream out("update_SSSP.txt");
        for (const auto& [u, dist, parent] : all_entries) {
            out << u << " " << dist << " " << parent << "\n";
        }
        
        out.close();
    }

    MPI_Barrier(MPI_COMM_WORLD);
}