//#include <iostream>
//#include <fstream>
//#include <vector>
//#include <queue>
//#include <unordered_set>
//#include <unordered_map>
//#include <limits>
//#include <sstream>
//#include <set>
//
//using namespace std;
//
//const int INF = numeric_limits<int>::max();
//
//struct Edge {
//    int to;
//    int weight;
//};
//
//using Graph = vector<vector<Edge>>;
//
//void loadGraph(Graph& graph, int& numNodes, const string& filename) {
//    ifstream file(filename);
//    unordered_set<int> nodes;
//    vector<pair<int, int>> edges;
//
//    int u, v;
//    while (file >> u >> v) {
//        edges.emplace_back(u, v);
//        nodes.insert(u);
//        nodes.insert(v);
//    }
//
//    numNodes = *max_element(nodes.begin(), nodes.end()) + 1;
//    graph.resize(numNodes);
//
//    for (auto& [src, dst] : edges) {
//        graph[src].push_back({ dst, 1 }); // weight = 1
//        // If undirected, also add reverse
//        graph[dst].push_back({src, 1});
//    }
//}
//
//void dijkstra(const Graph& graph, int source, vector<int>& dist, vector<int>& parent) {
//    int n = graph.size();
//    dist.assign(n, INF);
//    parent.assign(n, -1);
//    dist[source] = 0;
//
//    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
//    pq.emplace(0, source);
//
//    while (!pq.empty()) {
//        auto [d, u] = pq.top(); pq.pop();
//        if (d > dist[u]) continue;
//
//        for (const auto& edge : graph[u]) {
//            int v = edge.to;
//            int weight = edge.weight;
//            if (dist[u] + weight < dist[v]) {
//                dist[v] = dist[u] + weight;
//                parent[v] = u;
//                pq.emplace(dist[v], v);
//            }
//        }
//    }
//}
//
//// Apply edge insertions and update SSSP using Algorithm 2 & 3
//void processInsertions(Graph& graph, vector<pair<int, int>>& insertions,vector<int>& dist, vector<int>& parent, queue<int>& affected){
//    unordered_set<int> inQueue;
//
//    for (auto& [u, v] : insertions) {
//        graph[u].push_back({ v, 1 }); // insert new edge
//        graph[v].push_back({ u, 1 });
//        if (dist[u] != INF && dist[u] + 1 < dist[v]){
//            dist[v] = dist[u] + 1;
//            parent[v] = u;
//            affected.push(v);
//            inQueue.insert(v);
//        }
//    }
//
//    // Process affected subgraph (like Algorithm 3)
//    while (!affected.empty()) {
//        int u = affected.front(); affected.pop();
//        inQueue.erase(u);
//
//        if (dist[u] == INF) continue;
//
//        for (auto& edge : graph[u]) {
//            int v = edge.to;
//            if (dist[u] + 1 < dist[v]) {
//                dist[v] = dist[u] + 1;
//                parent[v] = u;
//                if (!inQueue.count(v)) {
//                    affected.push(v);
//                    inQueue.insert(v);
//                }
//            }
//        }
//    }
//}
//
//void loadInsertions(const string& filename, vector<pair<int, int>>& insertions) {
//    ifstream file(filename);
//    int u, v;
//    while (file >> u >> v) {
//        insertions.emplace_back(u, v);
//    }
//}
//
//void loadDeletions(const string& filename, vector<pair<int, int>>& deletions) {
//    ifstream file(filename);
//    int u, v;
//    while (file >> u >> v) {
//        deletions.emplace_back(u, v);
//    }
//}
//
//// Disconnect the subtree rooted at 'node'
//void markSubtree(int node, const Graph& graph, vector<int>& dist, vector<int>& parent, unordered_set<int>& affectedSet) {
//    queue<int> q;
//    q.push(node);
//    dist[node] = INF;
//    parent[node] = -1;
//    affectedSet.insert(node);
//
//    while (!q.empty()) {
//        int u = q.front(); q.pop();
//        for (auto& edge : graph[u]) {
//            int v = edge.to;
//            if (parent[v] == u && dist[v] != INF) { 
//                dist[v] = INF;
//                parent[v] = -1;
//                affectedSet.insert(v);
//                q.push(v);
//            }
//        }
//    }
//}
//
//// Apply deletions: remove edge + mark disconnected subtrees
//void processDeletions(Graph& graph, const vector<pair<int, int>>& deletions, vector<int>& dist, vector<int>& parent, queue<int>& affected) {
//    unordered_set<int> affectedSet;
//
//    for (auto& [u, v] : deletions) {
//        // Remove edge from both directions
//        graph[u].erase(remove_if(graph[u].begin(), graph[u].end(), [&](Edge e) { return e.to == v; }), graph[u].end());
//        graph[v].erase(remove_if(graph[v].begin(), graph[v].end(), [&](Edge e) { return e.to == u; }), graph[v].end());
//
//        // If the deleted edge was part of the SSSP tree
//        if (parent[v] == u) {
//            markSubtree(v, graph, dist, parent, affectedSet);
//        }
//        else if (parent[u] == v) {
//            markSubtree(u, graph, dist, parent, affectedSet);
//        }
//    }
//
//    // Add affected nodes to processing queue
//    for (int node : affectedSet) {
//        affected.push(node);
//    }
//}
//
//void writeSSSPToFile(const string& filename, const vector<int>& dist, const vector<int>& parent, int& source) {
//    ofstream out(filename);
//    if (!out.is_open()) {
//        cerr << "Error opening output file: " << filename << endl;
//        return;
//    }
//
//    for (int i = source; i < dist.size(); ++i) {
//        if (dist[i] != INF)
//            out << i << "\t" << dist[i] << "\t" << parent[i] << "\n";
//        else
//            out << i << "\tINF\t-\n";
//    }
//
//    out.close();
//    cout << "Results written to " << filename << endl;
//}
//
//
//int main() {
//    Graph graph;
//    int numNodes = 0;
//    string graphFile = "../../../dataset/dataset1/graph.txt";
//    string insertFile = "../../../dataset/dataset1/insertions.txt";
//    string deleteFile = "../../../dataset/dataset1/deletions.txt";
//    int source = 1;
//
//    loadGraph(graph, numNodes, graphFile);
//    cout << "Loaded graph with " << numNodes << " nodes.\n";
//
//
//    vector<int> dist, parent;
//    dijkstra(graph, source, dist, parent);
//    cout << "Initial SSSP computed from source " << source << ".\n";
//   
//    for (int i = source; i < dist.size(); ++i) {
//        if (dist[i] != INF)
//            cout << i << "\t" << dist[i] << "\t\t" << parent[i] << "\n";
//        else
//            cout << i << "\tINF\t\t-\n";
//    }
//
//    vector<pair<int, int>> deletions;
//    loadDeletions(deleteFile, deletions);
//    cout << "Loaded " << deletions.size() << " deletions.\n";
//
//    vector<pair<int, int>> insertions;
//    loadInsertions(insertFile, insertions);
//    cout << "Loaded " << insertions.size() << " insertions.\n";
//
//    queue<int> affected;
//    processDeletions(graph, deletions, dist, parent, affected);
//    processInsertions(graph, insertions, dist, parent, affected);
//
//    // Now propagate affected nodes (Step 2B)
//    unordered_set<int> inQueue;
//    while (!affected.empty()) {
//        int u = affected.front(); affected.pop();
//        inQueue.erase(u);
//        
//        if (dist[u] == INF) continue;
//
//        for (auto& edge : graph[u]) {
//            int v = edge.to;
//            if (dist[u] != INF && dist[u] + 1 < dist[v]) {
//                dist[v] = dist[u] + 1;
//                parent[v] = u;
//                if (!inQueue.count(v)) {
//                    affected.push(v);
//                    inQueue.insert(v);
//                }
//            }
//        }
//    }
//
//    cout << "Insertions + Deletions processed.\n";
//    cout << "Node\tDistance\tParent\n";
//
//    for (int i = source; i < dist.size(); ++i) {
//        if (dist[i] != INF)
//            cout << i << "\t" << dist[i] << "\t\t" << parent[i] << "\n";
//        else
//            cout << i << "\tINF\t\t-\n";
//    }
//
//    //writeSSSPToFile("../../../testing_visuals/output/output2.txt", dist, parent,source);
//
//    return 0;
//}
