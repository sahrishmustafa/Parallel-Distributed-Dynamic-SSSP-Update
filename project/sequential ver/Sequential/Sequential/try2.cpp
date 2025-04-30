#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <algorithm>

using namespace std;

const int INF = numeric_limits<int>::max();

struct Edge {
    int to;
    int weight;
};

unordered_map<int, vector<Edge>> graph;
unordered_map<int, int> dist;
unordered_map<int, int> parent;
unordered_map<int, bool> affected;
unordered_map<int, bool> affectedDel;

void addEdge(int u, int v) {
    graph[u].push_back({ v, 1 });
    graph[v].push_back({ u, 1 });
}

void removeEdge(int u, int v) {
    auto& edgesU = graph[u];
    edgesU.erase(remove_if(edgesU.begin(), edgesU.end(), [&](Edge e) { return e.to == v; }), edgesU.end());

    auto& edgesV = graph[v];
    edgesV.erase(remove_if(edgesV.begin(), edgesV.end(), [&](Edge e) { return e.to == u; }), edgesV.end());
}

void loadGraph(const string& filename) {
    ifstream fin(filename);
    int u, v;
    while (fin >> u >> v) {
        addEdge(u, v);
    }
}

void initializeSSSP(int src) {
    dist.clear();
    parent.clear();
    for (auto& [node, _] : graph) {
        dist[node] = INF;
        parent[node] = -1;
    }

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
    dist[src] = 0;
    pq.push({ 0, src });

    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        for (const Edge& edge : graph[u]) {
            int v = edge.to;
            if (dist[u] != INF && dist[v] > dist[u] + 1) {
                dist[v] = dist[u] + 1;
                parent[v] = u;
                pq.push({ dist[v], v });
            }
        }
    }
}

void processDeletions(const string& filename) {
    ifstream fin(filename);
    int u, v;
    while (fin >> u >> v) {
        if (parent[v] == u) {
            dist[v] = INF;
            parent[v] = -1;
            affectedDel[v] = true;
            affected[v] = true;
        }
        else if (parent[u] == v) {
            dist[u] = INF;
            parent[u] = -1;
            affectedDel[u] = true;
            affected[u] = true;
        }
        removeEdge(u, v);
    }
}

void processInsertions(const string& filename) {
    ifstream fin(filename);
    int u, v;
    while (fin >> u >> v) {
        addEdge(u, v);
        if (dist.find(u) != dist.end() && dist.find(v) != dist.end()) {
            if (dist[v] > dist[u] + 1) {
                dist[v] = dist[u] + 1;
                parent[v] = u;
                affected[v] = true;
            }
            else if (dist[v] != INF && dist[u] > dist[v] + 1) {
                dist[u] = dist[v] + 1;
                parent[u] = v;
                affected[u] = true;
            }
        }
    }
}

void disconnectSubtrees() {
    queue<int> q;
    for (auto& [v, isAffected] : affectedDel){
        if (isAffected) q.push(v);
    }

    while (!q.empty()) {
        int v = q.front(); q.pop();
        for (const Edge& edge : graph[v]) {
            int c = edge.to;
            if (parent[c] == v && dist[c] != INF) {
                dist[c] = INF;
                parent[c] = -1;
                affectedDel[c] = true;
                affected[c] = true;
                q.push(c);
            }
        }
    }
}

//void updateAffectedVertices() {
//    queue<int> q;
//    for (auto& [v, isAffected] : affected) {
//        if (isAffected) q.push(v);
//    }
//
//    while (!q.empty()) {
//        int v = q.front(); q.pop();
//        affected[v] = false;
//        for (const Edge& edge : graph[v]) {
//            int n = edge.to;
//            if (dist[v] != INF && dist[n] > dist[v] + 1){
//                dist[n] = dist[v] + 1;
//                parent[n] = v;
//                affected[n] = true;
//                q.push(n);
//            }
//            else if (dist[n] != INF && dist[v] > dist[n] + 1) {
//                dist[v] = dist[n] + 1;
//                parent[v] = n;
//                affected[v] = true;
//                q.push(v);
//            }
//        }
//    }
//}

void updateAffectedVertices() {
    queue<int> q;
    for (auto& [v, isAffected] : affected) {
        if (isAffected) q.push(v);
    }

    while (!q.empty()) {
        int v = q.front(); q.pop();
        affected[v] = false;

        // Only reconnect if v has been disconnected (dist[v] == INF)
        if (dist[v] == INF) {
            for (const Edge& edge : graph[v]) {
                int n = edge.to;
                if (dist[n] != INF && dist[n] + 1 < dist[v]) {
                    dist[v] = dist[n] + 1;
                    parent[v] = n;
                    affected[v] = true;
                    q.push(v);
                    break;  // stop as soon as you find a better path
                }
            }
        }
        else {
            for (const Edge& edge : graph[v]) {
                int n = edge.to;
                if (dist[n] > dist[v] + 1) {
                    dist[n] = dist[v] + 1;
                    parent[n] = v;
                    affected[n] = true;
                    q.push(n);
                }
                else if (dist[v] > dist[n] + 1) {
                    dist[v] = dist[n] + 1;
                    parent[v] = n;
                    affected[v] = true;
                    q.push(v);
                }
            }
        }
    }
}


void writeSSSPToFile(int src, const string& filename = "output2.txt") {
    ofstream fout(filename);
    if (!fout.is_open()) {
        cerr << "Failed to open " << filename << " for writing.\n";
        return;
    }

    for (const auto& entry : dist) {
        int node = entry.first;
        int d = (entry.second == INF) ? -1 : entry.second;
        fout << node << "\t" << d << "\t" << parent[node] << "\n";
    }

    fout.close();
    cout << "Results written to " << filename << endl;
}

void printSSSP(int src) {
    cout << "\nShortest paths from source " << src << ":\n";
    for (auto& [node, d] : dist) {
        cout << "Node: " << node << ", Distance: " << (d == INF ? -1 : d) << ", Parent: " << parent[node] << endl;
    }
}

int main() {
    int source = 1;
    string graphFile = "../../../dataset/dataset1/graph.txt";
    string insertFile = "../../../dataset/dataset1/insertions.txt";
    string deleteFile = "../../../dataset/dataset1/deletions.txt";
    string outputFile = "../../../testing_visuals/output/output2.txt";

    loadGraph(graphFile);
    initializeSSSP(source);

    processDeletions(deleteFile);
    processInsertions(insertFile);

    disconnectSubtrees();
    updateAffectedVertices();

    printSSSP(source);
    writeSSSPToFile(source, outputFile);

    return 0;
}