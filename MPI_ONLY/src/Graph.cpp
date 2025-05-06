#include "Graph.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

using namespace std;

Graph::Graph(int n) : numNodes(n) {
    adj.resize(n);
}

Graph::Graph(const string& filename) {
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening file: " << filename << endl;
        throw runtime_error("Failed to open graph file");
    }

    int maxNode = 0;
    vector<tuple<int, int, int>> edges;
    string line;

    while (getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        istringstream iss(line);
        int u, v, w;
        if (!(iss >> u >> v >> w)) continue;
        edges.emplace_back(u, v, w);
        maxNode = max({maxNode, u, v});
    }

    numNodes = maxNode + 1;
    adj.resize(numNodes);

    for (const auto& [u, v, w] : edges) addEdge(u, v, w);

    infile.close();
}

Graph::Graph(const Graph& global, const std::vector<int>& part, int my_rank) {
    numNodes = global.numNodes;
    adj.resize(numNodes);
    for (int u = 0; u < global.numNodes; ++u) {
        if (part[u] == my_rank) {
            for (const auto& [v, w] : global.adj[u]) {
                adj[u].push_back({v, w});
            }
        }
    }
}

void Graph::addEdge(int u, int v, int weight) {
    adj[u].push_back({v, weight});
}

void Graph::updateEdge(int u, int v, int newWeight) {
    for (auto& pair : adj[u]) {
        if (pair.first == v) {
            pair.second = newWeight;
            return;
        }
    }
    addEdge(u, v, newWeight);
}

void Graph::removeEdge(int u, int v) {
    adj[u].erase(
        remove_if(adj[u].begin(), adj[u].end(),
            [v](const pair<int, int>& p) { return p.first == v; }),
        adj[u].end()
    );
}

bool Graph::hasEdge(int u, int v) const {
    for (const auto& pair : adj[u]) {
        if (pair.first == v) return true;
    }
    return false;
}

int Graph::getEdgeWeight(int u, int v) const {
    for (const auto& pair : adj[u]) {
        if (pair.first == v) return pair.second;
    }
    return -1;
}