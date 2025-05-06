#include "Graph.h"
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

Graph::Graph(int n) : numNodes(n) {
    adj.resize(n);
}

Graph::Graph(const string& filename) {
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening graph file: " << filename << endl;
        exit(1);
    }

    int maxNode = 0;
    vector<tuple<int, int, int>> edges;
    string line;

    while (getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;  // skip comments
        istringstream iss(line);
        int u, v, w;
        if (!(iss >> u >> v >> w)) {
            cerr << "Invalid line in graph file: " << line << endl;
            continue;
        }
        edges.emplace_back(u, v, w);
        maxNode = max({maxNode, u, v});
    }

    numNodes = maxNode + 1;
    adj.resize(numNodes);

    for (const auto& [u, v, w] : edges) {
        addEdge(u, v, w);
    }

    infile.close();
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
    adj[u].erase(remove_if(adj[u].begin(), adj[u].end(),
        [v](const pair<int, int>& p) { return p.first == v; }), adj[u].end());
}

bool Graph::hasEdge(int u, int v) {
    for (auto& pair : adj[u]) {
        if (pair.first == v) return true;
    }
    return false;
}

int Graph::getEdgeWeight(int u, int v) {
    for (auto& pair : adj[u]) {
        if (pair.first == v) return pair.second;
    }
    return -1; // Not found
}
