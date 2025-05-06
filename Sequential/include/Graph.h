// Graph.h
#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <string>

using namespace std;

class Graph {
public:
    int numNodes;
    vector<vector<pair<int, int>>> adj;

    Graph(int n);
    Graph(const string& filename);
    void addEdge(int u, int v, int weight);
    void updateEdge(int u, int v, int newWeight);
    void removeEdge(int u, int v);
    bool hasEdge(int u, int v);
    int getEdgeWeight(int u, int v);
};

#endif 