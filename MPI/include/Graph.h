#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <utility>

class Graph {
public:
    int numNodes;
    std::vector<std::vector<std::pair<int, int>>> adj;

    Graph(int n);
    Graph(const std::string& filename);
    // New: construct subgraph from global graph and partition info
    Graph(const Graph& global, const std::vector<int>& part, int my_rank);

    void addEdge(int u, int v, int weight);
    void updateEdge(int u, int v, int newWeight);
    void removeEdge(int u, int v);
    bool hasEdge(int u, int v) const;
    int getEdgeWeight(int u, int v) const;
};

#endif