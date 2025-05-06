#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include <string>
#include <utility>

using namespace std;

class Graph {
public:
    int numNodes;
    vector<vector<pair<int, int>>> adj;

    Graph(int n);
    Graph(const string& filename);
    Graph(const Graph& global, const vector<int>& part, int my_rank);

    void addEdge(int u, int v, int weight);
    void updateEdge(int u, int v, int newWeight);
    void removeEdge(int u, int v);
    bool hasEdge(int u, int v) const;
    int getEdgeWeight(int u, int v) const;
};

#endif