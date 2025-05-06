#ifndef SSSP_H
#define SSSP_H

#include "Graph.h"
#include <vector>
#include <queue>
#include <limits>

using namespace std;

class SSSP {
public:
    Graph& graph;
    vector<int> dist;
    vector<int> parent;

    SSSP(Graph& g);
    void compute(int source);
    void updateAfterEdgeChange(int source, int u, int v);
    void printPaths(int source);
};

#endif