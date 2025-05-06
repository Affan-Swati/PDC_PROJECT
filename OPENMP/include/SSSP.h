#ifndef SSSP_H
#define SSSP_H

#include "Graph.h"
#include <vector>
#include <limits>

class SSSP {
public:
    std::vector<int> dist;
    std::vector<int> parent;
    Graph& graph;

    SSSP(Graph& g);
    void compute(int source);
    void computeWithGhosts(int source, const std::vector<int>& ghostVertices);
    void updateDistance(int vertex, int newDist);
    void printPaths(int source) const;
};

#endif
