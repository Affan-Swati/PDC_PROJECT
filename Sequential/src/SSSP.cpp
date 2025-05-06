#include "SSSP.h"
#include <iostream>
#include <set>

using namespace std;

SSSP::SSSP(Graph& g) : graph(g) {
    dist.resize(g.numNodes, numeric_limits<int>::max());
    parent.resize(g.numNodes, -1);
}

void SSSP::compute(int source) {
    // Don't reset ALL distances at start
    if (dist[source] > 0) {
        dist[source] = 0;
    }

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;
    pq.push({0, source});

    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;

        for (auto& [v, w] : graph.adj[u]) {
            int newDist = dist[u] + w;
            if (newDist < dist[v]) {
                dist[v] = newDist;
                parent[v] = u;
                pq.push({newDist, v});
            }
        }
    }
}

void SSSP::updateAfterEdgeChange(int source, int u, int v) {
    compute(source);
}

void SSSP::printPaths(int source) {
    for (int i = 0; i < graph.numNodes; ++i) {
        cout << "Node " << i << ": Distance = " << dist[i] << ", Parent = " << parent[i] << "\n";
    }
}
