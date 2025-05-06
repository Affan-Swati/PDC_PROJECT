#include "SSSP.h"
#include <iostream>
#include <set>

using namespace std;

SSSP::SSSP(Graph& g) : graph(g) {
    dist.resize(g.numNodes, numeric_limits<int>::max());
    parent.resize(g.numNodes, -1);
}

void SSSP::compute(int source) {
    dist.assign(graph.numNodes, numeric_limits<int>::max());
    parent.assign(graph.numNodes, -1);

    using P = pair<int, int>;
    priority_queue<P, vector<P>, greater<P>> pq;
    dist[source] = 0;
    pq.push({0, source});

    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;

        for (auto [v, w] : graph.adj[u]) {
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                parent[v] = u;
                pq.push({dist[v], v});
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
