#include "SSSP.h"
#include <queue>
#include <iostream>
#include <omp.h>

using namespace std;

SSSP::SSSP(Graph& g) : graph(g) {
    dist.resize(g.numNodes, numeric_limits<int>::max());
    parent.resize(g.numNodes, -1);
}

void SSSP::compute(int source) {
    // Only initialize source distance if not already set
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

void SSSP::computeWithGhosts(int source, const vector<int>& ghostVertices) {
    compute(source);
    
    // Update ghost vertices
    for (int v : ghostVertices) {
        if (dist[v] != numeric_limits<int>::max()) {
            for (auto& [u, w] : graph.adj[v]) {
                int newDist = dist[v] + w;
                if (newDist < dist[u]) {
                    updateDistance(u, newDist);
                }
            }
        }
    }
}

void SSSP::updateDistance(int vertex, int newDist) {
    if (newDist < dist[vertex]) {
        dist[vertex] = newDist;
        
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < graph.adj[vertex].size(); ++i) {
            auto [v, w] = graph.adj[vertex][i];
            int nextDist = newDist + w;
            if (nextDist < dist[v]) {
                updateDistance(v, nextDist);
            }
        }
    }
}

void SSSP::printPaths(int source) const {
    for (int i = 0; i < graph.numNodes; ++i) {
        cout << "Distance to node " << i << ": ";
        if (dist[i] == numeric_limits<int>::max())
            cout << "INF";
        else
            cout << dist[i];
        cout << endl;
    }
}
