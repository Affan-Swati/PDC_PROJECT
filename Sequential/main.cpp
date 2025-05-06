#include "Graph.h"
#include "SSSP.h"
#include <omp.h>
#include <iostream>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main(int argc, char* argv[]) {
    auto t0 = high_resolution_clock::now();

    auto t_graph_start = high_resolution_clock::now();
    string graph_file = "data/large_graph.txt";
    Graph g(graph_file);
    auto t_graph_end = high_resolution_clock::now();

    int source = 0;
    SSSP sssp(g);

    int num_updates = 1000; 
    auto t_update_start = high_resolution_clock::now();
    long long edge_update_time = 0;
    long long sssp_time = 0;

    for (int i = 0; i < num_updates; ++i) {
        auto t_edge_update_start = high_resolution_clock::now();
        g.updateEdge(i % g.numNodes, (i + 1) % g.numNodes, i + 10);
        auto t_edge_update_end = high_resolution_clock::now();
        edge_update_time += duration_cast<microseconds>(t_edge_update_end - t_edge_update_start).count();

        auto t_sssp_start = high_resolution_clock::now();
        sssp.compute(source);
        auto t_sssp_end = high_resolution_clock::now();
        sssp_time += duration_cast<microseconds>(t_sssp_end - t_sssp_start).count();

        t_edge_update_start = high_resolution_clock::now();
        g.addEdge((i + 2) % g.numNodes, (i + 3) % g.numNodes, i + 5);
        t_edge_update_end = high_resolution_clock::now();
        edge_update_time += duration_cast<microseconds>(t_edge_update_end - t_edge_update_start).count();

        t_sssp_start = high_resolution_clock::now();
        sssp.compute(source);
        t_sssp_end = high_resolution_clock::now();
        sssp_time += duration_cast<microseconds>(t_sssp_end - t_sssp_start).count();

        t_edge_update_start = high_resolution_clock::now();
        g.removeEdge((i + 4) % g.numNodes, (i + 5) % g.numNodes);
        t_edge_update_end = high_resolution_clock::now();
        edge_update_time += duration_cast<microseconds>(t_edge_update_end - t_edge_update_start).count();

        t_sssp_start = high_resolution_clock::now();
        sssp.compute(source);
        t_sssp_end = high_resolution_clock::now();
        sssp_time += duration_cast<microseconds>(t_sssp_end - t_sssp_start).count();
    }
    auto t_update_end = high_resolution_clock::now();
    auto t1 = high_resolution_clock::now();

    cout << "Graph load time: " << duration_cast<seconds>(t_graph_end - t_graph_start).count() << " seconds\n";
    cout << "Edge update total time: " << edge_update_time / 1e6 << " seconds\n";
    cout << "SSSP total time: " << sssp_time / 1e6 << " seconds\n";
    cout << "Total update loop time: " << duration_cast<seconds>(t_update_end - t_update_start).count() << " seconds\n";
    cout << "Total program time: " << duration_cast<seconds>(t1 - t0).count() << " seconds\n";

    // cout << "Final SSSP distances (Sequential):\n";
    // for (int u = 0; u < g.numNodes; ++u) {
    //     cout << u << ": " << sssp.dist[u] << "\n";
    // }

    return 0;
}