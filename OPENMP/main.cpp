#include "Graph.h"
#include "SSSP.h"
#include <metis.h>
#include <omp.h>
#include <iostream>
#include <vector>
#include <set>
#include <chrono>

using namespace std;
using namespace std::chrono;

// Partition the graph using METIS
void partitionGraph(Graph& g, int numParts, std::vector<int>& part) {
    idx_t nVertices = g.numNodes;
    idx_t nParts = numParts;
    idx_t objval;
    idx_t ncon = 1;
    std::vector<idx_t> xadj(nVertices + 1, 0);
    std::vector<idx_t> adjncy;
    std::vector<idx_t> adjwgt;
    for (int u = 0; u < nVertices; ++u) {
        for (auto& [v, w] : g.adj[u]) {
            adjncy.push_back(v);
            adjwgt.push_back(w);
        }
        xadj[u + 1] = adjncy.size();
    }
    part.resize(nVertices);
    std::vector<idx_t> part_metis(nVertices);
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;
    int status = METIS_PartGraphKway(
        &nVertices, &ncon, xadj.data(), adjncy.data(), nullptr, nullptr, adjwgt.data(),
        &nParts, nullptr, nullptr, options, &objval, part_metis.data()
    );
    if (status != METIS_OK) {
        cerr << "METIS partitioning failed!\n";
        exit(1);
    }
    for (int i = 0; i < nVertices; ++i) part[i] = part_metis[i];
}

int main(int argc, char* argv[]) {
    auto t0 = high_resolution_clock::now();

    auto t_graph_start = high_resolution_clock::now();
    string graph_file = "data/large_graph.txt";
    Graph g(graph_file);
    auto t_graph_end = high_resolution_clock::now();

    int num_threads = omp_get_max_threads();
    vector<int> part;
    auto t_partition_start = high_resolution_clock::now();
    partitionGraph(g, num_threads, part);
    auto t_partition_end = high_resolution_clock::now();

    int num_updates = 1000;
    int source = 0;
    auto t_update_start = high_resolution_clock::now();

    vector<int> global_dist(g.numNodes, std::numeric_limits<int>::max());
    long long edge_update_time = 0;
    long long sssp_time = 0;

    #pragma omp parallel num_threads(num_threads)
    {
        int tid = omp_get_thread_num();
        Graph localGraph(g.numNodes);
        for (int u = 0; u < g.numNodes; ++u) {
            if (part[u] == tid) {
                for (auto& [v, w] : g.adj[u]) {
                    localGraph.addEdge(u, v, w);
                }
            }
        }
        SSSP sssp(localGraph);

        long long local_edge_update_time = 0;
        long long local_sssp_time = 0;

        for (int i = 0; i < num_updates; ++i) {
            auto t_edge_update_start = high_resolution_clock::now();
            if (part[i % g.numNodes] == tid) {
                localGraph.updateEdge(i % g.numNodes, (i + 1) % g.numNodes, i + 10);
                auto t_edge_update_end = high_resolution_clock::now();
                local_edge_update_time += duration_cast<microseconds>(t_edge_update_end - t_edge_update_start).count();

                auto t_sssp_start = high_resolution_clock::now();
                sssp.compute(source);
                auto t_sssp_end = high_resolution_clock::now();
                local_sssp_time += duration_cast<microseconds>(t_sssp_end - t_sssp_start).count();
            }

            t_edge_update_start = high_resolution_clock::now();
            if (part[(i + 2) % g.numNodes] == tid) {
                localGraph.addEdge((i + 2) % g.numNodes, (i + 3) % g.numNodes, i + 5);
                auto t_edge_update_end = high_resolution_clock::now();
                local_edge_update_time += duration_cast<microseconds>(t_edge_update_end - t_edge_update_start).count();

                auto t_sssp_start = high_resolution_clock::now();
                sssp.compute(source);
                auto t_sssp_end = high_resolution_clock::now();
                local_sssp_time += duration_cast<microseconds>(t_sssp_end - t_sssp_start).count();
            }

            t_edge_update_start = high_resolution_clock::now();
            if (part[(i + 4) % g.numNodes] == tid) {
                localGraph.removeEdge((i + 4) % g.numNodes, (i + 5) % g.numNodes);
                auto t_edge_update_end = high_resolution_clock::now();
                local_edge_update_time += duration_cast<microseconds>(t_edge_update_end - t_edge_update_start).count();

                auto t_sssp_start = high_resolution_clock::now();
                sssp.compute(source);
                auto t_sssp_end = high_resolution_clock::now();
                local_sssp_time += duration_cast<microseconds>(t_sssp_end - t_sssp_start).count();
            }
        }

        // Merge results into global_dist (take minimum for each node)
        #pragma omp critical
        {
            edge_update_time += local_edge_update_time;
            sssp_time += local_sssp_time;
            for (int u = 0; u < g.numNodes; ++u) {
                if (sssp.dist[u] < global_dist[u]) {
                    global_dist[u] = sssp.dist[u];
                }
            }
        }
    }
    auto t_update_end = high_resolution_clock::now();
    auto t1 = high_resolution_clock::now();

    cout << "Graph load time: " << duration_cast<seconds>(t_graph_end - t_graph_start).count() << " seconds\n";
    cout << "Partitioning time: " << duration_cast<seconds>(t_partition_end - t_partition_start).count() << " seconds\n";
    cout << "Edge update total time: " << edge_update_time / 1e6 << " seconds\n";
    cout << "SSSP total time: " << sssp_time / 1e6 << " seconds\n";
    cout << "Total update loop time: " << duration_cast<seconds>(t_update_end - t_update_start).count() << " seconds\n";
    cout << "Total program time: " << duration_cast<seconds>(t1 - t0).count() << " seconds\n";
    return 0;
}