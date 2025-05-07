#include "Graph.h"
#include "SSSP.h"
#include <mpi.h>
#include <omp.h>
#include <metis.h>
#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <chrono>

using namespace std;
using namespace std::chrono;

// Partition the graph using METIS
void partitionGraph(Graph& g, int numParts, vector<int>& part) {
    idx_t nVertices = g.numNodes;
    idx_t nParts = numParts;
    idx_t objval;
    idx_t ncon = 1;
    vector<idx_t> xadj(nVertices + 1, 0);
    vector<idx_t> adjncy;
    vector<idx_t> adjwgt;
    for (int u = 0; u < nVertices; ++u) {
        for (auto& [v, w] : g.adj[u]) {
            adjncy.push_back(v);
            adjwgt.push_back(w);
        }
        xadj[u + 1] = adjncy.size();
    }
    part.resize(nVertices);
    vector<idx_t> part_metis(nVertices);
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NUMBERING] = 0;
    options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;  // Minimize edge cuts
    int status = METIS_PartGraphKway(
        &nVertices, &ncon, xadj.data(), adjncy.data(), nullptr, nullptr, adjwgt.data(),
        &nParts, nullptr, nullptr, options, &objval, part_metis.data()
    );
    if (status != METIS_OK) {
        cerr << "METIS partitioning failed!\n";
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (int i = 0; i < nVertices; ++i) part[i] = part_metis[i];
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    auto t0 = high_resolution_clock::now();

    auto t_graph_start = high_resolution_clock::now();
    Graph g(1);
    vector<int> part;
    int numNodes;

    if (rank == 0) {
        g = Graph("data/graph.txt");
        numNodes = g.numNodes;
        partitionGraph(g, size, part);
    }

    MPI_Bcast(&numNodes, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        g = Graph(numNodes);
        part.resize(numNodes);
    }
    MPI_Bcast(part.data(), numNodes, MPI_INT, 0, MPI_COMM_WORLD);

    // Broadcast complete graph structure
    if (rank == 0) {
        vector<tuple<int, int, int>> edge_list;
        for (int u = 0; u < numNodes; ++u) {
            for (auto& [v, w] : g.adj[u]) {
                edge_list.push_back({u, v, w});
            }
        }
        int edge_count = edge_list.size();
        MPI_Bcast(&edge_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        for (auto& [u, v, w] : edge_list) {
            MPI_Bcast(&u, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&v, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&w, 1, MPI_INT, 0, MPI_COMM_WORLD);
        }
    } else {
        int edge_count;
        MPI_Bcast(&edge_count, 1, MPI_INT, 0, MPI_COMM_WORLD);
        for (int i = 0; i < edge_count; ++i) {
            int u, v, w;
            MPI_Bcast(&u, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&v, 1, MPI_INT, 0, MPI_COMM_WORLD);
            MPI_Bcast(&w, 1, MPI_INT, 0, MPI_COMM_WORLD);
            g.addEdge(u, v, w);
        }
    }

    auto t_graph_end = high_resolution_clock::now();

    int source = 0;
    SSSP sssp(g);
    sssp.compute(source);

    int num_updates = 10;
    long long edge_update_time = 0;
    long long sssp_time = 0;

    auto t_update_start = high_resolution_clock::now();

    for (int i = 0; i < num_updates; ++i) {
        // Edge update: Update edge weight
        auto t_edge_update_start = high_resolution_clock::now();
        g.updateEdge(i % g.numNodes, (i + 1) % g.numNodes, i + 10);
        auto t_edge_update_end = high_resolution_clock::now();
        edge_update_time += duration_cast<microseconds>(t_edge_update_end - t_edge_update_start).count();

        // Compute SSSP
        auto t_sssp_start = high_resolution_clock::now();
        if(part[i % g.numNodes] == rank && part[(i + 1) % g.numNodes] == rank) {
            sssp.compute(source);
        }
        auto t_sssp_end = high_resolution_clock::now();
        sssp_time += duration_cast<microseconds>(t_sssp_end - t_sssp_start).count();

        // Edge update: Add new edge
        t_edge_update_start = high_resolution_clock::now();
        g.addEdge((i + 2) % g.numNodes, (i + 3) % g.numNodes, i + 5);
        t_edge_update_end = high_resolution_clock::now();
        edge_update_time += duration_cast<microseconds>(t_edge_update_end - t_edge_update_start).count();

        // Compute SSSP again after edge addition
        t_sssp_start = high_resolution_clock::now();
        if(part[(i + 2) % g.numNodes] == rank && part[(i + 3) % g.numNodes] == rank) {
            sssp.compute(source);
        }
        t_sssp_end = high_resolution_clock::now();
        sssp_time += duration_cast<microseconds>(t_sssp_end - t_sssp_start).count();

        // Edge update: Remove an edge
        t_edge_update_start = high_resolution_clock::now();
        g.removeEdge((i + 4) % g.numNodes, (i + 5) % g.numNodes);
        t_edge_update_end = high_resolution_clock::now();
        edge_update_time += duration_cast<microseconds>(t_edge_update_end - t_edge_update_start).count();

        // Compute SSSP again after edge removal
        t_sssp_start = high_resolution_clock::now();
        if(part[(i + 4) % g.numNodes] == rank && part[(i + 5) % g.numNodes] == rank) {
            sssp.compute(source);
        }
        t_sssp_end = high_resolution_clock::now();
        sssp_time += duration_cast<microseconds>(t_sssp_end - t_sssp_start).count();
    }

    auto t_update_end = high_resolution_clock::now();
    auto t1 = high_resolution_clock::now();

    if (rank == 0) {
        cout << "Graph load time: " << duration_cast<seconds>(t_graph_end - t_graph_start).count() << " seconds\n";
        cout << "Edge update total time: " << edge_update_time / 1e6 << " seconds\n";
        cout << "SSSP total time: " << sssp_time / 1e6 << " seconds\n";
        cout << "Total update loop time: " << duration_cast<seconds>(t_update_end - t_update_start).count() << " seconds\n";
        cout << "Total program time: " << duration_cast<seconds>(t1 - t0).count() << " seconds\n";

        cout << "Final SSSP distances (MPI + OPENMP):\n";
        for (int u = 0; u < g.numNodes; ++u) {
            cout << u << ": " << sssp.dist[u] << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
