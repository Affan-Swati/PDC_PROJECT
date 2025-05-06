#include "Graph.h"
#include "SSSP.h"
#include <mpi.h>
#include <omp.h>
#include <metis.h>
#include <iostream>
#include <vector>
#include <set>
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

    Graph g(1);
    vector<int> part;
    int numNodes;

    auto t_partition_start = high_resolution_clock::now();
    if (rank == 0) {
        g = Graph("data/large_graph.txt");
        numNodes = g.numNodes;
        partitionGraph(g, size, part);
    }
    MPI_Bcast(&numNodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank != 0) {
        g = Graph(numNodes);
        part.resize(numNodes);
    }
    MPI_Bcast(part.data(), numNodes, MPI_INT, 0, MPI_COMM_WORLD);
    auto t_partition_end = high_resolution_clock::now();

    auto t_subgraph_start = high_resolution_clock::now();
    Graph localGraph(numNodes);
    for (int u = 0; u < numNodes; ++u) {
        if (part[u] == rank) {
            for (auto& [v, w] : g.adj[u]) {
                localGraph.addEdge(u, v, w);
            }
        }
    }
    set<int> boundary_nodes;
    for (int u = 0; u < numNodes; ++u) {
        if (part[u] == rank) {
            for (auto& [v, w] : g.adj[u]) {
                if (part[v] != rank) boundary_nodes.insert(u);
            }
        }
    }
    auto t_subgraph_end = high_resolution_clock::now();

    SSSP sssp(localGraph);
    int source = 0;
    int num_updates = 1000;
    auto t_update_start = high_resolution_clock::now();

    long long edge_update_time = 0;
    long long sssp_time = 0;
    long long comm_time = 0;

    for (int i = 0; i < num_updates; ++i) {
        int u1 = i % numNodes;
        int u2 = (i + 1) % numNodes;

        auto t_edge_update_start = high_resolution_clock::now();
        bool did_update = false;
        if (part[u1] == rank) {
            localGraph.updateEdge(u1, u2, i + 10);
            did_update = true;
        }
        if (part[(i + 2) % numNodes] == rank) {
            localGraph.addEdge((i + 2) % numNodes, (i + 3) % numNodes, i + 5);
            did_update = true;
        }
        if (part[(i + 4) % numNodes] == rank) {
            localGraph.removeEdge((i + 4) % numNodes, (i + 5) % numNodes);
            did_update = true;
        }
        auto t_edge_update_end = high_resolution_clock::now();
        edge_update_time += duration_cast<microseconds>(t_edge_update_end - t_edge_update_start).count();

        auto t_sssp_start = high_resolution_clock::now();
        if (did_update) {
            sssp.compute(source);
        }
        auto t_sssp_end = high_resolution_clock::now();
        sssp_time += duration_cast<microseconds>(t_sssp_end - t_sssp_start).count();

        auto t_comm_start = high_resolution_clock::now();
        vector<int> local_boundary_dist(numNodes, -1);
        for (int u : boundary_nodes) local_boundary_dist[u] = sssp.dist[u];
        vector<int> global_boundary_dist(numNodes * size, -1);
        MPI_Allgather(local_boundary_dist.data(), numNodes, MPI_INT,
                      global_boundary_dist.data(), numNodes, MPI_INT, MPI_COMM_WORLD);

        // Update ghost nodes if needed (only for your local nodes)
        for (int p = 0; p < size; ++p) {
            if (p == rank) continue;
            for (int u = 0; u < numNodes; ++u) {
                if (part[u] == p && sssp.dist[u] > global_boundary_dist[p * numNodes + u] && global_boundary_dist[p * numNodes + u] != -1) {
                    sssp.dist[u] = global_boundary_dist[p * numNodes + u];
                }
            }
        }
        auto t_comm_end = high_resolution_clock::now();
        comm_time += duration_cast<microseconds>(t_comm_end - t_comm_start).count();
    }
    auto t_update_end = high_resolution_clock::now();

    auto t_reduce_start = high_resolution_clock::now();
    vector<int> globalDist(numNodes, std::numeric_limits<int>::max());
    MPI_Reduce(sssp.dist.data(), globalDist.data(), numNodes, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
    auto t_reduce_end = high_resolution_clock::now();

    auto t1 = high_resolution_clock::now();

    if (rank == 0) {
        cout << "Partitioning time: " << duration_cast<seconds>(t_partition_end - t_partition_start).count() << " seconds\n";
        cout << "Subgraph construction time: " << duration_cast<seconds>(t_subgraph_end - t_subgraph_start).count() << " seconds\n";
        cout << "Edge update total time: " << edge_update_time / 1e6 << " seconds\n";
        cout << "SSSP total time: " << sssp_time / 1e6 << " seconds\n";
        cout << "Communication total time: " << comm_time / 1e6 << " seconds\n";
        cout << "Reduce time: " << duration_cast<seconds>(t_reduce_end - t_reduce_start).count() << " seconds\n";
        cout << "Total update loop time: " << duration_cast<seconds>(t_update_end - t_update_start).count() << " seconds\n";
        cout << "Total program time: " << duration_cast<seconds>(t1 - t0).count() << " seconds\n";
    }
    MPI_Finalize();
    return 0;
}