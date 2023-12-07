#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <vector>
#include <array>
#include <omp.h>
#include <iostream>

#include <Rmath.h>
#include "rgwish.h"
#include <Rcpp.h>
#include "rgwish_optim.h"
#include "util.h"

extern "C" void test_rgwish();

void set_seed(unsigned int seed)
{
    Rcpp::Environment base_env("package:base");
    Rcpp::Function set_seed_r = base_env["set.seed"];
    set_seed_r(seed);
}

void test_rgwish()
{
    // Testing rgwish_sigma
    // All arrays are column-major
    // M is matrix N*M, then the array representing M is A
    // A is array N*M, and M[i, j] is element at row i, col j is at index A[j * N + i]
    // R is 1-based index_in_array
    // We'll use 0-based index_in_array here.
    int dim = 6;

    // For reproducibility, no parallel computing
    omp_set_num_threads(1);

    // G is graph with 6 nodes
    // Edges (2,3), (2, 4), (3, 4), (3, 5)
    int *G = (int *)calloc(dim * dim, sizeof(int));
    std::list<std::array<int, 2>> edges = {{2, 3}, {2, 4}, {3, 4}, {3, 5}};

    for (std::array<int, 2> edge : edges)
    {
        G[index_in_array(edge[1], edge[0], dim)] = 1;
        G[index_in_array(edge[0], edge[1], dim)] = 1;
    }

    int *size_node = (int *)calloc(dim, sizeof(int));
    for (int j = 0; j < dim; j++)
    {
        // The size of node J: the number of neighbors in the graph
        // It's the sum of values in column J
        int nb_neighbors = 0;
        for (int i = 0; i < dim; i++)
        {
            nb_neighbors += G[index_in_array(i, j, dim)];
        }

        size_node[j] = nb_neighbors;
    }

    // Copy-Pasted from the value obtained after sim
    double Ts[36] = {
        0.44634516019764375, 0, 0, 0, 0, 0, -0.31033079003969677, 0.32620258455677181, 0, 0, 0, 0, 0.033466329196080419, 0.055322323558783706, 0.41061830086713985, 0, 0,
        0, -0.028951167486129126, -0.13820837175960685, 0.17411887254534961, 0.3762944457024982, 0, 0, 0.10873283169019725, -0.086249845545013082,
        -0.38590660775221153, -0.15336507588866624, 0.36109356855701669, 0, -0.09374389754174213, -0.10957031308262671, -0.030969127202672163,
        -0.069314528765121758, -0.27830829752347497, 0.16274398975626489};

    // Constants
    int b_star = 21;
    int p = dim;
    double threshold = 1e-8;

    set_seed(42);
    // Now the output-only variables
    double *K = (double *)calloc(p * p, sizeof(double));
    rgwish_c(G, Ts, K, &b_star, &p, &threshold);

    set_seed(42);
    // Now the output-only variables
    double *K1 = (double *)calloc(p * p, sizeof(double));
    double *sigma1 = (double *)calloc(p * p, sizeof(double));
    GWishart gw(b_star, p, threshold);
    gw.rgwish_c(G, Ts, K1, sigma1);

    std::vector<double> diff(dim * dim);
    for (int i = 0; i < dim * dim; i++)
    {
        diff[i] = fabs(K[i] - K1[i]);
    }

    // The test is successful when the result returned by the 2 methods is similar
    // We consider a maximum difference of 1e-2
    double max_diff = *max_element(diff.begin(), diff.end());
    if (max_diff > 1e-2)
    {
        Rprintf("Test Failed. Too large difference.");
    }
    else
    {
        Rprintf("Test is OK.");
    }

    return;
}