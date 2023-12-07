// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |
//     Copyright (C) 2012 - 2020  Reza Mohammadi                                                   |
//                                                                                                 |
//     This file is part of BDgraph package.                                                       |
//                                                                                                 |
//     BDgraph is free software: you can redistribute it and/or modify it under                    |
//     the terms of the GNU General Public License as published by the Free                        |
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.                   |
//                                                                                                 |
//     Maintainer: Reza Mohammadi <a.mohammadi@uva.nl>                                             |
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - |

#include <string>
#include <math.h>

#include <R_ext/Lapack.h>

#include "rgwish_optim.h"
#include "matrix.h"
#include "rgwish.h"
#include "util.h"

using namespace std;

/**
 * Encapsulates the operations of

*/
GWishart::GWishart(int b, int p, double threshold) : p(p), pxp(p * p), b(b), threshold(threshold)
{
    this->sigma_start = vector<double>(this->pxp);
    this->sigma_last = vector<double>(this->pxp);
    this->beta_star = vector<double>(p);
    this->sigma_start_i = vector<double>(p);
    this->sigma_start_N_i = vector<double>(p);
    this->N_i = vector<int>(p);
    this->sigma_N_i = vector<double>(this->pxp);
    this->tmpSigma = vector<double>(this->pxp);
}

/*
 * @param G: (INPUT) Graph G adjacency matrix p-by-p with binary values in {0, 1} as COLUMN-MAJOR flat array
 * @param Ts: (INPUT) what is it exactly ?
 * @param K: (OUTPUT) square matrix p-by-p as COLUMN-MAJOR flat array
 * @param sigma: (OUTPUT) square matrix p-by-p as COLUMN-MAJOR flat array, inverse of K
 */
void GWishart::rgwish_c(int *G, double *Ts, double *K, double *sigma)
{
    rwish_c(&Ts[0], K, &b, &p);                           // Random-draw of K
    inverse(K, &sigma_start[0], &p);                      // Sigma_start = inverse of K
    memcpy(sigma, &sigma_start[0], pxp * sizeof(double)); // sigma = sigma_start

    // From Lenkoski, "A Direct Sampler for G-Wishart Variates" (2013) page 6
    // Code provided by the author. Thanks !!
    int counter = 0;
    double mean_diff = 1.0;
    while ((mean_diff > threshold) && (counter < 5000))
    {
        counter++;
        memcpy(&sigma_last[0], sigma, pxp * sizeof(double)); // sigma_last = sigma

        for (int i = 0; i < p; i++)
        {
            // i designates 1 node in graph G
            // Count size of node = #neighbors in the graph
            int size_node = 0;
            for (int j = 0; j < p; j++)
                size_node += G[index_in_array(i, j, p)];

            // If this node has neighbors in G
            if (size_node > 0)
            {
                int l = 0;
                for (int j = 0; j < p; j++)
                {
                    // Scanning the other nodes in graph
                    if (G[index_in_array(i, j, p)])
                    {
                        // Node j is the l-th neighbor of node i
                        sigma_start_N_i[l] = sigma_start[index_in_array(j, i, p)];
                        N_i[l] = j;
                        l++;
                    }
                    else
                        beta_star[j] = 0.0;
                }

                // assert l == size_node
                sub_matrix(&sigma[0], &sigma_N_i[0], &N_i[0], &size_node, &p); // sigma_N_i = sigma[N_i, N_i]

                // A * X = B   for   sigma_start_N_i := (sigma_N_i)^{-1} * sigma_start_N_i
                int info = 0;
                int one = 1;
                char uplo = 'U';
                F77_NAME(dposv)
                (&uplo, &size_node, &one, &sigma_N_i[0], &size_node, &sigma_start_N_i[0], &size_node, &info FCONE);

                if (info != 0)
                {
                    // ERROR
                    exit(1);
                }

                for (int j = 0; j < size_node; j++)
                    beta_star[N_i[j]] = sigma_start_N_i[j];

                char transN = 'N';
                double alpha = 1.0;
                double beta_zero = 0.0;
                F77_NAME(dgemm)
                (&transN, &transN, &p, &one, &p, &alpha, &sigma[0], &p, &beta_star[0], &p, &beta_zero, &sigma_start_i[0], &p FCONE FCONE);

                // Update sigma
                for (int j = 0; j < p; j++)
                {
                    if (i == j)
                        continue;

                    sigma[index_in_array(i, j, p)] = sigma_start_i[j];
                    sigma[index_in_array(j, i, p)] = sigma_start_i[j];
                }
            }
            else // No neighbor in G for node i
            {
                for (int j = 0; j < p; j++)
                {
                    if (i == j)
                        continue;
                    sigma[index_in_array(i, j, p)] = 0.0;
                    sigma[index_in_array(j, i, p)] = 0.0;
                }
            }
        }

        // Mean diff between sigma sigma_last
        // mean(updated sigma after iteration - sigma at start of iteration)
        mean_diff = 0.0;
        for (int i = 0; i < pxp; i++)
            mean_diff += fabs(sigma[i] - sigma_last[i]);
        mean_diff /= pxp;
    }

    // K is the inverse of sigma
    memcpy(&tmpSigma[0], &sigma[0], pxp * sizeof(double));  // tmpSigma = sigma
    inverse(&tmpSigma[0], K, &p);                           // inverse 
}
