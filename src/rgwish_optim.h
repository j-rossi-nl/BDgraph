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

#ifndef rgwish_optim_H
#define rgwish_optim_H

#include <vector>

class GWishart
{
private:
    int p;
    int pxp;
    int b;
    double threshold;

    std::vector<double> sigma_start;
    std::vector<double> sigma_last;
    std::vector<double> beta_star;
    std::vector<double> sigma_start_i;
    std::vector<double> sigma_start_N_i;
    std::vector<int> N_i;
    std::vector<double> sigma_N_i;
    std::vector<double> tmpSigma;

public:
    GWishart(int b, int p, double threshold);

    void rgwish_c(int *G, double *Ts, double *K, double *sigma);
};

#endif
