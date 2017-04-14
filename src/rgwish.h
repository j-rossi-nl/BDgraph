// ----------------------------------------------------------------------------|
//     Copyright (C) 2012-2016 Mohammadi A. and Wit C. E.
//
//     This file is part of BDgraph package.
//
//     BDgraph is free software: you can redistribute it and/or modify it under 
//     the terms of the GNU General Public License as published by the Free 
//     Software Foundation; see <https://cran.r-project.org/web/licenses/GPL-3>.
//
//     Maintainer:
//     Abdolreza Mohammadi: a.mohammadi@rug.nl or a.mohammadi@uvt.nl
// ----------------------------------------------------------------------------|
  
#ifndef rgwish_H
#define rgwish_H

#include <R.h>
#include <Rmath.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#include <R_ext/Utils.h>
#include <sstream>
#include <string>        // std::string, std::to_string
#include <vector>        // for using vector
#include <limits>        // for numeric_limits<double>::max()
#include "matrix.h"

using namespace std;

extern "C" {
	void rwish_c( double Ts[], double K[], int *b, int *p );

	void rgwish_c( int G[], double Ts[], double K[], int *b, int *p );

	void rgwish_sigma( int G[], int size_node[], double Ts[], double K[], double sigma[], int *b_star, int *p,
					double sigma_start[], double inv_C[], double beta_star[], double sigma_i[], 
					vector<double> &sigma_start_N_i, vector<double> &sigma_N_i, vector<int> &N_i );

	void log_exp_mc( int G[], int nu[], int *b, double H[], int *check_H, int *mc, int *p, double f_T[] );
}

#endif