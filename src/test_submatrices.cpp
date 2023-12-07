#include <string>
#include <vector>

#include <R.h>
#include <Rcpp.h>

#include "util.h"
#include "matrix.h"

extern "C" void test_submat();

using namespace Rcpp;

void test_submat()
{
    // Testing sub_cols_mins
    // A is matrix 6x6 stored in col-major format
    // A is similar to R:matrix(1:36, 6, 6)
    int p = 6;
    NumericVector v(36);
    Range r(1, 36);
    v.assign(r.begin(), r.end());
    NumericMatrix A(p, p, v.begin());

    std::vector<double> juju = as<std::vector<double>>(A);
    double x = A(2, 4);  // should be 27
    double xx = A(1, 5); // should be 32
    

    return;
}