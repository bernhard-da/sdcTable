#include <Rcpp.h>
using namespace Rcpp;
RcppExport SEXP normQuader2(SEXP indices, SEXP nDims, SEXP lVec) {
  BEGIN_RCPP
    int i, j, index;
    int lV = as<int> ( lVec ) ;
    int nD = as<int> ( nDims );
    int N = lV / nD;

    NumericVector indices2;
    indices2 = (clone(indices));

    for (i=2; i <= N; i++) {
      for (j=0; j < nD; j++) {
        index = (i-1)*nD+j;
        if(indices2[index] == indices2[j]) {
          indices2[index] = 0;
        }
        else {
          indices2[index] = 1;
        }
      }
    }
    for (j=0; j < nD; j++) {
      indices2[j] = 0;
    }
    return indices2;
  END_RCPP
}
