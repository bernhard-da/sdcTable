#include <Rcpp.h>
using namespace Rcpp;
RcppExport SEXP normQuader2(SEXP indices, SEXP nDims, SEXP lVec) {
	BEGIN_RCPP
		int i, j, index;
		int lV = as<int> ( lVec ) ;
		int nD = as<int> ( nDims );
		int N = lV / nD;
		
		NumericVector indices (clone(indices));
		
		for (i=2; i <= N; i++) {
   		for (j=0; j < nD; j++) {
      		index = (i-1)*nD+j;
            if(indices[index] == indices[j]) {
            	indices[index] = 0;
            }
            else {
            	indices[index] = 1;
            }
  			}
  		}
  		for (j=0; j < nD; j++) {
   		indices[j] = 0;
  		}
   return indices ;
	END_RCPP
}