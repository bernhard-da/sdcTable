#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP mySplitFn(SEXP stringvec, SEXP indices) {
  BEGIN_RCPP
    Rcpp::CharacterVector origStringvec(stringvec);       
    Rcpp::NumericVector origIndices(indices);    
    Rcpp::CharacterVector result(origStringvec.size());   

    int lenStr = result.size();
    int lenIndices = origIndices.size();

    std::string str;
    std::string a;
    for ( int i=0; i < lenStr; i++ ) {
      str.clear();
      a.clear();
      a = origStringvec[i];
      for( int j=0; j < lenIndices; j++ ) {
        str.append(a.substr(origIndices[j], 1));			
      }
      result[i] = str;
    }
    return result;
  END_RCPP
}
