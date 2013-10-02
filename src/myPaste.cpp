#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP myPaste(SEXP stringvec, SEXP nrKeyVars) {
  BEGIN_RCPP
    CharacterVector stringVec(stringvec);
    int keyVars = as<int>( nrKeyVars ) ;
    int nrRows = stringVec.size();
    int by = nrRows / keyVars;

    CharacterVector outVec( by ) ;
    std::string str;
    for (int i=0; i < by; i++) {
      str.clear() ;
      for(int j=0; j < keyVars; j++) {
        str.append(stringVec[i+by*j]);			
      }
      outVec[i] = str;
    }
    return outVec ;
  END_RCPP
}
