#include <Rcpp.h>
using namespace Rcpp;

RcppExport SEXP myPasteWithSep(SEXP stringvec, SEXP nrKeyVars, SEXP seperator) {
  BEGIN_RCPP
    CharacterVector stringVec(stringvec);
    int keyVars = as<int>( nrKeyVars ) ;
    int nrRows = stringVec.size();
    int by = nrRows / keyVars;
    Rcpp::CharacterVector sepOrig(seperator);

    CharacterVector outVec( by ) ;
    std::string str;
    std::string sep;
    sep = sepOrig[0];
    for (int i=0; i < by; i++) {
      str.clear() ;
      for(int j=0; j < keyVars; j++) {
        str.append(stringVec[i+by*j]);
        if ( j < (keyVars-1) ) {
          str.append(sep);
        }
      }
      outVec[i] = str;
    }
    return outVec ;
  END_RCPP
}
