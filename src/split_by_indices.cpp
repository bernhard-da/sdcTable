#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
/*
  strInfo should be a List of Start/End-indices
  as available in dimInfo-objects in slot 'strInfo'
*/
List cpp_splitByIndices(std::vector<std::string> strings, List strInfo) {
  int nDims=strInfo.size();
  int nrElements=strings.size();
  List out;
  IntegerVector intv;
  for ( int j=0; j<nDims; j++ ) {
    intv=strInfo[j];
    int cur_from=intv[0]-1;
    int cur_len=intv[1]-intv[0]+1;
    std::vector<std::string> tmp;
    for ( int i=0; i<nrElements; i++ ) {
      tmp.push_back(strings[i].substr(cur_from, cur_len));
    }
    out.push_back(tmp);
  }
  return(out);
}

