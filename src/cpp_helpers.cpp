#include <Rcpp.h>
using namespace Rcpp;

/*
  strInfo should be a List of Start/End-indices
  as available in dimInfo-objects in slot 'strInfo'

  The idea is to split a 'strID' into the default codes that
  are specified by the indices in input strInfo.

  The function returns a list of n-Elements (length of strInfo)
  with each element holding default codes of one dimension.
*/
// [[Rcpp::export]]
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

/*
  This function takes a CharacterVector and combines it, eventually
  using a seperator. If seperator equals NA, no seperator is used
*/
// [[Rcpp::export]]
CharacterVector cpp_myPaste(CharacterVector stringvec, int nrKeyVars, CharacterVector seperator) {
  CharacterVector stringVec(stringvec);
  int nrRows = stringVec.size();
  int by = nrRows / nrKeyVars;
  CharacterVector outVec(by);
  std::string str;
  LogicalVector na_sep=is_na(seperator);
  bool have_sep=true;
  if ( na_sep[0]==true ) {
    have_sep=false;
  }
  std::string sep;
  if ( have_sep ) {
    sep=seperator[0];
  }
  for (int i=0; i<by; i++) {
    str.clear() ;
    for(int j=0; j<nrKeyVars; j++) {
      str.append(stringVec[i+by*j]);
      if ( j < (nrKeyVars-1) and (have_sep) ) {
        str.append(sep);
      }
    }
    outVec[i]=str;
  }
  return outVec;
}

/*
  This function takes a CharacterVector returns a CharacterVector containing
  only the chars specified in input-argument 'indices'
*/

// [[Rcpp::export]]
CharacterVector cpp_mySplit(CharacterVector stringvec, IntegerVector indices) {
  CharacterVector result(stringvec.size());

  int lenStr = result.size();
  int lenIndices = indices.size();

  std::string str;
  std::string a;
  for (int i=0; i<lenStr; i++) {
    str.clear();
    a.clear();
    a = stringvec[i];
    for (int j=0; j<lenIndices; j++) {
      str.append(a.substr(indices[j]-1, 1));
    }
    result[i]=str;
  }
  return result;
}
