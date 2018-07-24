#include <Rcpp.h>
using namespace Rcpp;

/*
  * codes in sdcStatus:
  * 0: possible suppression partner ("s")
  * 1: primary suppression ("u")
  * 2: secondary suppression ("x")
  * 3: non-applicable ("z")
  * 4: dummy-cells ("w")
*/
IntegerVector convert_sdcStatus_to_num(CharacterVector sdcStatus) {
  IntegerVector result(sdcStatus.size());
  for (int i=0; i<sdcStatus.size(); i++) {
    if (sdcStatus[i]=="s") {
      result[i] = 0;
    }
    if (sdcStatus[i]=="u") {
      result[i] = 1;
    }
    if (sdcStatus[i]=="x") {
      result[i] = 2;
    }
    if (sdcStatus[i]=="z") {
      result[i] = 3;
    }
    if (sdcStatus[i]=="w") {
      result[i] = 4;
    }
  }
  return(result);
}

// [[Rcpp::export]]
List greedyMultDimSuppression(DataFrame dat, List indices, List subIndices, IntegerVector dimVars, bool verbose) {
  bool debug=false;
  Function cpp_print("print");

  /* start protection of data() */
  int nrDims=dimVars.size();
  if (verbose == true) {
    Rcout << "We have to protect an " << nrDims << " dimensional dataset!" << std::endl;
  }
  /* extracting input data from list */
  IntegerVector freq=dat["freq"];
  NumericVector weights=dat["weights"];
  CharacterVector sdcStatus=dat["sdcStatus"];
  IntegerVector sdcStatus_num=convert_sdcStatus_to_num(sdcStatus);
  IntegerVector id=dat["id"];

  // is_real_z is true for all cells that are initially "z"
  LogicalVector is_real_z = sdcStatus_num==3;

  IntegerVector current_indices,st_freq,st_sdcStatus_num,st_id;
  NumericVector st_weights;
  CharacterVector st_sdcStatus;
  LogicalVector st_is_real_z;
  List st_subIndices, subList;
  int nrGroups=indices.size();
  int ind_x=-1;
  bool runInd=true;
  int counter=1;
  int total_new_supps=0; /* total number of required secondary supps */
  while (runInd==true ) {
    if (verbose) {
      Rcout << "Start of run " << counter << std::endl;
    }
    bool override=false;
    bool final_ok=true;
    for (int group=0; group<nrGroups; group++) {
      List subList=indices[group];
      int nrTabs=subList.size();
      /* tab is iterating over all subtables within a given group */
      for (int tab=0; tab<nrTabs; tab++) {
        /* subsetting data to current subtable */
        List subList=indices[group];
        current_indices=subList[tab];
        int n_st=current_indices.size();

        /* Select freqs and suppression pattern of current subtable */
        st_freq=freq[current_indices-1];
        st_weights=weights[current_indices-1];
        st_sdcStatus=sdcStatus[current_indices-1];
        st_sdcStatus_num=sdcStatus_num[current_indices-1];
        st_id=id[current_indices-1];
        st_is_real_z=is_real_z[current_indices-1]; /* ids in entire dataset, c-style (starting with 0!) */
        st_id=st_id-1;

        IntegerVector sub_ids=seq_along(st_id);
        sub_ids=sub_ids-1;

        /* create list for indices defining subtables */
        List tmp_subIndices=subIndices[group];
        st_subIndices=tmp_subIndices[tab];

        /* iterativly protecting the current simple tab */
        bool newSuppsAdded=true;
        int suppsAdded=0;
        while (newSuppsAdded==true) {
          int additional_supps=0;
          for (int i=0; i<nrDims; i++) {
            IntegerVector gr=st_subIndices[i];
            int nrSimpleTables=max(gr);
            IntegerVector c_indices=seq_along(gr)-1; // indices relative to current simple table, c-style!

            for (int j=0; j<nrSimpleTables; j++) {
              /* subset original vectors to current subtable */
              /* ind_c: index relative to entire data, starting with 0 */
              IntegerVector ind_c=c_indices[gr==j+1];
              int nCells=ind_c.size();
              if (nCells>1) {
                /* extract values of current simple table */
                IntegerVector cur_freq=st_freq[ind_c];
                NumericVector cur_weights=st_weights[ind_c];
                NumericVector cur_weights_o=clone(cur_weights);
                IntegerVector cur_sdcStatus=st_sdcStatus_num[ind_c];
                IntegerVector cur_id=st_id[ind_c];
                IntegerVector cur_sub_ids=sub_ids[ind_c];
                LogicalVector cur_is_real_z=st_is_real_z[ind_c];
                LogicalVector isCandidate(cur_freq.size());

                if (debug==true) {
                  Rcout << "curSTDF" << std::endl;
                  DataFrame curSTDF = DataFrame::create(
                    Named("cur_freq")=cur_freq,
                    Named("cur_weights")=cur_weights,
                    Named("cur_weights_o")=cur_weights_o,
                    Named("cur_sdcStatus")=cur_sdcStatus,
                    Named("cur_id")=cur_id,
                    Named("cur_sub_ids")=cur_sub_ids,
                    Named("cur_is_real_z")=cur_is_real_z);
                  print(curSTDF);
                }

                /* calculating the number of suppressed cells */
                int nr_supps=0;
                int upVal=max(cur_weights)+1.0;
                int nr_dummycells=0;
                int nr_zcells=0;
                for (int kk=0; kk<nCells; kk++) {
                  isCandidate[kk]=false;
                  if ((cur_sdcStatus[kk] == 1) or (cur_sdcStatus[kk]==2) ) {
                    nr_supps = nr_supps+1;
                    cur_weights[kk]=upVal;
                  }
                  /* we only count 'z' cells as a candiate that were not initially z */
                  if (cur_sdcStatus[kk] == 3) {
                    nr_zcells=nr_zcells+1;
                  }
                  if (cur_sdcStatus[kk] == 4) {
                    nr_dummycells=nr_dummycells+1;
                  }
                  if ((cur_sdcStatus[kk]==0) & (cur_freq[kk]>0)) {
                    isCandidate[kk]=true;
                  }
                }

                /*
                  we have at least 1 observation in the current
                  simple table but only a single suppressed cell and
                  the number of dummy-cells (which should never be published) is 0
                */
                if ((nr_supps == 1) & (nr_dummycells==0)) {
                  int nrCandidates=sum(isCandidate);
                  if (nrCandidates==0) {
                    if (nr_zcells==0) {
                      if (debug) {
                        DataFrame debugDF = DataFrame::create(
                          Named("nCells")=nCells,
                          Named("nr_supps")=nr_supps,
                          Named("nr_dummycells")=nr_dummycells,
                          Named("cur_freq")=cur_freq,
                          Named("cur_weights")=cur_weights,
                          Named("cur_weights_o")=cur_weights_o,
                          Named("cur_sdcStatus")=cur_sdcStatus);
                        print(debugDF);
                      }
                      stop("Unfortunately, it is not possible to find a suppression pattern!");
                    }

                    /*
                      Relaxation case
                      In this case, it is not possible to find a pattern with only 's'-cells.
                      we need to relax 'z' (code 3) cells to 's' (code 0) and try again.
                      However, we do not allow z-cells that were initially "z"
                    */
                    if (debug==true) {
                      Rcout << "we need to set 'z'-cells to 's'!" << std::endl;
                      Rcout << "the current subtable (debugDF):" << std::endl;
                      DataFrame debugDF = DataFrame::create(
                        Named("nCells")=nCells,
                        Named("nr_supps")=nr_supps,
                        Named("nr_dummycells")=nr_dummycells,
                        Named("cur_freq")=cur_freq,
                        Named("cur_weights")=cur_weights,
                        Named("cur_weights_o")=cur_weights_o,
                        Named("cur_sdcStatus")=cur_sdcStatus,
                        Named("cur_real_z")=cur_is_real_z);
                      print(debugDF);
                    }
                    override=true;
                    LogicalVector ii=(cur_sdcStatus==3) & (cur_is_real_z==false);

                    // check if we have candidates; if not we stop with an error
                    // as we do not want to change pre-existing z-cells
                    if (sum(ii)==0) {
                      DataFrame debugDF = DataFrame::create(
                        Named("cur_id")=cur_id,
                        Named("cur_freq")=cur_freq,
                        Named("cur_weights")=cur_weights,
                        Named("cur_weights_o")=cur_weights_o,
                        Named("cur_sdcStatus")=cur_sdcStatus,
                        Named("cur_real_z")=cur_is_real_z);
                      print(debugDF);
                      stop("No valid suppression pattern can be found (too many initial z-cells?)");
                    }

                    cur_weights[ii] = cur_weights_o[ii]; // reset to original value
                    cur_sdcStatus[ii] = 0; // set cur_sdcStatus temporarily to 's' (0)
                    if (debug==true) {
                      Rcout << "relaxdebugDF:" << std::endl;
                      DataFrame relaxdebugDF = DataFrame::create(
                        Named("ii")=ii,
                        Named("cur_freq")=cur_freq,
                        Named("cur_weights")=cur_weights,
                        Named("cur_weights_o")=cur_weights_o,
                        Named("cur_sdcStatus")=cur_sdcStatus,
                        Named("cur_is_real_z")=cur_is_real_z,
                        Named("cur_sub_ids")=cur_sub_ids);
                      print(relaxdebugDF);
                    }
                    isCandidate=(cur_sdcStatus==0);
                  }

                  /* show output on current simple table */
                  if (debug==true) {
                    Rcout << "dfST:" << std::endl;
                    DataFrame dfST = DataFrame::create(
                      Named("cur_id")=cur_id,
                      Named("cur_freq")=cur_freq,
                      Named("cur_sub_ids")=cur_sub_ids,
                      Named("isCandidate")=isCandidate,
                      Named("cur_sdcStatus")=cur_sdcStatus,
                      Named("is_real_z")=cur_is_real_z);
                    print(dfST);
                    R_FlushConsole();
                  }

                  /* Restrict to candidate cells */
                  cur_id=cur_id[isCandidate];
                  cur_freq=cur_freq[isCandidate];
                  cur_weights=cur_weights[isCandidate];
                  cur_sub_ids=cur_sub_ids[isCandidate];
                  cur_sdcStatus=cur_sdcStatus[isCandidate];
                  cur_is_real_z=cur_is_real_z[isCandidate];

                  /*
                    find the index with lowest weights but take the
                    highest index if multiple obs exist!
                    This helps to prevent suppressing marginal cells
                  */
                  double minVal=min(cur_weights);
                  IntegerVector v=seq(0, cur_weights.size()-1);
                  IntegerVector tmpres=v[cur_weights==minVal];
                  ind_x=max(tmpres);
                  int finalIndex=cur_id[ind_x]; // c-indices
                  int final_stIndex=cur_sub_ids[ind_x];

                  /* check if suppressed cell is > 0! */
                  if (freq[finalIndex]==0) {
                    stop("frequency of suppressed cell is 0!\n");
                  }

                  /* print info about newly suppressed cell */
                  if (debug==true) {
                    Rcout << "--> Suppression found" << std::endl;
                    DataFrame dfSupp = DataFrame::create(
                      Named("id")=finalIndex,
                      Named("freq")=freq[finalIndex],
                      Named("weights")=weights[finalIndex],
                      Named("sdcStatus_num")=sdcStatus_num[finalIndex]);
                    print(dfSupp);
                    R_FlushConsole();
                  }
                  additional_supps=additional_supps+1;

                  /* add suppressions to required vectors (also to current subtable!) */
                  sdcStatus[finalIndex]="x";
                  sdcStatus_num[finalIndex]=2;
                  st_sdcStatus_num[final_stIndex]=2;
                  suppsAdded=suppsAdded+1;
                } else {
                  // the single suppressed cell is protected by a dummy-cell that
                  // must never be published!
                }
              } else {
                // nothing to do for subtables with only one cell
              }
            } /* end j-loop (simple tables) */
          } /* end i-loop for groups */

          if (additional_supps == 0) {
            /* we can stop the inner loop and move to the next subtable */
            if (verbose == true) {
              Rcout << "Check of group " << group+1 <<"|" << nrGroups << ": ";
              Rcout << "subTable " << tab+1 <<"|" << nrTabs;
              Rcout << " | nrCells: " << n_st;
              if (suppsAdded == 0) {
                Rcout << " | everything ok/nothing todo." << std::endl;
              } else {
                Rcout << " | additionally suppressed cells: " << suppsAdded << std::endl;
              }
              R_FlushConsole();
            }
            newSuppsAdded=false;
            total_new_supps=total_new_supps+suppsAdded;
          } else {
            final_ok=false;
          }
        } /* inner while()-loop */
      } /* end for-loop (tabs) */

      /* set all 's' cells in this subtable to 'z' */
      for (int x1=0; x1<st_freq.size(); x1++) {
        if ((st_sdcStatus_num[x1]==0) and st_freq[x1] > 0) {
          int ind=st_id[x1];
          sdcStatus_num[ind]=3;
          sdcStatus[ind]="z";
        }
      }
    } /* end for-loop (groups) */

    /*
      check if we went through all subtables in all groups
      without setting any 'z' cells to 's'.
      In this case, we can stop the outer while-loop!
    */
    if ((override==false) & (final_ok==true)) {
      runInd=false;
    } else {
      counter=counter+1;
      /* we need to set all 'z' cells back to to 's' */
      LogicalVector vv1=((sdcStatus_num==0) | (sdcStatus_num==3)) & (is_real_z==false);
      sdcStatus_num[vv1]=0;
      sdcStatus[vv1]="s";

      LogicalVector vv2=((sdcStatus_num==0) | (sdcStatus_num==3)) & (is_real_z==false);
      sdcStatus_num[vv2]=3;
      sdcStatus[vv2]="z";
    }
  }

  if (verbose==true) {
    Rcout << "Finished. Total number of new suppressions: " << total_new_supps << std::endl;
  }
  IntegerVector total_new_suppsv(1);
  total_new_suppsv[0]=total_new_supps;
  return Rcpp::List::create(
    Rcpp::Named("id")=id,
    Rcpp::Named("freq")=freq,
    Rcpp::Named("sdcStatus")=sdcStatus,
    Rcpp::Named("total_new_supps")=total_new_suppsv);
}
