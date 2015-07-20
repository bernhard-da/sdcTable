#include <Rcpp.h>
using namespace Rcpp;

/*
 * codes in sdcStatus:
 * 0: possible suppression partner ("s")
 * 1: primary suppression ("u")
 * 2: secondary suppression ("x")
 * 3: non-applicable ("z")
 */
IntegerVector convert_sdcStatus_to_num(CharacterVector sdcStatus) {
  int n=sdcStatus.size();
  IntegerVector result(n);
  for ( int i=0; i<n; i++ ) {
    if ( sdcStatus[i]=="s" ) {
      result[i] = 0;
    }
    if ( sdcStatus[i]=="u" ) {
      result[i] = 1;
    }
    if ( sdcStatus[i]=="x" ) {
      result[i] = 2;
    }
    if ( sdcStatus[i]=="z" ) {
      result[i] = 3;
    }
  }
  return(result);
}

// [[Rcpp::export]]
List greedyMultDimSuppression(DataFrame dat, List indices, List subIndices, IntegerVector dimVars, bool verbose) {
  bool debug=false;
  Function cpp_print("print");

  int nrDims=dimVars.size();
  Rcout << "we are going to protect an " << nrDims << " dimensional dataset!" << std::endl;

  List subList;

  /* extracting input data from list */
  IntegerVector freq=dat["freq"];
  CharacterVector sdcStatus=dat["sdcStatus"];
  IntegerVector sdcStatus_num=convert_sdcStatus_to_num(sdcStatus);
  IntegerVector id=dat["id"];

  IntegerVector current_indices;

  IntegerVector st_freq;
  CharacterVector st_sdcStatus;
  IntegerVector st_sdcStatus_num;
  IntegerVector st_id;
  List st_subIndices;
  int nrGroups=indices.size();
  IntegerVector ind_x(1);
  bool runInd=true;
  int total_new_supps=0; /* total number of required secondary supps */
  while (runInd==true ) {
    bool override=false;
    for ( int group=0; group<nrGroups; group++ ) {
      List subList=indices[group];
      int nrTabs=subList.size();
      /* tab is iterating over all subtables within a given group */
      for ( int tab=0; tab<nrTabs; tab++ ) {
        /* subsetting data to current subtable */
        List subList=indices[group];
        current_indices=subList[tab];
        int n_st=current_indices.size();
        //if ( verbose ) {
        //  Rcout << "Group " << group+1 <<"|" << nrGroups;
        //  Rcout << " | subTable " << tab+1 <<"|" << nrTabs;
        //  Rcout << " | nrElements: " << n_st;
        //  R_FlushConsole();
        //}

        /* Select freqs and suppression pattern of current subtable */
        st_freq=freq[current_indices-1];
        st_sdcStatus=sdcStatus[current_indices-1];
        st_sdcStatus_num=sdcStatus_num[current_indices-1];
        st_id=id[current_indices-1]; /* ids in entire dataset, c-style (starting with 0!) */
        st_id=st_id-1;

        IntegerVector sub_ids=seq_along(st_id);
        sub_ids=sub_ids-1;

        /* create list for indices defining subtables */
        List tmp_subIndices=subIndices[group];
        st_subIndices=tmp_subIndices[tab];

        /* iterativly protecting the current simple tab */
        bool newSuppsAdded=true;
        int suppsAdded=0;
        while ( newSuppsAdded==true ) {
          int additional_supps=0;
          for ( int i=0; i<nrDims; i++ ) {
            IntegerVector gr=st_subIndices[i];
            int nrSimpleTables=max(gr);
            IntegerVector c_indices=seq_along(gr)-1; // indices relative to current simple table, c-style!

            for ( int j=0; j<nrSimpleTables; j++) {
              /* subset original vectors to current subtable */
              /* ind_c: index relative to entire data, starting with 0 */
              IntegerVector ind_c=c_indices[gr==j+1];

              /* extract values of current simple table */
              IntegerVector cur_freq=st_freq[ind_c];
              IntegerVector cur_freq_o=clone(cur_freq);
              IntegerVector cur_sdcStatus=st_sdcStatus_num[ind_c];
              IntegerVector cur_id=st_id[ind_c];
              IntegerVector cur_sub_ids=sub_ids[ind_c];

              /* calculating the number of suppressed cells */
              int nr_supps=0;
              int nCells=cur_freq.size();
              int upVal=max(cur_freq)+1;
              for ( int kk=0; kk<nCells; kk++ ) {
                if ( (cur_sdcStatus[kk] == 1) or (cur_sdcStatus[kk]==2) ) {
                  nr_supps = nr_supps+1;
                  cur_freq[kk]=upVal;
                }
                if ( cur_sdcStatus[kk] == 3 ) {
                  cur_freq[kk]=upVal;
                }
              }

              /*
                we have at least 1 observation in the current
                simple table but only a single suppressed cell
              */
              ind_x[0]=-1;
              if ( (nCells > 1) & (nr_supps == 1) ) {
                /* show output on current simple table */
                if ( debug == true ) {
                  Rcout << "--> one suppression in simpleTable " << j+1 << "|" << nrSimpleTables << " and dim " << i+1 << std::endl;
                  for ( int m=0; m<nCells; m++) {
                    Rcout << "--> index: " << cur_id[m];
                    Rcout << " (freq:" << cur_freq_o[m];
                    Rcout << " | status: " << cur_sdcStatus[m];
                    Rcout << " | sub_ids: " << cur_sub_ids[m] << ")" << std::endl;
                    R_FlushConsole();
                  }
                }

                ind_x=which_min(cur_freq);
                if ( cur_sdcStatus[ind_x[0]]!=0 ) {
                  /*
                    In this case, it is not possible to find a suppression
                    pattern with only 's'-cells.
                    we need to relax 'z' (code 3) cells to 's' (code 0) and try again
                  */
                  override=true;
                  IntegerVector xx=seq_along(cur_sdcStatus);
                  IntegerVector ind_x2=xx[(cur_sdcStatus==3) & (cur_freq>0)];
                  if ( ind_x2.size()==0 ) {
                    stop("cannot find suitable suppression pattern! (no cells with status=3 > freq>0 available!)\n");
                  }
                  int ii=0;
                  for ( int z=0; z<ind_x2.size(); z++ ) {
                    ii = ind_x2[z]-1;
                    if ( cur_freq_o[ii] > 0 ) {
                      cur_freq[ii] = cur_freq_o[ii]; // reset st_freq to original value
                      cur_sdcStatus[ii] = 0; // set cur_sdcStatus temporarily to 's' (0)
                    }

                    /*
                    if ( debug ) {
                      Rcout << "ii: " << ii;
                      Rcout << " | cur_freq: " << cur_freq[ii];
                      Rcout << " | cur_freq_o: " << cur_freq_o[ii];
                      Rcout << " | cur_sdcStatus: " << cur_sdcStatus[ii];
                      Rcout << " | sub_ids: " << cur_sub_ids[ii] << std::endl;
                      R_FlushConsole();
                    }
                    */
                  }
                  ind_x = which_min(cur_freq);
                  if ( ind_x.size()==0 ) {
                    stop("Something went horribly wrong!\n");
                  }
                }
              }
              if ( ind_x[0]!=-1 ) {
                int finalIndex=cur_id[ind_x[0]]; // c-indices
                int final_stIndex=cur_sub_ids[ind_x[0]];

                /* print info about newly suppressed cell */
                if ( debug ) {
                  Rcout << "--> Suppression found | x_ind: " << ind_x[0];
                  Rcout << " | freq: " << freq[finalIndex];
                  Rcout << " | finalIndex: " << finalIndex;
                  Rcout << " | sdcStatus_num: " << sdcStatus_num[finalIndex];
                  Rcout << " | final_stIndex: " << final_stIndex << std::endl;
                }
                additional_supps=additional_supps+1;

                /* check if suppressed cell is > 0! */
                if ( freq[finalIndex]==0 ) {
                  stop("frequency of suppressed cell is 0!\n");
                }

                /* add suppressions to required vectors (also to current subtable!) */
                sdcStatus[finalIndex]="x";
                sdcStatus_num[finalIndex]=2;
                st_sdcStatus_num[final_stIndex]=2;
                suppsAdded=suppsAdded+1;
              }
            } /* end j-loop (simple tables) */
          } /* end i-loop for groups */
          if ( additional_supps == 0 ) {
            /* we can stop the inner loop and move to the next subtable */
            if ( (verbose==true) & (suppsAdded >0 ) ) {
              Rcout << "Group " << group+1 <<"|" << nrGroups;
              Rcout << " | subTable " << tab+1 <<"|" << nrTabs;
              Rcout << " | nrElements: " << n_st;
              Rcout << " | additional supp: " << suppsAdded << std::endl;
              R_FlushConsole();
            }
            newSuppsAdded=false;
            total_new_supps=total_new_supps+suppsAdded;
          }
        } /* inner while()-loop */
      } /* end for-loop (tabs) */

      /* set all 's' cells in this subtable to 'z' */
      for ( int x1=0; x1<st_freq.size(); x1++ ) {
        if ( (st_sdcStatus_num[x1]==0) and st_freq[x1] > 0 ) {
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
    if ( override==false ) {
      runInd=false;
    } else {
      /* we need to set all 'z' cells back to to 's' */
      for ( int i=0; i<freq.size(); i++) {
        if ( (sdcStatus_num[i]==0) or (sdcStatus_num[i]==3) ) {
          if ( freq[i]>0 ) {
            sdcStatus_num[i]=0;
            sdcStatus[i]="s";
          } else {
            sdcStatus_num[i]=3;
            sdcStatus[i]="z";
          }
        }
      }
    }
  }

  if ( verbose ) {
    Rcout << "Finished. Total Suppressions: " << total_new_supps << std::endl;
  }
  IntegerVector total_new_suppsv(1);
  total_new_suppsv[0]=total_new_supps;
  return Rcpp::List::create(
    Rcpp::Named("id")=id,
    Rcpp::Named("freq")=freq,
    Rcpp::Named("sdcStatus")=sdcStatus,
    Rcpp::Named("total_new_supps")=total_new_suppsv);
 }
