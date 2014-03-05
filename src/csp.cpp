#include <glpk.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <list>
#include <iterator>
#include <numeric>
using namespace std;

/* structs to get index or sorted value|index pair */
/* see also, http://bit.ly/1kfn8pR */
struct val_index_double { 
  double number;
  size_t index;
};
struct val_index_int { 
  int number;
  size_t index;
};
struct sort_by_number_double { 
  bool operator()(val_index_double const &left, val_index_double const &right) { 
    return left.number < right.number;
  }
};
struct sort_by_number_int { 
  bool operator()(val_index_int const &left, val_index_int const &right) { 
    return left.number < right.number;
  }
};
struct branchnode {
  vector<int> indices;
  vector<double> values;
};
struct mprob_constraint {
  vector<int> indices;
  vector<double> values;
  int type;
  double low;
  double up;
  bool is_active;
};

/* write the constraint pool in a file */
void write_constraint_pool(list<mprob_constraint>& constraint_pool) {
  ofstream f ("pool.txt");
  f << "we have a total of " << constraint_pool.size() << " constraints!\n\n";

  int nr = 1;
  for(std::list<mprob_constraint>::iterator iter=constraint_pool.begin(); iter!=constraint_pool.end(); ++iter) {
    f << "\nconstraint nr " << nr << ":\n";

    f << "index | values: ";
    for ( int i = 0; i < iter->indices.size(); ++i ) {
      f << iter->indices[i] << " (" << iter->values[i] << ") ";
    }
    f << "\n";
    nr = nr + 1;  
  }   
  f.close();  
}

void clean_up_constraints(glp_prob *mprob) {
  bool is_slack;
  vector<int> remove_indices;
  for( int i=1; i<=glp_get_num_rows(mprob); ++i ) {
    is_slack = false;
    /* 
      http://bit.ly/1nxC5Rm
      basic row <--> inactive constraint
    */    
    if ( glp_get_row_stat(mprob, i) == 1 ) {
      is_slack = true;
      remove_indices.push_back(i);
    }
  }

  int nr_elements = remove_indices.size();  
  if ( nr_elements > 0 ) {
    int remove_indices2[nr_elements+1];
    remove_indices2[0] = -1;
    for ( int k=0; k<nr_elements; ++k ) {
      remove_indices2[k+1] = remove_indices[k];
    }    
    glp_del_rows(mprob, remove_indices.size(), remove_indices2); 
    //printf("\n--> we have removed %d constraints from mprob!\n", nr_elements);
  }
}

/* we remove invalid constraints due to fixing variables! */
/*
void remove_invalid_constraints(glp_prob *mprob) {
   bool is_invalid;
   vector<int> remove_indices;
   double lhs;
   int inds[glp_get_num_cols(mprob)+1];
   double vals[glp_get_num_cols(mprob)+1];  
   int len;
   for( int i=1; i<=glp_get_num_rows(mprob); ++i ) {
    is_invalid = false;

    len = glp_get_mat_row(mprob, i, inds, vals);
    lhs = 0;
    for ( int k=1; k<=len; ++k ) {
      lhs += vals[k] * glp_get_col_prim(mprob, inds[k]);
    }

    if ( lhs < glp_get_row_lb(mprob, i) ) {
      printf("lhs=%g | lb=%g\n", lhs, glp_get_row_lb(mprob, i));

      is_invalid = true;
      remove_indices.push_back(i);
    }
  }

  int nr_elements = remove_indices.size();  
  if ( nr_elements > 0 ) {
    int remove_indices2[nr_elements+1];
    remove_indices2[0] = -1;
    for ( int k=0; k<nr_elements; ++k ) {
      remove_indices2[k+1] = remove_indices[k];
    }    
    glp_del_rows(mprob, remove_indices.size(), remove_indices2); 
    printf("\n--> we have removed %d constraints from mprob (bounding)!\n", nr_elements);
  }
}
*/

void insert_violated_constraints(glp_prob *mprob, list<mprob_constraint>& constraint_pool, vector<double> &xi) {
  /* 
     for each constraint check if opt. solution violates it
     if so, set it to active and add it to master problem   
  */
  bool is_violated;
  int nr_elements;
  int nr_added=0;
  double lhs;
  bool first=true;

  int run_ind = 0;
  for(std::list<mprob_constraint>::iterator iter=constraint_pool.begin(); iter!=constraint_pool.end(); ++iter) {
    /* calculate inner product */
    lhs = 0;
    for ( int i=0; i<iter->indices.size(); ++i ) {
      lhs = lhs + glp_get_col_prim(mprob, iter->indices[i]) * iter->values[i];
    }    

    is_violated = false;

    // GLP_LO
    if ( (iter->type == 2) and (lhs < iter->low) ) { 
      is_violated = true;

      //if ( first == true ) {
      //  printf("run_ind=%d\n", run_ind);
      //  write_constraint_pool(constraint_pool);
      //  first = false;
      //}
    }   
    /*
    // GLP_UP
    if ( iter->type == 3 && lhs > iter->up ) {
    is_violated = true;
    }
    // GLP_DB
    if ( iter->type == 4 ) {
    if ( lhs < iter->low || lhs > iter->up  ) {
    is_violated = true;
    }        
    }   
    // GLP_FX
    if ( iter->type == 5 && lhs != iter->up ) {
    is_violated = true;
    }         
    */
    /* current solution violates this constraint, we are adding it! */

    if ( is_violated == true ) {
      nr_added += 1;
      glp_add_rows(mprob, 1); // add one constraint
      int lastrow = glp_get_num_rows(mprob);

      /*
         if ( first == true ) {
         printf("run_ind=%d\n", run_ind);
         write_constraint_pool(constraint_pool);
         glp_write_lp(mprob, NULL, "mprob_debug1.txt"); 
         }
      */
      if ( iter->type == 2 ) {
        glp_set_row_bnds(mprob, lastrow, GLP_LO, iter->low, iter->up);
      }

      nr_elements = iter->indices.size();
      int indices2[nr_elements+1];
      double vals2[nr_elements+1];
      indices2[0] = -1;
      vals2[0] = -1;
      for ( int k=0; k<nr_elements; ++k ) {
        indices2[k+1] = iter->indices[k];
        vals2[k+1] = iter->values[k];
      }

      /*
      if ( first == true ) {
        for ( int k=0; k <nr_elements+1; ++k ) {
          if ( vals2[k] > 0 ) {
            printf("ind2=%d | vals2=%g\n", indices2[k], vals2[k]);
          }
        }
      }
      */
      glp_set_mat_row(mprob, lastrow, nr_elements, indices2, vals2);        
      iter->is_active == true;

      /*
      if ( first == true ) {
        printf("lastrow=%d\n", lastrow);
        glp_write_lp(mprob, NULL, "mprob_debug2.txt"); 
        first = false;
      } 
      */
    }      
    run_ind += 1;
  }
  //printf("--> we have added %d constraints to the master problem!\n", nr_added);
}

void update_constraint_pool(list<mprob_constraint>& constraint_pool, int *constraint_indices, double *constraint_values, double bound, int type, int nr_vars) {
  mprob_constraint con;

  vector<int> ind;
  vector<double> vals;
  for ( int i=1; i<=nr_vars; ++i) {
    if ( constraint_values[i] != 0 ) {
      ind.push_back(constraint_indices[i]);
      vals.push_back(constraint_values[i]);      
    }
  }

  if ( vals.size() > 1 ) {
    con.indices = ind;
    con.values = vals;
    con.low = bound;
    con.up = bound;
    con.is_active = true;
    con.type = type;
    constraint_pool.push_back(con);    
  }
}

/* structure to hold information on bounds and constraint matrix M of given sdcproblem instance */
struct sdcinfo {
  vector<int> vals;
  vector<int> LPL;
  vector<int> UPL;
  vector<int> SPL;
  vector<double> UB;
  vector<double> LB;
  int *ia;
  int *ja;
  double *ar;
  int *ind_prim;
  int len_prim;
  int *ind_fixed;
  int len_fixed;
  int cells_mat;
  int nr_vars;
  int nr_rows;
  double upper_bound;
  vector<int> current_best_solution;
  vector<int> branchvars;
  double tol;
  bool verbose;
};

bool is_integer(double n, double tol) {
  double absd = abs(n);
  if( absd - floor(absd) > 0.5 ) {
    return (ceil(absd) - absd) < tol;
  }      
  return (n - floor(absd)) < tol;
}

bool solution_is_integer(glp_prob *linprob, double tol) {
  bool result = true;
  for (	int i=1; i<=glp_get_num_cols(linprob); ++i ) {
    if ( is_integer(glp_get_col_prim(linprob, i), tol) == false ) {
      result = false;
      break;
    }    
  }
  return result;
}

void update_master_problem(glp_prob *mprob, int *constraint_indices, double *constraint_values, double bound) {
  glp_add_rows(mprob, 1); // add one constraint
  int lastrow = glp_get_num_rows(mprob);
  int nr_vars = glp_get_num_cols(mprob);
  glp_set_mat_row(mprob, lastrow, nr_vars, &constraint_indices[0], &constraint_values[0]);
  glp_set_row_bnds(mprob, lastrow, GLP_LO, bound, 0.0); // lower bound, ub ignored
}

void solve_master_problem(glp_prob *mprob, vector<double> &xi, sdcinfo *info) {  
  glp_simplex(mprob, NULL);
  clean_up_constraints(mprob); 
  /* updating the solution */
  for ( int j=0; j < info->nr_vars; ++j ) {  
    xi[j] = glp_get_col_prim(mprob, j+1); 
  }
}

/* add a constraint to a problem */
/*
void add_constraint(glp_prob *linprob, int type, double lb, double ub, int len, int *indices, double *vals) {
  glp_add_rows(linprob, 1); // add one constraint
  int lastrow = glp_get_num_rows(linprob);
  glp_set_row_bnds(linprob, lastrow, type, lb, ub);
  glp_set_mat_row(linprob, lastrow, len, indices, vals);
}
*/

glp_prob * setup_attacker_problem(sdcinfo *info, vector<double> &xi) {        
  /* we are setting the dual problem, we have to transpose the matrix */
  int rr = info->nr_vars;
  int cc = info->nr_rows;

  glp_prob *aprob;
  aprob = glp_create_prob();
  glp_set_prob_name(aprob, "attackers_problem2");

  glp_add_cols(aprob, cc);	
  glp_add_rows(aprob, rr);
  glp_load_matrix(aprob, info->cells_mat-1, info->ja, info->ia, info->ar);  

  /* row bounds eqal zero when initializing the problem */
  int bound, index;
  for ( int j=1; j<=rr; ++j ) {
    glp_set_row_bnds(aprob, j, GLP_FX, 0.0, 0.0);
  }

  for ( int j=1; j<=cc; ++j ) {
    if ( j <= 2*info->nr_vars ) {
      glp_set_col_bnds(aprob, j, GLP_LO, 0.0, 0.0);
    } else {
      glp_set_col_bnds(aprob, j, GLP_FR, 0.0, 0.0);
    }     
  }
  return aprob;
}

/*
aprob: pointer to attacker's problem
mprob: pointer to master problem
indexvar: cell that should be attacked
sdcinfo: sdcinfo object containing bounds,...
bridgeless=0 -> normal attackerProblem
bridgeless=1 -> calculate bridgeless inequalities
*/
int solve_att_prob(glp_prob *aprob, glp_prob *mprob, list<mprob_constraint>& constraint_pool, int indexvar, sdcinfo *info, vector<double> &xi, int bridgeless, bool verbose) {
  int o_LPL = info->LPL[indexvar-1];
  int o_UPL = info->UPL[indexvar-1];
  int o_SPL = info->SPL[indexvar-1];
  if ( bridgeless == 1 ) {
    info->LPL[indexvar-1] = 0;
    info->UPL[indexvar-1] = 0;
    info->SPL[indexvar-1] = 1;
  }

  int nr_additional_constraints = 0;

  /* update objective function */  
  int nr_vars = glp_get_num_cols(aprob);
  int nr_rows = glp_get_num_rows(aprob);

  int nr_real_variables = nr_rows;

  double obj_up=0; 
  double obj_low=0;

  for ( int j=1; j <= nr_real_variables; ++j ) {
    obj_up = info->vals[j-1] + xi[j-1]*info->UB[j-1];
    obj_low = (-1)*(info->vals[j-1] - xi[j-1]*info->LB[j-1]);       
    glp_set_obj_coef(aprob, j, obj_up);
    glp_set_obj_coef(aprob, j+nr_real_variables, obj_low);
  }    

  // define vectors for possible insert into master problem
  double constraint_values_min[nr_real_variables+1];
  double constraint_values_max[nr_real_variables+1];
  double constraint_values_tot[nr_real_variables+1];

  double alphas_min[nr_real_variables+1];
  double alphas_max[nr_real_variables+1];
  double betas_min[nr_real_variables+1];
  double betas_max[nr_real_variables+1]; 
  double alphas_tot[nr_real_variables+1]; 
  double betas_tot[nr_real_variables+1]; 

  int constraint_indices[nr_real_variables+1];
  double constraint_val, len, zmin, zmax;

  //int stat_up, stat_down;

  /* initialize constraints index vector (1:nr of variables) */
  for ( int i=0; i <= nr_real_variables; ++i ) {
    constraint_indices[i] = i;
  }
  glp_set_obj_dir(aprob, GLP_MIN);

  /* 1) maximize the problem */
  if ( info->UPL[indexvar-1] > 0 | info->SPL[indexvar-1] > 0 ) {
    /* set rowbound to 1 for the active cell cell */
    glp_set_row_bnds(aprob, indexvar, GLP_FX, 1, 1);        

    glp_simplex(aprob, NULL);		
    zmax = glp_get_obj_val(aprob);
    //stat_up = glp_get_status(aprob);
    //glp_write_lp(aprob, NULL, "aprob-max.txt");

    // adding a constraint to the master problem, if necessary
    if ( zmax < info->vals[indexvar-1] + info->UPL[indexvar-1] ) {
      constraint_val = (double)info->UPL[indexvar-1];
      for ( int j=1; j<=nr_real_variables; ++j) {
        alphas_max[j] = glp_get_col_prim(aprob, j);
        betas_max[j] = glp_get_col_prim(aprob, j+nr_real_variables);
        double v = glp_get_col_prim(aprob, j)*info->UB[j-1] + glp_get_col_prim(aprob,j+nr_real_variables)*info->LB[j-1];
        if ( abs(v) < info->tol ) { v = 0.0; }
        constraint_values_max[j] = fmin(v, constraint_val);
      }
      if ( bridgeless == 0 ) {
        update_master_problem(mprob, constraint_indices, constraint_values_max, constraint_val);
        update_constraint_pool(constraint_pool, constraint_indices, constraint_values_max, constraint_val, 2, info->nr_vars);
        nr_additional_constraints = nr_additional_constraints + 1;
      }
    }
  }

  /* 2) minimize the problem */
  if ( info->LPL[indexvar-1] > 0 | info->SPL[indexvar-1] > 0 ) {
    /* set rowbound to -1 for the active cell cell */
    glp_set_row_bnds(aprob, indexvar, GLP_FX, -1, -1);        
    glp_simplex(aprob, NULL);
    zmin = (-1)*glp_get_obj_val(aprob);
    //stat_down = glp_get_status(aprob);

    //glp_write_lp(aprob, NULL, "aprob-min.txt");  

    // adding a constraint to the master problem, if necessary
    if ( zmin > info->vals[indexvar-1] - info->LPL[indexvar-1] ) {
      constraint_val = (double)info->LPL[indexvar-1];
      for ( int j=1; j<=nr_real_variables; ++j) {
        alphas_min[j] = glp_get_col_prim(aprob, j);
        betas_min[j] = glp_get_col_prim(aprob, j+nr_real_variables);
        double v = glp_get_col_prim(aprob, j)*info->UB[j-1] + glp_get_col_prim(aprob,j+nr_real_variables)*info->LB[j-1];
        if ( abs(v) < info->tol ) { v = 0.0; }
        constraint_values_min[j] = fmin(v, constraint_val);
      }
      if ( bridgeless == 0 ) {
        update_master_problem(mprob, constraint_indices, constraint_values_min, constraint_val);
        update_constraint_pool(constraint_pool, constraint_indices, constraint_values_min, constraint_val, 2, info->nr_vars);
        nr_additional_constraints = nr_additional_constraints + 1;        
      }
    }      
  }        

  /* check sliding protection level  */
  len = zmax - zmin;
  if ( info->SPL[indexvar-1] > 0 && len < info->SPL[indexvar-1] ) {
    if ( bridgeless == 0 ) {
      constraint_val = (double)info->SPL[indexvar-1];
    }
    if ( bridgeless == 1 ) {
      constraint_val = 0;
    }

    for ( int j=1; j<=nr_real_variables; ++j) {
      alphas_tot[j] = alphas_min[j] + alphas_max[j];
      betas_tot[j] = betas_min[j] + betas_max[j];  
      double v = alphas_tot[j]*info->UB[j-1] + betas_tot[j]*info->LB[j-1];
      if ( abs(v) < info->tol ) { v = 0.0; }
      if ( bridgeless == 0 ) {
        constraint_values_tot[j] = fmin(v, constraint_val);
      }
      if ( bridgeless == 1 ) {
        if ( j == indexvar ) {
          constraint_values_tot[j] = -1;
        } else {
          if ( v > 0 ) {
            constraint_values_tot[j] = 1;
          } else {
            constraint_values_tot[j] = 0;
          }          
        }
      }      
    }
    update_master_problem(mprob, constraint_indices, constraint_values_tot, constraint_val);
    update_constraint_pool(constraint_pool, constraint_indices, constraint_values_tot, constraint_val, 2, info->nr_vars);
    nr_additional_constraints = nr_additional_constraints + 1;      
  }

  /* reset rowbound for cell that should be maxed/minimized */
  glp_set_row_bnds(aprob, indexvar, GLP_FX, 0, 0);

  /* reset original bounds for bridgeless vars */
  if ( bridgeless == 1 ) {
    info->LPL[indexvar-1] = o_LPL;
    info->UPL[indexvar-1] = o_UPL;
    info->SPL[indexvar-1] = o_SPL;
  }
  /*  
  if ( info->verbose == true && bridgeless == 0 ) {
    printf("prim_supp=%d (%d): [[ %g | %g ]]\n", indexvar, info->vals[indexvar-1], zmin, zmax);
  }
  */
  /*
  if ( info->verbose && bridgeless == 1 ) {
    printf("bridgevar=%d (%d): [[ %g | %g ]]\n", indexvar, info->vals[indexvar-1], zmin, zmax);
  }
  */
  return(nr_additional_constraints);
}

glp_prob * setup_incprob(sdcinfo *info, vector<double> &xi) {  
  glp_prob *incprob;
  incprob = glp_create_prob();
  glp_set_prob_name(incprob, "incprob");

  glp_set_obj_dir(incprob, GLP_MIN);

  int nr_vars = info->LB.size();
  glp_add_cols(incprob, 2*nr_vars);

  /* set objective coefficients */
  for ( int i=1; i <= nr_vars; ++i ) {
    glp_set_obj_coef(incprob, i, 0.0);

    glp_set_col_bnds(incprob, i, GLP_DB, 0.0, info->UB[i-1]);
    glp_set_col_bnds(incprob, i+nr_vars, GLP_DB, 0.0, info->LB[i-1]);
  }

  int nr_constraints = (info->cells_mat - 2*nr_vars -1);
  int ia2[1+nr_constraints*2];
  int ja2[1+nr_constraints*2];
  double ar2[1+nr_constraints*2];

  ia2[0] = 0;
  ja2[0] = 0;
  ar2[0] = 0;

  /* fixme: not optimal, can be done better, but still working */
  int ind, j;
  int nr_rows = 0;

  for (int i=1; i<=nr_constraints; ++i) {
    j = 2*nr_vars + i;
    ind = nr_constraints+i;
    ia2[i] = info->ia[j]-2*nr_vars; 
    ja2[i] = info->ja[j]; 
    ar2[i] = info->ar[j];    

    nr_rows = max(nr_rows, ia2[i]);

    ia2[ind] = info->ia[j]-2*nr_vars; 
    ja2[ind] = nr_vars+info->ja[j]; 
    ar2[ind] = (-1.0)*info->ar[j];

    nr_rows = max(nr_rows, ia2[ind]);
  }

  /* adding rows with bounds equal to 0 */
  glp_add_rows(incprob, info->nr_rows);
  for ( int i=1; i<=info->nr_rows; ++i ) {
    glp_set_row_bnds(incprob, i, GLP_FX, 0.0, 0.0);
  }  

  /* calculate length of array including leading 0 */
  int s_ia2 = (sizeof(ia2)/sizeof(*ia2));
  glp_load_matrix(incprob, s_ia2-1, ia2, ja2, ar2);
  //glp_write_lp(incprob, NULL, "incprob.txt");
  return incprob;
}

void heuristic_solution_primitive(sdcinfo *info) {
  double bound;  
  for ( int i=0; i < info->nr_vars; ++i ) {
    info->current_best_solution[i] = 1;
    bound += info->vals[i];
  }  
  info->upper_bound = bound;
}

void heuristic_solution(glp_prob *incprob, sdcinfo *info, vector<double> &xi, int use_existing_solution) {
  int ik, nr;
  int heuristic_solution[info->nr_vars];
  for ( int i=0; i<info->nr_vars; ++i) {
    heuristic_solution[i] = 0;
  }

  /* set obj coefficients */
  double v;
  for ( int j=1; j<=info->nr_vars; ++j ) {
    if ( use_existing_solution == 1 ) {
      v = (double)info->vals[j-1] * (1-xi[j-1]);
      glp_set_obj_coef(incprob, j, v);
      glp_set_obj_coef(incprob, info->nr_vars+j, v);
    } 
    if ( use_existing_solution == 0 ) {
      v = (double)info->vals[j-1];
      glp_set_obj_coef(incprob, j, v);
      glp_set_obj_coef(incprob, j+info->nr_vars, v);
    }
    
    if ( xi[j-1] == 1 ) {
      glp_set_obj_coef(incprob, j, 0.0);
      glp_set_obj_coef(incprob, j+info->nr_vars, 0.0);
    }
  } 
  
  int ind[3] = {0}; // initialize all elements to zero  
  double val[3] = {0.0,1.0,1.0};

  for ( int z=1; z<=info->len_prim; ++z ) {
    ik = info->ind_prim[z-1];
    if ( info->UPL[ik-1] > 0 || info->SPL[ik-1] > 0 ) {
      glp_set_col_bnds(incprob, ik, GLP_FX, info->UPL[ik-1], info->UPL[ik-1]);
      glp_set_col_bnds(incprob, ik+info->nr_vars, GLP_FX, 0.0, 0.0);
      glp_simplex(incprob, NULL); 
      for ( int j=1; j <= info->nr_vars; ++j ) {
        if ( glp_get_col_prim(incprob, j) + glp_get_col_prim(incprob, info->nr_vars+j) > 0 ) {
          heuristic_solution[j-1] = 1;
        } 
      }      
    }

    if ( info->LPL[ik-1] > 0 || info->SPL[ik-1] > 0 ) {
      glp_set_col_bnds(incprob, ik, GLP_FX, 0.0, 0.0);
      glp_set_col_bnds(incprob, ik+info->nr_vars, GLP_FX, info->LPL[ik-1], info->LPL[ik-1]);      
      glp_simplex(incprob, NULL); 
      for ( int j=1; j<=info->nr_vars; ++j ) {
        if ( glp_get_col_prim(incprob, j) + glp_get_col_prim(incprob, info->nr_vars+j) > 0 ) {
          heuristic_solution[j-1] = 1;
        } 
      }         
    }

    if ( info->SPL[ik-1] > 0 ) {
      glp_set_col_bnds(incprob, ik, GLP_DB, 0, info->UB[ik-1]);
      glp_set_col_bnds(incprob, ik+info->nr_vars, GLP_DB, 0, info->LB[ik-1]);
      nr = glp_add_rows(incprob, 1);

      ind[1] = ik;
      ind[2] = ik + info->nr_vars;

      glp_set_mat_row(incprob, nr, 2, ind, val);
      glp_set_row_bnds(incprob, nr, GLP_FX, info->SPL[ik-1], info->SPL[ik-1]); 
      glp_simplex(incprob, NULL);
      for ( int j=1; j<=info->nr_vars; ++j ) {
        if ( glp_get_col_prim(incprob, j) + glp_get_col_prim(incprob, info->nr_vars+j) > 0 ) {
          heuristic_solution[j-1] = 1;
        } 
      }     
    }
    
    /* set coefficients to zero for supressed cells */
    for ( int i=1; i<=info->nr_vars; ++i ) {
      if ( heuristic_solution[i-1] == 1 ) {
        glp_set_obj_coef(incprob, i, 0);
        glp_set_obj_coef(incprob, i+info->nr_vars, 0);
      }
    }
  }

  /* reset obj coefficients */
  for ( int j=1; j<=glp_get_num_cols(incprob); ++j ) {
    glp_set_obj_coef(incprob, j, 0.0);
  }

  /* calculate upper bound */
  double bound = 0.0;  
  for ( int j=1; j<=info->nr_vars; ++j ) {
    bound += (double)heuristic_solution[j-1] * (double)info->vals[j-1];
  }  
  if ( use_existing_solution == 0 ) {
    info->upper_bound = bound;
    for ( int i=0; i < info->nr_vars; ++i ) {
      info->current_best_solution[i] = heuristic_solution[i];
    }       
  }
  if ( use_existing_solution == 1 ) {
    if ( bound < info->upper_bound ) {
      if ( info->verbose == true ) {
        printf("improved heuristic solution was found: bound=%f!\n", bound);
      }
      for ( int i=0; i < info->nr_vars; ++i ) {
        info->current_best_solution[i] = heuristic_solution[i];
      }
      info->upper_bound = bound;
    }      
  }
}

void preprocess(glp_prob *aprob, glp_prob *mprob, sdcinfo *info, vector<double> &xi) {
  int nr_vars = glp_get_num_cols(aprob);  
  int nr_real_variables = info->nr_vars;
  vector<double> HIGH(info->len_prim);
  vector<double> LOW(info->len_prim);
  vector<int> vMAX(info->len_prim);
  vector<int>ind_sorted(info->len_prim);

  vector<val_index_int> index_pair;

  /* 1) initialize HIGH and LOW */  
  for ( int i=0; i < info->len_prim; ++i ) {    

    HIGH[i] = info->vals[info->ind_prim[i]-1];
    LOW[i] = info->vals[info->ind_prim[i]-1];

    /* instead of vMax and ind_sorted */
    val_index_int pair;
    pair.index = i;
    pair.number = max(info->SPL[info->ind_prim[i]-1], info->UPL[info->ind_prim[i]-1]+info->LPL[info->ind_prim[i]-1]);;
    index_pair.push_back(pair);
  }

  /* 2) order HIGH | LOW according to vMAX */
  sort(index_pair.begin(), index_pair.end(), sort_by_number_int());

  /* 3) setup attacker problems */
  /* update objective function */  
  double obj_up=0; 
  double obj_low=0;

  for ( int j=1; j <= nr_real_variables; ++j ) {
    obj_up = info->vals[j-1] + xi[j-1]*info->UB[j-1];
    obj_low = (-1)*(info->vals[j-1] - xi[j-1]*info->LB[j-1]);
    glp_set_obj_coef(aprob, j, obj_up);
    glp_set_obj_coef(aprob, j+nr_real_variables, obj_low);
  }       

  // define vectors for possible insert into master problem
  double constraint_values_min[nr_real_variables+1];
  double constraint_values_max[nr_real_variables+1];
  double constraint_values_tot[nr_real_variables+1];  

  int constraint_indices[nr_real_variables+1];
  double constraint_val, len, zmin, zmax;

  /* initialize constraints index vector (1:nr of variables) */
  for ( int i=0; i <= nr_real_variables; ++i ) {
    constraint_indices[i] = i;
  }
  glp_set_obj_dir(aprob, GLP_MIN);	

  int indexvar, arrval, sortvar;
  double v;
  for ( int i=0; i<info->len_prim; ++i ) {
    sortvar = index_pair[i].index;
    indexvar = info->ind_prim[sortvar];
    arrval = indexvar-1;    

    /* 1) maximize the problem */    
    if ( HIGH[sortvar] < info->vals[arrval] + info->UPL[arrval] || HIGH[sortvar] - LOW[sortvar] < info->SPL[arrval] ) {
      glp_set_row_bnds(aprob, indexvar, GLP_FX, 1, 1);  
      glp_simplex(aprob, NULL);		
      zmax = glp_get_obj_val(aprob);
      // adding a constraint to the master problem, if necessary      
      if ( zmax < info->vals[indexvar-1] + info->UPL[indexvar-1] ) {
        constraint_val = (double)info->UPL[indexvar-1];
        for ( int j=1; j<=nr_real_variables; ++j) {
          v = glp_get_col_prim(aprob, j)*info->UB[j-1] + glp_get_col_prim(aprob,j+nr_real_variables)*info->LB[j-1];
          if ( abs(v) < 1e-9 ) { v = 0.0; }
          constraint_values_max[j] = fmin(v, constraint_val);
        }
        update_master_problem(mprob, constraint_indices, constraint_values_max, constraint_val);
      }
      /* update arrays */
      HIGH[sortvar] = max(zmax, HIGH[sortvar]);
      LOW[sortvar] = min(zmax, HIGH[sortvar]);
    }

    /* 2) minimize the problem */
    if ( LOW[sortvar] > info->vals[arrval] - info->LPL[arrval] || HIGH[sortvar] - LOW[sortvar] < info->SPL[arrval] ) {
      glp_set_row_bnds(aprob, indexvar, GLP_FX, -1, -1);        
      glp_simplex(aprob, NULL);
      zmin = (-1)*glp_get_obj_val(aprob);      
      // adding a constraint to the master problem, if necessary
      if ( zmin > info->vals[indexvar-1] - info->LPL[indexvar-1] ) {        
        constraint_val = (double)info->LPL[indexvar-1];
        for ( int j=1; j<=nr_real_variables; ++j) {
          v = glp_get_col_prim(aprob, j)*info->UB[j-1] + glp_get_col_prim(aprob,j+nr_real_variables)*info->LB[j-1];
          if ( abs(v) < 1e-9 ) { v = 0.0; }          
          constraint_values_min[j] = fmin(v, constraint_val);
        }
        update_master_problem(mprob, constraint_indices, constraint_values_min, constraint_val);
      }
      /* update arrays */
      HIGH[sortvar] = fmax(zmin, HIGH[sortvar]);
      LOW[sortvar] = fmin(zmin, HIGH[sortvar]);
    }    
    /* reset rowbound for cell that should be maxed/minimized */
    glp_set_row_bnds(aprob, indexvar, GLP_FX, 0, 0);    
  }

  /* set UPL | LPL | SPL to 0, if possible depending on HIGH|LOW */
  int index;
  for ( int i=0; i < info->len_prim; ++i ) {
    index = info->ind_prim[i]-1;
    if ( HIGH[i] > info->vals[index] + info->UPL[index] ) {
      info->UPL[index] = 0;
    }
    if ( LOW[i] < info->vals[index] - info->LPL[index] ) {
      info->LPL[index] = 0;
    }  
    if ( HIGH[i] - LOW[i] >= info->SPL[index] ) {
      info->SPL[index] = 0;
    }      
  }
}

int calculate_branching_variable(glp_prob *mprob, vector<double> &xi, sdcinfo *info) {  
  vector<val_index_double> index_pair;  
  bool is_present;
  for ( int i=0; i < glp_get_num_cols(mprob); ++i ) {

    is_present = (std::find(info->branchvars.begin(), info->branchvars.end(), i+1) != info->branchvars.end());
    if ( is_present == false && is_integer(xi[i], info->tol) == false ) { /* new */
      val_index_double pair;
      pair.index = i+1;
      pair.number = (double)abs(xi[i]-0.5);
      index_pair.push_back(pair);   
    }
  }

  if ( index_pair.size() == 0 ) {
    return(0);
  }

  sort(index_pair.begin(), index_pair.end(), sort_by_number_double());
  int tmp_max = fmin(index_pair.size(), 10);

  double bmin, bmax;
  double totmax=0.0;
  int bv;
  for ( int i=0; i < tmp_max; ++i ) {
    // set branching variable to 0
    glp_set_col_bnds(mprob, index_pair[i].index, GLP_FX, 0, 0);
    glp_simplex(mprob, NULL);
    bmin = glp_get_obj_val(mprob);    

    // set branching variable to 1
    glp_set_col_bnds(mprob, index_pair[i].index, GLP_FX, 1, 1);
    glp_simplex(mprob, NULL);
    bmax = glp_get_obj_val(mprob);    

    // delete last constraint
    glp_set_col_bnds(mprob, index_pair[i].index, GLP_DB, 0, 1);    

    if ( i == 0 ) {
      bv = index_pair[i].index;
      totmax = (bmin+bmax)/2;
    } else {
      if ( (bmin+bmax)/2 > totmax ) {
        totmax = (bmin+bmax)/2;
        bv = index_pair[i].index;
      }      
    }
  }  
  return(bv);
}

/* delete all (row)-constraints from a problem */
void delete_all_constraints(glp_prob *p) {
  int xx[glp_get_num_rows(p)];
  int nrs = glp_get_num_rows(p);
  for ( int i=1; i<=nrs; ++i ) {
    xx[i] = i;
  }
  glp_del_rows(p, nrs, xx);
}

bool solve_relaxation(glp_prob *mprob, glp_prob *aprob, list<mprob_constraint>& constraint_pool, sdcinfo *info, vector<double> &xi) {
  //char mprob_name[1000]; char attack_name[1000];
  int run_ind=0;
  int nr_additional_constraints1, nr_additional_constraints2;
  double lower_bound;

  /* [TESTING]: remove all constraints from mprob */
  delete_all_constraints(mprob);

  //int solution_status;
  do {
    run_ind += 1;
    nr_additional_constraints1 = 0;
    nr_additional_constraints2 = 0;

    /* xi is updated by reference in the function! */
    solve_master_problem(mprob, xi, info);
    if ( glp_get_prim_stat(mprob) != GLP_FEAS ) {  
      return(false);
    }
    lower_bound = glp_get_obj_val(mprob);
    if ( lower_bound > info->upper_bound ) {
      break;
    }        

    /* fix variables with large reduced costs! */
    int nr_removed_vars = 0;
    for ( int j=1; j <= info->nr_vars; ++j ) {
      if ( glp_get_col_type(mprob, j) != GLP_FX ) {
        if ( glp_get_col_dual(mprob, j) > (info->upper_bound - lower_bound) ) {
          //printf("--> removing var %d: dc[%d]=%g > %f\n", j, j, glp_get_col_dual(mprob, j), (upper_bound - lower_bound));
          glp_set_col_bnds(mprob, j, GLP_FX, 0, 0);
          nr_removed_vars += 1;
        }	        
      }
    }
    /*
    if ( info->verbose == true ) {
      printf("%d variables were removed due to reduced costs!\n", nr_removed_vars);
    }
    */
    //sprintf(mprob_name, "master_problem_%d.txt", run_ind);
    //glp_write_lp(mprob, NULL, mprob_name);

    /* solve attacker problems for all primary sensitive cells */
    int bridgeless=0;  
    for ( int k=0; k < info->len_prim; ++k ) {
      nr_additional_constraints1 += solve_att_prob(aprob, mprob, constraint_pool, info->ind_prim[k], info, xi, bridgeless, false);

      //sprintf(attack_name, "attackers_problem_%d.txt", run_ind);
      //glp_write_lp(aprob, NULL, attack_name);
    }
    //printf("----------------------\n");

    /* generating bridgeless inequalities */
    if ( nr_additional_constraints1 == 0 ) {
      /* remove non-binding constraints */
      clean_up_constraints(mprob);

      bridgeless = 1;
      for ( int k=1; k <= info->nr_vars; ++k ) {
        double tmp = xi[k-1];

        if ( abs(glp_get_col_prim(mprob, k)) > info->tol ) {
          nr_additional_constraints2 += solve_att_prob(aprob, mprob, constraint_pool, k, info, xi, bridgeless, false);
        }
      }
      //printf("----------------------\n");
    }
    //printf("add_const1=%d | add_const2=%d\n", nr_additional_constraints1,nr_additional_constraints2);
  }
  while ( nr_additional_constraints1 + nr_additional_constraints2 > 0 );  

  /* is it a valid integer solution? */ 
  bool is_integer = solution_is_integer(mprob, info->tol);

  /* finished and integer, print solution */  
  int sol_val;
  if ( is_integer == true && glp_get_obj_val(mprob) < info->upper_bound ) {
    info->upper_bound = glp_get_obj_val(mprob);
    for ( int i=1; i <= info->nr_vars; ++i ) {
      sol_val = lround(glp_get_col_prim(mprob, i));
      xi[i-1] = sol_val;
      info->current_best_solution[i-1] = sol_val;      
    }
  }
  return(is_integer);  
}

void branch_and_bound(glp_prob *mprob, glp_prob *aprob, glp_prob *incprob, list<mprob_constraint>& constraint_pool, sdcinfo *info, vector<double> &xi) {
  vector<int> cur_indices;
  vector<double> cur_vals;
  bool is_int;
  int nr_additional_constraints;

  if ( info->verbose == true ) {
    printf("\nbranch and bound algorithm is starting!\n");
    printf("--> starting to calculate branching variable.\n");
  }

  /* calculate first branching variable (with index startig from 1!) */
  int bvar = calculate_branching_variable(mprob, xi, info); 

  /* define a queue to my branch_objects */ 
  if ( info->verbose == true ) {
    printf("--> initializing pool with first branchvar=%d...", bvar);
  }
  list<branchnode> pool;
  if ( info->verbose == true ) {
    printf("[done]\n");
  }
  branchnode n1, n2;
  cur_indices.push_back(bvar);
  cur_vals.push_back(0);

  n1.indices = cur_indices;
  n1.values = cur_vals;
  pool.push_back(n1);
  cur_indices.clear(); cur_vals.clear();
  cur_indices.push_back(bvar);
  cur_vals.push_back(1);  
  n2.indices = cur_indices;
  n2.values = cur_vals;  
  pool.push_back(n2);  
  int stop_ind = 0;
  while ( pool.empty() == false && stop_ind != 1 ) {        
    /* take first node and apply bounds to master problem */
    //branchnode cur_node = pool.front();
    branchnode cur_node = pool.back();

    //remove_invalid_constraints(mprob);

    //printf("... current indices that will be fixed ... \n");
    //for ( int i=0; i<cur_node.indices.size(); ++i ) {
    //  printf("fixing variable %d to %g!\n", cur_node.indices[i], cur_node.values[i]);
    //}   
    
    if ( info->verbose == true ) {
      printf("current node has %d elements!\n", cur_node.indices.size());
    }
    for ( int i=1; i<=info->nr_vars; ++i ) {
      glp_set_col_bnds(mprob, i, GLP_DB, 0, 1);
    }       

    /* reset bounds */
    if ( info->verbose == true ) {
      printf("fixing %d variables in the master problem with %d constraints ...\n", cur_node.indices.size(), glp_get_num_rows(mprob));
    }
    for ( int i=0; i<cur_node.indices.size(); ++i ) {
      glp_set_col_bnds(mprob, cur_node.indices[i], GLP_FX, cur_node.values[i], cur_node.values[i]);
    }
    for ( int i=0; i < info->len_prim; ++i ) {
      glp_set_col_bnds(mprob, info->ind_prim[i], GLP_FX, 1, 1);
    } 
    if ( info->len_fixed > 0 ) {
      for ( int i=0; i < info->len_fixed; ++i ) {
        glp_set_col_bnds(mprob, info->ind_fixed[i], GLP_FX, 0, 0);
      }
    }
    
    insert_violated_constraints(mprob, constraint_pool, xi);

    /* solve relaxed master problem with added constraints */
    is_int = solve_relaxation(mprob, aprob, constraint_pool, info, xi); 

    /* reset bounds in mprob */        
    for ( int i=0; i<cur_node.indices.size(); ++i ) {
      glp_set_col_bnds(mprob, cur_node.indices[i], GLP_DB, 0, 1);
    }        

    /* constraint added is not valid -> we have to remove the node from the pool */
    if ( glp_get_status(mprob) == GLP_NOFEAS ) {
      if ( info->verbose == true ) {
        printf("current node is not solvable. we are removing it!\n");
      }
      pool.pop_back();      
    } else {
      if ( glp_get_obj_val(mprob) > info->upper_bound ) {
        if ( info->verbose == true ) {
          printf("current node has obj > current best solution! --> pruning!\n");
        }
        pool.pop_back();
      } else {
        if ( is_int == true ) {
          if ( info->verbose == true ) {
            printf("we are finished and have a new valid best solution!\n"); 
          }
          // info->upper_bound and info->current_best_solution are updated in solve_relaxation()
          pool.pop_back();          
        } else {          
          pool.pop_back();        
          bvar = calculate_branching_variable(mprob, xi, info); 
          if ( bvar == 0 ) {
            if ( info->verbose == true ) {
              printf("no more branching possible!\n");
            }
          } else {
            if ( info->verbose == true ) {
              printf("we got another non-integer solution and ");
              printf("are branching again on variable %d!\n", bvar);               
            }
            
            /* trying to improve current best solution using our heuristics */
            if ( info->verbose == true ) {
              printf("trying to improve solution with our heuristic approach ...");
            }
            heuristic_solution(incprob, info, xi, 1);
            if ( info->verbose == true ) {
              printf("[done] (upper_bound=%g)\n", info->upper_bound);
            }
            /* add constraints to pool */        
            branchnode a1, a2;

            cur_indices.clear(); cur_vals.clear();
            cur_indices = cur_node.indices;        
            cur_vals = cur_node.values;      
            cur_indices.push_back(bvar);
            cur_vals.push_back(0);   
            a1.indices = cur_indices;
            a1.values = cur_vals;
            pool.push_back(a1);        

            cur_indices.clear(); cur_vals.clear();
            cur_indices = cur_node.indices;        
            cur_vals = cur_node.values;           
            cur_indices.push_back(bvar);
            cur_vals.push_back(1);  
            a2.indices = cur_indices;
            a2.values = cur_vals;  
            pool.push_back(a2);            
          }
        }
      }       
      /* removing current node from pool */
      if ( info->verbose == true ) {
        printf("poolsize=%d | upper_bound=%g | nr_constraints in pool=%d\n", pool.size(), info->upper_bound, constraint_pool.size());
      }
      //if ( pool.size() >= 200 ) {
      //  stop_ind = 1;
      //  //glp_write_lp(aprob, NULL, "latest_attacker_problem.txt");
      //}      
    }
  }
}

bool is_valid_solution(glp_prob *aprob, glp_prob *mprob, list<mprob_constraint>& constraint_pool, sdcinfo *info, vector<double> &xi) {
  //if ( info->verbose == true ) {
  //  printf("checking if current best solution is valid! --> ");
  //}
  bool ok = false;
  int nr_additional_constraints = 0;

  /* delete all constraints from mprob */
  int del_rows[glp_get_num_rows(mprob)];
  int nrs = glp_get_num_rows(mprob);
  for ( int i=1; i<=nrs; ++i ) {
    del_rows[i] = i;
  }
  glp_del_rows(mprob, nrs, del_rows);  

  for ( int k=0; k < info->len_prim; ++k ) {
    nr_additional_constraints += solve_att_prob(aprob, mprob, constraint_pool, info->ind_prim[k], info, xi, 0, false);
  }
  int bridgeless = 1;
  for ( int k=1; k <= info->nr_vars; ++k ) {
    double tmp = xi[k-1];
    if ( glp_get_col_type(mprob, k) != GLP_FX && abs(tmp) > 1e-9 ) {
      nr_additional_constraints += solve_att_prob(aprob, mprob, constraint_pool, k, info, xi, bridgeless, false);
    }
  }

  if ( nr_additional_constraints== 0 ) {
    ok = true;
  }
  //glp_write_lp(aprob, NULL, "aprob_validation.txt");
  return(ok);
}

extern "C" { 
  void csp(int *ind_prim, int *len_prim, int *ind_fixed, int *len_fixed, int *ia, int *ja, double *ar, 
      int *cells_mat, int *nr_vars, int *nr_rows, int *vals, double *lb, double *ub,
      int *LPL, int *UPL, int *SPL, int *final_pattern, int *verbose) {

    glp_term_out(GLP_OFF);

    /* constraint pool */
    list<mprob_constraint> constraint_pool;

    int bridgeless, use_existing_solution;
    bool is_int;

    /* converting to std::vectors */
    vector<int> v_ai(nr_vars[0]);
    vector<int> v_LPL(nr_vars[0]);
    vector<int> v_UPL(nr_vars[0]);
    vector<int> v_SPL(nr_vars[0]);
    vector<double> v_LB(nr_vars[0]);
    vector<double> v_UB(nr_vars[0]);

    v_ai.assign(vals, vals + nr_vars[0]);
    v_LPL.assign(LPL, LPL + nr_vars[0]);
    v_UPL.assign(UPL, UPL + nr_vars[0]);
    v_SPL.assign(SPL, SPL + nr_vars[0]);

    /* calculate relative external bounds LB, UB from lb, ub */
    for ( int i=0; i < nr_vars[0]; ++i ) {
      v_LB[i] = vals[i] - lb[i];			
      v_UB[i] = ub[i] - vals[i];
    }
    double upper_bound, lower_bound;

    /* 
       solution vectors for heuristic and 
       possible fractional solution of master problem 
    */
    vector<double> xi(nr_vars[0]);
    for ( int k=0; k<len_prim[0]; ++k ) {
      xi[k] = 1;
    }
    vector<int> xi_heur(nr_vars[0]);  

    /* set up master problem */
    glp_prob *mprob;
    mprob = glp_create_prob();
    glp_set_prob_name(mprob, "mprob");
    glp_add_cols(mprob, nr_vars[0]);
    for ( int i=1; i <= nr_vars[0]; ++i ) {
      glp_set_obj_coef(mprob, i, vals[i-1]);
    }
    glp_set_obj_dir(mprob, GLP_MIN); 

    for ( int j=1; j <= nr_vars[0]; ++j) {
      glp_set_col_bnds(mprob, j, GLP_DB, 0, 1);
    }
    /* set primary suppressions to 1! */
    for ( int k=0; k < len_prim[0]; ++k ) {
      glp_set_col_bnds(mprob, ind_prim[k], GLP_FX, 1, 1);
    }

    /* if such cells exist, fix them to zero */
    if ( len_fixed[0] > 0 ) {
      for ( int k=0; k < len_fixed[0]; ++k ) {
        glp_set_col_bnds(mprob, ind_fixed[k], GLP_FX, 0, 0);
      }
    }

    /* create new structure holding information on sdc problem */
    sdcinfo info;
    info.vals = v_ai; info.ia = ia; info.ja = ja; info.ar = ar;
    info.LPL = v_LPL; info.UPL = v_UPL; info.SPL = v_SPL;
    info.UB = v_UB; info.LB = v_LB; 
    info.ind_prim = ind_prim; info.len_prim = len_prim[0];
    info.ind_fixed= ind_fixed; info.len_fixed = len_fixed[0];
    info.cells_mat = cells_mat[0]; info.nr_vars = nr_vars[0]; info.nr_rows = nr_rows[0];
    info.upper_bound = -1; info.current_best_solution = xi_heur;
    info.tol = 1e-9;

    if ( verbose[0] == 0 ) {
      info.verbose = false;
    } else {
      info.verbose = true;
    }
    vector<int> bv; /* new */
    info.branchvars = bv; /* new */

    /* pointer to this object */
    sdcinfo * pinfo; pinfo = &info;      

    /* set up attackers problem with given solution of master problem */
    if ( info.verbose == true ) { 
      printf("setting up attacker's problem...");
    }
    glp_prob *aprob = setup_attacker_problem(pinfo, xi);
    if ( info.verbose == true ) {
      printf("[done]\n");
    }
    /* perform pre-processing and initialize constraint-pool */
    if ( info.verbose == true ) {
      printf("performing preprocessing phase...");
    }
    preprocess(aprob, mprob, pinfo, xi); 
    if ( info.verbose == true ) {
      printf("[done]\n");
    }
    
    /* set up incremental attacker problem */
    if ( info.verbose == true ) {
      printf("setup incremental attacker's problem...");
    }
    glp_prob *incprob = setup_incprob(pinfo, xi);
    if ( info.verbose == true ) {
      printf("[done]\n");
    }
    
    /* calculate a heuristic solution and return an upper_bound */
    use_existing_solution = 0;
    if ( info.verbose == true ) {
      printf("calculate a heuristic solution...");
    }
    heuristic_solution(incprob, pinfo, xi, use_existing_solution);
    //heuristic_solution_primitive(pinfo);
    if ( info.verbose == true ) {
      printf("[done] (upper_bound=%g)\n", info.upper_bound);
    }

    int upper_bound_original = info.upper_bound;

    /* solve relaxation of master problem */
    if ( info.verbose == true ) {
      printf("solve the relaxation of the master problem (%d constraints)...", glp_get_num_rows(mprob));
    }
    is_int = solve_relaxation(mprob, aprob, constraint_pool, pinfo, xi);    
    if ( info.verbose == true ) {
      printf("[done]\n");      
    }
    if ( is_int == true ) {
      if ( info.verbose == true ) {
        printf("we are finished and have a valid integer solution!\n");
      }
      if ( glp_get_obj_val(mprob) < info.upper_bound ) {
        info.upper_bound = glp_get_obj_val(mprob);
        /* updating current best solution */
        for ( int k=0; k<info.nr_vars; ++k ) {
          info.current_best_solution[k] = xi[k];
        }
      }
    } else {
      /* we have to do branch and bound to get an integer solution */
      if ( info.verbose == true ) {
        printf("we are not finished and have to enforce integrality on current solution using branch/bound!\n");
      }
      branch_and_bound(mprob, aprob, incprob, constraint_pool, pinfo, xi);    
    }
    if ( info.verbose == true ) {
      printf("algorithm has terminated!\n");
      printf("best solution has an objective val of %g.\n", info.upper_bound);
    }
    /* update final solution */
    for ( int i=0; i<info.nr_vars; ++i ) {
      final_pattern[i] = info.current_best_solution[i];
      xi[i] = final_pattern[i];
    }    
    bool res = is_valid_solution(aprob, mprob, constraint_pool, &info, xi);

    //if ( info.verbose == true ) {
    //  if ( res == true ) {
    //    printf("ok!\n");
    //  } else {
    //    printf("no valid solution!\n");
    //  }
    //}
  }
}
