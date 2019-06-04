# sdcTable 0.28
- new function `contributing_indices()`
- call `.Defunct` for old exported functions
- make default bounds (`ub`, `lb`) depend on costs
- allow to specify numVarInd also by name via `numVarName` in primarySupression() for dominance rules (`p`, `nk` and `pq`)
- make use of sampling weights (by replicating) values for dominance rules.
- don't export unused functionality
- `createJJFormat()` and `writeJJFormat()` can be used to create text files with
sdcProblems in "JJs Format"
- allow to specify variable names too in `makeProblem`

# sdcTable 0.27
- use functionality from [`sdcHierarchies`](https://cran.r-project.org/package=sdcHierarchies) to build hierarchies.     - removed `create_node()`; please use `hier_create()` instead
    - removed `add_nodes()`; please use `hier_add()` instead
    - removed `delete_nodes()`; please use `hier_delete()` instead
    - removed `rename_nodes()`; please use `hier_rename()` instead

# sdcTable 0.26
* bugfix when computing indices for contributing units used in some primary suppression methods

# sdcTable 0.25
* better error messages in case invalid hierarchies have been specified

# sdcTable 0.24
* update in SIMPLEHEURISTIC-algorithm that now really respects cells (status `"z"`) that must not be suppressed

# sdcTable 0.23
* bugfix in SIMPLEHEURISTIC-algorithm
* code cleanup
* updated package vignette
* new functions `create_node()`, `add_nodes()` and `delete_nodes()` that allow a dynamic generation of hierarchies. For an example have a look at `?makeProblem`

# sdcTable 0.22.x
* Feature to create BATCH files suitable for [**tau-argus**](https://github.com/sdcTools/tauargus)
