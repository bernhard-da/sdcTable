#' query \code{sdcProblem}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{sdcProblem}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item dataObj: a list containing the (raw) input data
#' \item problemInstance: return the current problem instance
#' \item partition: a list containing information on the subtables that are required to be protected as well as information on the processing order of the subtables
#' \item elapsedTime: the elapsed time of the protection algorithm so far
#' \item dimInfo: information on the variables defining the hierarchical table
#' \item indicesDealtWith: a set of indices that have already been dealt with during the protection algorithmus
#' \item startI: current level at which subtables need to be protected (useful when restarting HITAS|HYPERCUBE)
#' \item startJ: current number of the subtable within a given level that needs to be protected (useful when restarting HITAS|HYPERCUBE)
#' \item innerAndMarginalCellInfo: for a given problem, get indices of inner- and marginal table cells
#'
#' @return information from objects of class \code{sdcProblem} depending on argument \code{type}
#' \itemize{
#' \item an object of class \code{dataObj} (or NULL) if \code{type} matches 'dataObj'
#' \item an object of class \code{problemInstance} (or NULL) if \code{type} matches 'problemInstance'
#' \item a list (or NULL) if argument \code{type} matches 'partition' containing the following elements:
#' \itemize{
#' \item element 'groups': list with each list-element being a character vector specifying a specific level-group
#' \item element 'indices': list with each list-element being a numeric vector defining indices of a subtable
#' \item element 'strIDs': list with each list-element being a character vector defining IDs of a subtable
#' \item element 'nrGroups': numeric vector of length 1 defining the total number of groups that have to be considered
#' \item element 'nrTables': numeric vector of length 1 defining the total number of subtables that have to be considered}
#' \item a list (or NULL) if argument \code{type} matches 'innerAndMarginalCellInfo' containing the following elements:
#' \itemize{
#' \item element 'innerCells': character vector specifying ID's of inner cells
#' \item element 'totCells': character vector specifying ID's of marginal cells
#' \item element 'indexInnerCells': numeric vector specifying indices of inner cells
#' \item element 'indexTotCells': numeric vector specifying indices of marginal cells}
#' \item an object of class \code{dimInfo} (or NULL) if \code{type} matches 'dimInfo'
#' \item numeric vector if argument \code{type} matches 'elapsedTime'
#' \item numeric vector of length 1 if argument \code{type} matches 'startI' or 'startJ'
#' }
#'
#' @export
#' @docType methods
#' @rdname get.sdcProblem-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('get.sdcProblem', function(object, type) {standardGeneric('get.sdcProblem')})

#' modify \code{sdcProblem}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{sdcProblem}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:
#' \itemize{
#' \item problemInstance: set|modify slot 'problemInstance' of argument \code{object}
#' \item partition: set|modify slot 'partition' of argument \code{object}
#' \item startI: set|modify slot 'startI' of argument \code{object}
#' \item startJ: set|modify slot 'startJ' of argument \code{object}
#' \item indicesDealtWith: set|modify slot 'indicesDealtWith' of argument \code{object}
#' \item elapsedTime: set|modify slot 'elapsedTime' of argument \code{object}}
#' @param input a list with elements depending on argument \code{type}.
#' \itemize{
#' \item an object of class \code{problemInstance} if argument \code{type} matches 'problemInstance'
#' \item a numeric vector of length 1 if argument \code{type} matches 'startI', 'startJ' or 'elapsedTime'
#' \item a numeric vector if argument \code{type} matches 'indicesDealtWith'}
#'
#' @return an object of class \code{sdcProblem}
#'
#' @export
#' @docType methods
#' @rdname set.sdcProblem-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('set.sdcProblem', function(object, type, input) {standardGeneric('set.sdcProblem')})

#' perform calculations on \code{sdcProblem}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{sdcProblem}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item rule.freq: modify suppression status within \code{object} according to frequency suppression rule
#' \item rule.nk: modify sdcStatus of \code{object} according to nk-dominance rule
#' \item rule.p: modify sdcStatus of \code{object} according to p-percent rule
#' \item rule.pq: modify sdcStatus of \code{object} according to pq-rule
#' \item cellID: find index of cell defined by information provided with argument \code{input}
#' \item finalize: create an object of class \code{safeObj}
#' \item contributingIndices: calculate indices within the current problem that contribute to a given cell
#' @param input a list depending on argument \code{type}.}
#' \itemize{
#' \item a list (typically generated using genParaObj()) specifying parameters for primary cell suppression if argument \code{type} matches 'rule.freq', 'rule.nk' or 'rule.p'
#' \item a list of length 3 if argument \code{type} matches 'cellID' having following elements
#' \itemize{
#' \item first element: character vector specifying variable names that need to exist in slot 'dimInfo' of \code{object}
#' \item second element: character vector specifying codes for each variable that define a specific table cell
#' \item third element: logical vector of length 1 with TRUE setting verbosity and FALSE to turn verbose output off}
#' \item a list of length 1 if argument \code{type} matches 'contributingIndices' having following element
#' \itemize{
#' \item first element: character vector of length 1 being an ID for which contributing indices should be calculated }
#' @return information from objects of class \code{sdcProblem} depending on argument \code{type}
#' \itemize{
#' \item an object of class \code{sdcProblem} if argument \code{type} matches 'rule.freq', 'rule.nk' or 'rule.p'
#' \item a numeric vector of length 1 specifying the index of the cell of interest if argument \code{type} matches 'cellID'
#' \item an object of class \code{safeObj} if argument \code{type} matches 'finalize'
#' \item a numeric vector with indices that contribute to the desired table cell if argument \code{type} matches 'contributingIndices'
#' }
#' @export
#' @docType methods
#' @rdname calc.sdcProblem-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('calc.sdcProblem', function(object, type, input) {standardGeneric('calc.sdcProblem')})

# get-methods
setGeneric("g_problemInstance", function(object) { standardGeneric("g_problemInstance") })
setGeneric("g_dimInfo", function(object) { standardGeneric("g_dimInfo") })
setGeneric("g_partition", function(object) { standardGeneric("g_partition") })
setGeneric("g_elapsedTime", function(object) { standardGeneric("g_elapsedTime") })
setGeneric("g_dataObj", function(object) { standardGeneric("g_dataObj") })
setGeneric("g_startI", function(object) { standardGeneric("g_startI") })
setGeneric("g_startJ", function(object) { standardGeneric("g_startJ") })
setGeneric("g_indicesDealtWith", function(object) { standardGeneric("g_indicesDealtWith") })
setGeneric("g_df", function(object, ...) { standardGeneric("g_df") })

# set-methods
setGeneric("s_problemInstance<-", function(object, value) standardGeneric("s_problemInstance<-"))
setGeneric("s_partition<-", function(object, value) standardGeneric("s_partition<-"))
setGeneric("s_startI<-", function(object, value) standardGeneric("s_startI<-"))
setGeneric("s_startJ<-", function(object, value) standardGeneric("s_startJ<-"))
setGeneric("s_indicesDealtWith<-", function(object, value) standardGeneric("s_indicesDealtWith<-"))
setGeneric("s_elapsedTime<-", function(object, value) standardGeneric("s_elapsedTime<-"))

# calc-methods
setGeneric("c_rule_freq", function(object, input) { standardGeneric("c_rule_freq") })
setGeneric("c_rule_nk", function(object, input) { standardGeneric("c_rule_nk") })
setGeneric("c_rule_nk", function(object, input) { standardGeneric("c_rule_nk") })
setGeneric("c_rule_p", function(object, input) { standardGeneric("c_rule_p") })
setGeneric("c_rule_pq", function(object, input) { standardGeneric("c_rule_pq") })
setGeneric("c_quick_suppression", function(object, input) { standardGeneric("c_quick_suppression") })
setGeneric("c_cellID", function(object, input) { standardGeneric("c_cellID") })
setGeneric("c_finalize", function(object, input) { standardGeneric("c_finalize") })
setGeneric("c_contributing_indices", function(object, input) { standardGeneric("c_contributing_indices") })
