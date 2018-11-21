##' @name manage_hierarchies
##' @rdname manage_hierarchies
##'
##' @title Create and modify the structure of hierarchies
##'
##' @description Functions \code{create_node()}, \code{add_nodes()}
##' and \code{delete_nodes()} allow to define and modify hierarchical structures
##' represented as trees. These objects can be used in \code{\link{makeProblem}} to
##' define the (hierarchical) structure of tables.
##'
##' @param total_lab the name of the overall total (summation over all contributing)
##' @param node a \code{} object as created in \code{create_node()} or returned from
##' \code{add_nodes()} or \code{delete_nodes()}.
##' @param node_labs character name(s) of new elements that should be inserted to or
##' deleted from a hierarchical structure
##' @param node_labs_new a character vector of new node names replacing the existing node names given in
##' argument \code{node_labs} when using \code{rename_nodes}.
##' @param reference_node character name of an existing node in the hierarchical
##' structure. When using \code{add_nodes()}, the new elements are created as children
##' of the reference node. In \code{delete_nodes()}, all children of the reference node that
##' match the names with argument \code{node_labs} are deleted from the hierarchy.
##' @examples
##' dim <- create_node(total_lab="Total")
##' dim <- add_nodes(dim, reference_node="Total", node_labs=LETTERS[1:4])
##' print(dim)
##'
##' ## add some levels below "A" and "C"
##' dim <- add_nodes(dim, reference_node="A", node_labs=paste0("a", 1:5))
##' dim <- add_nodes(dim, reference_node="C", node_labs=paste0("c", 1:5))
##' print(dim)
##'
##' ## delete some specific levels
##' dim <- delete_nodes(dim, node_labs=c("a1", "a4"))
##' print(dim)
##'
##' ## delete entire subtree
##' dim <- delete_nodes(dim, node_labs=c("C"))
##' print(dim)
##' # plot(dim)
##' @return a \code{Node} object that can be used to specify a hierarchy used
##' to define table inputs in \code{\link{protectTable}}.
NULL

##' @rdname manage_hierarchies
##' @export
create_node <- function(total_lab="Total", node_labs=NULL) {
  return(sdcHier_create(tot_lab=total_lab, node_labs=node_labs))
}

##' @rdname manage_hierarchies
##' @export
add_nodes <- function(node, node_labs, reference_node) {
  return(sdcHier_add(h=node, node_labs=node_labs, refnode=reference_node))
}

##' @rdname manage_hierarchies
##' @export
delete_nodes <- function(node, node_labs) {
  return(sdcHier_delete(h=node, node_labs=node_labs))
}

##' @rdname manage_hierarchies
##' @export
rename_nodes <- function(node, node_labs, node_labs_new) {
  return(sdcHier_rename(h=node, node_labs=node_labs, node_labs_new=node_labs_new))
}

# internal function used in hierarchies before sdcTable
# --> structure required for sdcTable (already works)
node_to_sdcinput <- function(inp, addNumLevels=FALSE) {
  dt <- data.table(sdcHier_convert(h=inp, format="data.frame"))
  setnames(dt, c("levels","codes"))
  if (addNumLevels==TRUE) {
    num_levels <- NULL
    dt[,num_levels:=nchar(dt[,levels])]
  }
  return(dt[])
}
