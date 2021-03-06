% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generics_simpleTriplet.r,
%   R/methods_class_simpleTriplet.r
\docType{methods}
\name{get.simpleTriplet}
\alias{get.simpleTriplet}
\alias{get.simpleTriplet,simpleTriplet,character,list-method}
\title{query \code{simpleTriplet}-objects depending on argument \code{type}}
\usage{
get.simpleTriplet(object, type, input)

\S4method{get.simpleTriplet}{simpleTriplet,character,list}(object, type,
  input)
}
\arguments{
\item{object}{an object of class \code{simpleTriplet}}

\item{type}{a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
\itemize{
\item rowInd: extract all row-indices
\item colInd: extract all column-indices
\item values: extract all values
\item nrRows: return the number of rows of the input object
\item nrCols: return the number of columns of the input object
\item nrCells: return the number of cells (different from 0!)
\item duplicatedRows: return a numeric vector showing indices of duplicated rows
\item transpose: transpose input \code{object} and return the transposed matrix
\item getRow: return a specific row of input \code{object}
\item getCol: return a specific column of input \code{object}}

\item{input}{a list depending on argument \code{type}.}
\itemize{
\item type == 'getRow': input is a list of length 1
\itemize{
\item first element: numeric vector of length 1 defining index of row that is to be returned }
\item type == 'getCol': input is a list of length 1
\itemize{
\item first element: numeric vector of length 1 defining index of column that is to be returned }
\item else: input is not used at all (empty list)}
}
\value{
information from \code{object} depending on \code{type}
\itemize{
\item a numeric vector if type matches 'rowInd', 'colInd', 'values', 'nrRows', 'nrCols', 'nrCells' or 'duplicatedRows'
\item an object of class \code{simpleTriplet} if type matches 'transpose', 'getRow' or 'getCol'
}
}
\description{
query \code{simpleTriplet}-objects depending on argument \code{type}
}
\note{
internal function
}
\author{
Bernhard Meindl \email{bernhard.meindl@statistik.gv.at}
}
\keyword{internal}
