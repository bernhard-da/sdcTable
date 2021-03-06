% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tauBatch_helpers.R
\name{argusVersion}
\alias{argusVersion}
\title{argusVersion}
\usage{
argusVersion(exe, verbose = FALSE)
}
\arguments{
\item{exe}{a path to a tau-argus executable}

\item{verbose}{(logical) if \code{TRUE}, the version info and build number
of the given tau-argus executable will be printed.}
}
\value{
a list with two elements being the tau-argus version and the build-number.
}
\description{
returns the version and build number of a given tau-argus executable
specified in argument \code{exe}.
}
\examples{
\dontrun{
argusVersion(exe="C:\\\\Tau\\\\TauArgus.exe", verbose=TRUE)
}
}
