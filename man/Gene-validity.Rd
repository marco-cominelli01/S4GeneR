% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_Gene.R
\name{Gene-validity}
\alias{Gene-validity}
\title{Validity Check for Gene Virtual Class}
\arguments{
\item{object}{Object of class \code{Gene} (or derived classes).}
}
\value{
\code{TRUE} if the object is valid, otherwise the function
throws an error via \code{stop()}.
}
\description{
The function checks that: each slot, except the special ones,
is of length 1 (or empty); each ID is a valid ENSEMBL ID; all coordinates
are consistent w.r.t. each other and that gene sequence is consistent
with the provided coordinates.
}
\seealso{
\linkS4class{Gene}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
