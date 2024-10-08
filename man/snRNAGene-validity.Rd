% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/12_snRNAGene.R
\name{snRNAGene-validity}
\alias{snRNAGene-validity}
\title{Validity Check for snRNAGene Class}
\arguments{
\item{object}{Object of class \code{snRNAGene}.}
}
\value{
\code{TRUE} if the object is valid, otherwise the function
throws an error via \code{stop()}.
}
\description{
The function checks that: each slots has a length of 1 (or
empty), the inserted spliceosome complex is valid and the mature sequence
is consistent with the gene sequence and the defined coordinates.
}
\seealso{
\link{RegulatoryNonCodingGene-validity}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
