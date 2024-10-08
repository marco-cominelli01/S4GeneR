% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07_RegulatoryNonCodingGene.R
\name{targets_id}
\alias{targets_id}
\title{Get Targets}
\usage{
targets_id(object)
}
\arguments{
\item{object}{Object of class \code{RegulatoryNonCodingGene}
(or derived classes).}
}
\value{
List of targets ENSEMBL IDs.
}
\description{
Getter for targets_ensembl_id attribute of
\code{RegulatoryNonCodingGene} class (and derived classes).
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
mir1 <- miRNAGene(targets_ensembl_id = list("ENST12341234090"))
targets_id(mir1)
}
\seealso{
\link{targets_id,RegulatoryNonCodingGene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
