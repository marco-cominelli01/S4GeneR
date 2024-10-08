% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/07_RegulatoryNonCodingGene.R
\name{regulatory_functions,RegulatoryNonCodingGene-method}
\alias{regulatory_functions,RegulatoryNonCodingGene-method}
\title{Get Regulatory Functions}
\usage{
\S4method{regulatory_functions}{RegulatoryNonCodingGene}(object)
}
\arguments{
\item{object}{Object of class \code{RegulatoryNonCodingGene}
(and derived classes).}
}
\value{
List of regulatory functions.
}
\description{
Getter for regulatory_functions attribute of
\code{RegulatoryNonCodingGene} class (and derived classes).
}
\examples{
# Since RegulatoryNonCodingGene is a virtual class, the example is made
# using a derived class
mir1 <- miRNAGene(regulatory_functions = list("Transcription modulation"))
regulatory_functions(mir1)
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
