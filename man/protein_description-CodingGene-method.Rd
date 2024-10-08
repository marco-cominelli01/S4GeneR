% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_CodingGene.R
\name{protein_description,CodingGene-method}
\alias{protein_description,CodingGene-method}
\title{Get Protein Description}
\usage{
\S4method{protein_description}{CodingGene}(object)
}
\arguments{
\item{object}{Object of class \code{CodingGene}.}
}
\value{
Protein description.
}
\description{
Getter for protein_description attribute of \code{CodingGene}
class.
}
\examples{
cg1 <- CodingGene(protein_description = "Transcription factor")
protein_description(cg1)
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
