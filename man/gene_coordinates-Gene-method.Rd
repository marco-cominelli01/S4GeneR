% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_Gene.R
\name{gene_coordinates,Gene-method}
\alias{gene_coordinates,Gene-method}
\title{Get Gene Genomic Coordinates}
\usage{
\S4method{gene_coordinates}{Gene}(object)
}
\arguments{
\item{object}{Object of class \code{Gene} (and derived classes).}
}
\value{
\code{GRanges} object representing gene coordinates.
}
\description{
Getter for gene_coordinates attribute of \code{Gene} class
(and derived classes).
}
\examples{
# Since Gene is a virtual class, the example is made using a derived
# class.
cg1 <- CodingGene(chromosome = 'chr1', strand = '+', gene_start = 210,
gene_end = 800)
gene_coordinates(cg1)
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
