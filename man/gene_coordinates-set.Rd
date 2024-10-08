% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_Gene.R
\name{gene_coordinates<-}
\alias{gene_coordinates<-}
\title{Set Gene Genomic Coordinates}
\usage{
gene_coordinates(object) <- value
}
\arguments{
\item{object}{Object of class \code{Gene} (and derived classes).}

\item{value}{Named list with chromosome, strand, start and end values
to build the new gene coordinates.}
}
\value{
The modified object.
}
\description{
Setter for gene_coordinates attribute of \code{Gene} class
(and derived classes).
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
cg1 <- CodingGene(chromosome = 'chr1', strand = '+', gene_start = 20,
gene_end = 300)
gene_coordinates(cg1)
gene_coordinates(cg1) <- list(chromosome = 'chr2', strand = '+',
start = 30, end = 50)
gene_coordinates(cg1)
}
\seealso{
\link{gene_coordinates<-,Gene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
