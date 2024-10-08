% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_Gene.R
\name{gene_sequence<-}
\alias{gene_sequence<-}
\title{Set Gene Sequence}
\usage{
gene_sequence(object) <- value
}
\arguments{
\item{object}{Object of class \code{Gene} (and derived classes).}

\item{value}{New value for gene_sequence.}
}
\value{
The modified object.
}
\description{
Setter for gene_sequence attribute of \code{Gene} class
(and derived classes).
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
cg1 <- CodingGene(gene_sequence = 'ATCCAT')
gene_sequence(cg1)
gene_sequence(cg1) <- 'AAAAAAAA'
gene_sequence(cg1)
}
\seealso{
\link{gene_sequence<-,Gene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
