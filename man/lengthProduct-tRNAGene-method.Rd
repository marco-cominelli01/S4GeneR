% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/15_lengthProduct.R
\name{lengthProduct,tRNAGene-method}
\alias{lengthProduct,tRNAGene-method}
\title{Get Mature tRNA Length}
\usage{
\S4method{lengthProduct}{tRNAGene}(object)
}
\arguments{
\item{object}{Object of class \code{tRNAGene}.}
}
\value{
Mature tRNA length.
}
\description{
Compute length of tRNA_mature_sequence attribute of \code{tRNAGene}
class.
}
\examples{
t1 <- tRNAGene(tRNA_mature_sequence = 'UAUUAUAUAA')
lengthProduct(t1)
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
