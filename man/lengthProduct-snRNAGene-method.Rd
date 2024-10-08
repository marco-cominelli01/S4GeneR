% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/15_lengthProduct.R
\name{lengthProduct,snRNAGene-method}
\alias{lengthProduct,snRNAGene-method}
\title{Get Mature snRNA Length}
\usage{
\S4method{lengthProduct}{snRNAGene}(object)
}
\arguments{
\item{object}{Object of class \code{snRNAGene}.}
}
\value{
Mature snRNA length.
}
\description{
Compute length of snRNA_mature_sequence attribute of
\code{snRNAGene} class.
}
\examples{
sn1 <- snRNAGene(snRNA_mature_sequence = 'UAUUAUAUAA')
lengthProduct(sn1)
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
