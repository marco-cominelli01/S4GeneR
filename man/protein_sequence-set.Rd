% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_CodingGene.R
\name{protein_sequence<-}
\alias{protein_sequence<-}
\title{Set Protein Sequence}
\usage{
protein_sequence(object) <- value
}
\arguments{
\item{object}{Object of class \code{CodingGene}.}

\item{value}{New value for protein_sequence.}
}
\value{
The modified object.
}
\description{
Setter for protein_sequence attribute of \code{CodingGene} class.
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
cg1 <- CodingGene(protein_sequence = 'KRT')
protein_sequence(cg1)
protein_sequence(cg1) <- 'KKK'
protein_sequence(cg1)
}
\seealso{
\link{protein_sequence<-,CodingGene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
