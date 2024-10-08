% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05_tRNAGene.R
\name{amino_acid}
\alias{amino_acid}
\title{Get tRNA Aminoacid}
\usage{
amino_acid(object)
}
\arguments{
\item{object}{Object of class \code{tRNAGene}.}
}
\value{
Aminoacid sequence (one letter, or empty).
}
\description{
Getter for amino_acid attribute of \code{tRNAGene} class.
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
t1 <- tRNAGene(amino_acid = 'K')
amino_acid(t1)
}
\seealso{
\link{amino_acid,tRNAGene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
