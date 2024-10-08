% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02_CodingGene.R
\name{protein_description<-}
\alias{protein_description<-}
\title{Set Protein Description}
\usage{
protein_description(object) <- value
}
\arguments{
\item{object}{Object of class \code{CodingGene}.}

\item{value}{New value for protein_description.}
}
\value{
The modified object.
}
\description{
Setter for protein_description attribute of \code{CodingGene}
class.
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
cg1 <- CodingGene(protein_description = 'Transcription factor')
protein_description(cg1)
protein_description(cg1) <- 'Chromatin remodeling factor'
protein_description(cg1)
}
\seealso{
\link{protein_description<-,CodingGene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
