% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_rRNAGene.R
\name{mature_sequence<-,rRNAGene-method}
\alias{mature_sequence<-,rRNAGene-method}
\title{Set Mature rRNA Sequence}
\usage{
\S4method{mature_sequence}{rRNAGene}(object) <- value
}
\arguments{
\item{object}{Object of class \code{rRNAGene}.}

\item{value}{New value for rRNA_mature_sequence.}
}
\value{
The modified object.
}
\description{
Setter for rRNA_mature_sequence attribute of \code{rRNAGene} class.
}
\examples{
r1 <- rRNAGene(rRNA_mature_sequence = 'AUUA')
mature_sequence(r1)
mature_sequence(r1) <- 'AAAGCA'
mature_sequence(r1)
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
