% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/05_tRNAGene.R
\name{anticodon<-}
\alias{anticodon<-}
\title{Set tRNA Anticodon}
\usage{
anticodon(object) <- value
}
\arguments{
\item{object}{Object of class \code{tRNAGene}.}

\item{value}{New value for anticodon.}
}
\value{
The modified object.
}
\description{
Setter for anticodon attribute of \code{tRNAGene} class.
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
t1 <- tRNAGene(anticodon = 'GGG')
anticodon(t1)
anticodon(t1) <- 'UUU'
anticodon(t1)
}
\seealso{
\link{anticodon<-,tRNAGene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
