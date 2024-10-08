% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/06_rRNAGene.R
\name{rRNA_type<-}
\alias{rRNA_type<-}
\title{Set rRNA Type}
\usage{
rRNA_type(object) <- value
}
\arguments{
\item{object}{Object of class \code{rRNAGene}.}

\item{value}{New value for rRNA_type}
}
\value{
The modified object.
}
\description{
Setter for rRNA_type attribute of \code{rRNAGene} class.
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
r1 <- rRNAGene(rRNA_type = '5S')
rRNA_type(r1)
rRNA_type(r1) <- '18S'
rRNA_type(r1)
}
\seealso{
\link{rRNA_type<-,rRNAGene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
