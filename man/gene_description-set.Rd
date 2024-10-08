% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_Gene.R
\name{gene_description<-}
\alias{gene_description<-}
\title{Set Gene Description}
\usage{
gene_description(object) <- value
}
\arguments{
\item{object}{Object of class \code{Gene} (and derived classes).}

\item{value}{New value for gene_description.}
}
\value{
The modified object.
}
\description{
Setter for gene_description attribute of \code{Gene} class
(and derived classes).
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
cg1 <- CodingGene(gene_description = 'TP53 is a tumor suppressor gene')
gene_description(cg1)
gene_description(cg1) <- 'BRCA1 is involved in breast cancer onset'
gene_description(cg1)
}
\seealso{
\link{gene_description<-,Gene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
