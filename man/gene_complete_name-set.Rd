% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_Gene.R
\name{gene_complete_name<-}
\alias{gene_complete_name<-}
\title{Set Gene Complete Name}
\usage{
gene_complete_name(object) <- value
}
\arguments{
\item{object}{Object of class \code{Gene} (and derived classes).}

\item{value}{New value for gene_complete_name.}
}
\value{
The modified object.
}
\description{
Setter for gene_complete_name attribute of \code{Gene} class
(and derived classes).
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
cg1 <- CodingGene(gene_complete_name = 'Tumor Protein 53')
gene_complete_name(cg1)
gene_complete_name(cg1) <- 'Breast cancer type 1 susceptibility protein'
gene_complete_name(cg1)
}
\seealso{
\link{gene_complete_name<-,Gene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
