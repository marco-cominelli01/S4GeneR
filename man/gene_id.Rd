% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_Gene.R
\name{gene_id}
\alias{gene_id}
\title{Get Gene ENSEMBL ID}
\usage{
gene_id(object)
}
\arguments{
\item{object}{Object of class \code{Gene} (and derived classes).}
}
\value{
Value of gene_ensembl_id.
}
\description{
Getter for gene_ensembl_id attribute of \code{Gene} class
(and derived classes).
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
cg1 <- CodingGene(gene_ensembl_id = 'ENSG12312312300')
gene_id(cg1)
}
\seealso{
\link{gene_id,Gene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
