% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/03_NonCodingGene.R
\name{pre_ncRNA_id}
\alias{pre_ncRNA_id}
\title{Get pre-ncRNA ENSEMBL ID}
\usage{
pre_ncRNA_id(object)
}
\arguments{
\item{object}{Object of class \code{NonCodingGene} (and derived classes).}
}
\value{
Value of pre_ncRNA_ensembl_id.
}
\description{
Getter for pre_ncRNA_ensembl_id attribute of \code{NonCodingGene}
class (and derived classes).
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
mir1 <- miRNAGene(pre_ncRNA_ensembl_id = 'ENST12312312300')
pre_ncRNA_id(mir1)
}
\seealso{
\link{pre_ncRNA_id,NonCodingGene-method}
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
