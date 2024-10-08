% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/15_lengthProduct.R
\name{lengthProduct}
\alias{lengthProduct}
\title{Get Gene Product Length}
\usage{
lengthProduct(object)
}
\arguments{
\item{object}{Object of a class derived from \code{Gene}.
The specific type of object depends on the method implementation.}
}
\value{
Gene product length.
}
\description{
Compute length of the gene product.
The specific attribute depends on the method implementation.
}
\examples{
# This is only the generic function, look in 'See Also' section
# for the method implementation.
# Example on a CodingGene object
cg1 <- CodingGene(protein_sequence = 'KRTA')
lengthProduct(cg1)
}
\seealso{
\link{lengthProduct,CodingGene-method},
\link{lengthProduct,tRNAGene-method},
\link{lengthProduct,tRNAGene-method},
\link{lengthProduct,miRNAGene-method},
\link{lengthProduct,siRNAGene-method},
\link{lengthProduct,piRNAGene-method},
\link{lengthProduct,snRNAGene-method},
\link{lengthProduct,snoRNAGene-method},
\link{lengthProduct,LongNonCodingGene-method},
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
