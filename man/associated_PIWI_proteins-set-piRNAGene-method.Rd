% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/11_piRNAGene.R
\name{associated_PIWI_proteins<-,piRNAGene-method}
\alias{associated_PIWI_proteins<-,piRNAGene-method}
\title{Add/Remove Associated PIWI Proteins}
\usage{
\S4method{associated_PIWI_proteins}{piRNAGene}(object, action = "add") <- value
}
\arguments{
\item{object}{Object of class \code{piRNAGene}.}

\item{action}{Action to be done: can be 'add' or 'remove'.}

\item{value}{PIWI protein to add/remove.}
}
\value{
The modified object.
}
\description{
Setter for piRNA_associated_PIWI_proteins attribute of
\code{piRNAGene} class. By choosing an action, which can be
'add' or 'remove', it's possible to add and remove, respectively,
PIWI proteins from piRNA_associated_PIWI_proteins attribute.
}
\examples{
pir1 <- piRNAGene()
associated_PIWI_proteins(pir1)
associated_PIWI_proteins(pir1, action = 'add') <- 'PIWIL2'
associated_PIWI_proteins(pir1)
associated_PIWI_proteins(pir1, action = 'remove') <- 'PIWIL2'
associated_PIWI_proteins(pir1)
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
