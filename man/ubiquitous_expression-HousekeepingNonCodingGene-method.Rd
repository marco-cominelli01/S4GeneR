% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04_HousekeepingNonCodingGene.R
\name{ubiquitous_expression,HousekeepingNonCodingGene-method}
\alias{ubiquitous_expression,HousekeepingNonCodingGene-method}
\title{Get Ubiquitous Expression}
\usage{
\S4method{ubiquitous_expression}{HousekeepingNonCodingGene}(object)
}
\arguments{
\item{object}{Object of class \code{HousekeepingNonCodingGene}
(and derived classes).}
}
\value{
Value of ubiquitous_expression.
}
\description{
Getter for ubiquitous_expression attribute of
\code{HousekeepingNonCodingGene} class (and derived classes).
}
\examples{
# Since HousekeepingNonCodingGene is a virtual class, the example is
# made using a derived class
t1 <- tRNAGene(ubiquitous_expression = FALSE)
ubiquitous_expression(t1)
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
