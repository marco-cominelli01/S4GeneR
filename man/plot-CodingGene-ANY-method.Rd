% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/16_PlotGene.R
\name{plot,CodingGene,ANY-method}
\alias{plot,CodingGene,ANY-method}
\title{Plot CodingGene Objects}
\usage{
\S4method{plot}{CodingGene,ANY}(x, y = 0)
}
\arguments{
\item{x}{Object of class \code{CodingGene}.}

\item{y}{Required parameter to extend the generic function 'plot()'.
Its value doesn't influence the final plot.}
}
\value{
Gene plot.
}
\description{
The function plots the gene's exons and, if provided,
the mRNA and CDS. To plot the CDS, mRNA coordinates are necessary,
otherwise only the gene's exons are plotted.
}
\examples{
cg1 <- CodingGene(chromosome = 'chrM', strand = '-',
                exons_starts = c(30,60,90),
                exons_ends = c(40,70,100),
                mRNA_start = 35, mRNA_end = 67,
                cds_start = 37, cds_end = 62,
                gene_ensembl_id = "ENSG12312300000",
                mRNA_ensembl_id = "ENST00011133300")
plot(cg1)
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
