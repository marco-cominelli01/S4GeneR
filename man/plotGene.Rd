% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/16_PlotGene.R
\name{plotGene}
\alias{plotGene}
\title{Plot Gene Exons Structure}
\usage{
plotGene(x, transcr_label, transcr_id)
}
\arguments{
\item{x}{An object of class \code{Gene} (and derived classes).}

\item{transcr_label}{Name of the transcript to plot. It's used to name the
respective y-tick (e.g. "mRNA").}

\item{transcr_id}{Slot containing the transcript ENSEMBL ID.}
}
\value{
The gene plot (without x-ticks: this function is not meant to be
used alone).
}
\description{
The function, relying on ggplot2, plots the gene's exons and
prepares the plot for further additions (mRNA or pre-nc RNA).
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
