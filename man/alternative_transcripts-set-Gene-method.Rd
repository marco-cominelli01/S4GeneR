% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/01_Gene.R
\name{alternative_transcripts<-,Gene-method}
\alias{alternative_transcripts<-,Gene-method}
\title{Add/Remove Alternative Transcripts}
\usage{
\S4method{alternative_transcripts}{Gene}(object, action = "add") <- value
}
\arguments{
\item{object}{Object of class \code{Gene} (and derived classes).}

\item{action}{Action to be done: can be 'add' or 'remove'.}

\item{value}{Named list with transcript_id of the transcript to add/remove
and protein_coding, a \code{logical} to indicate the nature of the transcript.}
}
\value{
The modified object.
}
\description{
Setter for alternative_transcripts attribute of \code{Gene} class
(and derived classes). By choosing an action, which can be 'add' or
'remove', it's possible to add and remove, respectively, transcripts from
alternative_transcripts attribute.
}
\examples{
# Since Gene is a virtual class, the example is made using a derived
# class.
cg1 <- CodingGene()
alternative_transcripts(cg1)
alternative_transcripts(cg1, action = 'add') <- list(transcript_id =
'ENST12312312312', protein_coding = TRUE)
alternative_transcripts(cg1)
alternative_transcripts(cg1, action = 'remove') <- list(transcript_id =
'ENST12312312312', protein_coding = TRUE)
alternative_transcripts(cg1)
}
\author{
Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
}
