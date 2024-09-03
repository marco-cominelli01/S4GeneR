#' @title SmallNonCodingGene Virtual Class
#' @description This virtual class is the class from which
#' miRNAGene, siRNAGene, piRNAGene, snRNAGene and snoRNAGene classes are
#' derived.
#' The class has no additional slots w.r.t. RegulatoryNonCodingGene class
#' because it was created just for the sake of clarity.
#' @return No objects can be created directly from this class since
#' it's virtual.
#' @seealso
#' \linkS4class{RegulatoryNonCodingGene}
#' @examples
#' # The SmallNonCodingGene class is virtual, so no objects can
#' # be created directly from this class. It needs to be extended by
#' # other classes.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("SmallNonCodingGene",
        contains = c("RegulatoryNonCodingGene","VIRTUAL"))
