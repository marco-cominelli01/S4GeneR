#' @title HousekeepingNonCodingGene Virtual Class
#' @description This virtual class is the class from which
#' tRNAGene and rRNAGene classes are derived.
#' @slot essentiality_score `numeric`. Gene essentiality score.
#' @slot ubiquitous_expression `logical`.
#' Gene ubiquitously expressed or not.
#' @return No objects can be created directly from this class since
#' it's virtual.
#' @seealso
#' \linkS4class{NonCodingGene}, \linkS4class{tRNAGene},
#' \linkS4class{rRNAGene}
#' @examples
#' # The HousekeepingNonCodingGene class is virtual, so no objects can
#' # be created directly from this class. It needs to be extended by
#' # other classes.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("HousekeepingNonCodingGene",
        contains = c("NonCodingGene","VIRTUAL"),
        slots = c(
            essentiality_score = 'numeric',
            ubiquitous_expression = 'logical'
        ),
        prototype = list(
            essentiality_score = NA_real_,
            ubiquitous_expression = NA
        ) )


#' @title Validity Check for HousekeepingNonCodingGene Virtual Class
#' @description The function checks that: each slot, except the special ones,
#' is of length 1 (or empty).
#' @param object Object of class `HousekeepingNonCodingGene`
#' (and derived classes).
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \link{NonCodingGene-validity}
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name HousekeepingNonCodingGene-validity

setValidity("HousekeepingNonCodingGene", function(object) {

    ### 1. Check on slots length (they must contain 1 element, except
    ###    the special ones).
    check_slots_length(obj = object, class_obj = "HousekeepingNonCodingGene",
                        mother_class_obj = "NonCodingGene")

    return(TRUE) })


# GETTERs for HousekeepingNonCodingGene

#' @title Get Gene Essentiality Score
#' @description Getter for essentiality_score attribute of
#' `HousekeepingNonCodingGene` class (and derived classes).
#' @param object Object of class `HousekeepingNonCodingGene`
#' (and derived classes).
#' @return Value of essentiality_score.
#' @seealso
#' \link{essentiality_score,HousekeepingNonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' t1 <- tRNAGene(essentiality_score = 12)
#' essentiality_score(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("essentiality_score",
            function(object) standardGeneric("essentiality_score"))


#' @title Get Gene Essentiality Score
#' @description Getter for essentiality_score attribute of
#' `HousekeepingNonCodingGene` class (and derived classes).
#' @param object Object of class `HousekeepingNonCodingGene`
#' (and derived classes).
#' @return Value of essentiality_score.
#' @examples
#' # Since HousekeepingNonCodingGene is a virtual class, the example is made
#' # using a derived class
#' t1 <- tRNAGene(essentiality_score = 12)
#' essentiality_score(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("essentiality_score", "HousekeepingNonCodingGene", function(object){
    object@essentiality_score })


#' @title Get Ubiquitous Expression
#' @description Getter for ubiquitous_expression attribute of
#' `HousekeepingNonCodingGene` class (and derived classes).
#' @param object Object of class `HousekeepingNonCodingGene`
#' (and derived classes).
#' @return Value of ubiquitous_expression.
#' @seealso
#' \link{ubiquitous_expression,HousekeepingNonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' t1 <- tRNAGene(ubiquitous_expression = FALSE)
#' ubiquitous_expression(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("ubiquitous_expression",
            function(object) standardGeneric("ubiquitous_expression"))

#' @title Get Ubiquitous Expression
#' @description Getter for ubiquitous_expression attribute of
#' `HousekeepingNonCodingGene` class (and derived classes).
#' @param object Object of class `HousekeepingNonCodingGene`
#' (and derived classes).
#' @return Value of ubiquitous_expression.
#' @examples
#' # Since HousekeepingNonCodingGene is a virtual class, the example is
#' # made using a derived class
#' t1 <- tRNAGene(ubiquitous_expression = FALSE)
#' ubiquitous_expression(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("ubiquitous_expression",
            "HousekeepingNonCodingGene", function(object) {
                object@ubiquitous_expression })


# SETTERs for HousekeepingNonCodingGene

#' @title Set Gene Essentiality Score
#' @description Setter for essentiality_score attribute of
#' `HousekeepingNonCodingGene` class (and derived classes).
#' @param object Object of class `HousekeepingNonCodingGene`
#' (and derived classes).
#' @param value New value for essentiality_score.
#' @return The modified object.
#' @seealso
#' \link{essentiality_score<-,HousekeepingNonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' t1 <- tRNAGene(essentiality_score = 12)
#' essentiality_score(t1)
#' essentiality_score(t1) <- 20
#' essentiality_score(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("essentiality_score<-",
            function(object, value) standardGeneric("essentiality_score<-"))


#' @title Set Gene Essentiality Score
#' @description Setter for essentiality_score attribute of
#' `HousekeepingNonCodingGene` class (and derived classes).
#' @param object Object of class `HousekeepingNonCodingGene`
#' (and derived classes).
#' @param value New value for essentiality_score.
#' @return The modified object.
#' @examples
#' # Since HousekeepingNonCodingGene is a virtual class, the example is made
#' # using a derived class
#' t1 <- tRNAGene(essentiality_score = 12)
#' essentiality_score(t1)
#' essentiality_score(t1) <- 20
#' essentiality_score(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("essentiality_score<-",
            "HousekeepingNonCodingGene", function(object, value) {
                object@essentiality_score <- value
                validObject(object)
                return(object) })


#' @title Set Ubiquitous Expression
#' @description Setter for ubiquitous_expression attribute of
#' `HousekeepingNonCodingGene` class (and derived classes).
#' @param object Object of class `HousekeepingNonCodingGene`
#' (and derived classes).
#' @param value New value for ubiquitous_expression.
#' @return The modified object.
#' @seealso
#' \link{ubiquitous_expression<-,HousekeepingNonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' t1 <- tRNAGene(ubiquitous_expression = FALSE)
#' ubiquitous_expression(t1)
#' ubiquitous_expression(t1) <- TRUE
#' ubiquitous_expression(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("ubiquitous_expression<-",
            function(object, value)
                standardGeneric("ubiquitous_expression<-"))


#' @title Set Ubiquitous Expression
#' @description Setter for ubiquitous_expression attribute of
#' `HousekeepingNonCodingGene` class (and derived classes).
#' @param object Object of class `HousekeepingNonCodingGene`
#' (and derived classes).
#' @param value New value for ubiquitous_expression.
#' @return The modified object.
#' @examples
#' # Since HousekeepingNonCodingGene is a virtual class, the example is made
#' # using a derived class
#' t1 <- tRNAGene(ubiquitous_expression = FALSE)
#' ubiquitous_expression(t1)
#' ubiquitous_expression(t1) <- TRUE
#' ubiquitous_expression(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("ubiquitous_expression<-",
            "HousekeepingNonCodingGene", function(object, value) {
                object@ubiquitous_expression <- value
                validObject(object)
                return(object) })
