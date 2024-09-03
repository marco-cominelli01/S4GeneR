#' @title RegulatoryNonCodingGene Virtual Class
#' @description This virtual class is the class from which
#' SmallNonCodingGene and LongNonCodingGene classes are derived.
#' @slot regulatory_functions `list`.
#' Regulatory functions of the gene.
#' @slot targets_ensembl_id `list`. ENSEMBL IDs of the gene targets.
#' @return No objects can be created directly from this class since
#' it's virtual.
#' @seealso
#' \linkS4class{NonCodingGene}, \linkS4class{SmallNonCodingGene},
#' \linkS4class{LongNonCodingGene}
#' @examples
#' # The RegulatoryNonCodingGene class is virtual, so no objects can
#' # be created directly from this class. It needs to be extended by
#' # other classes.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("RegulatoryNonCodingGene",
        contains = c("NonCodingGene", "VIRTUAL"),
        slots = c(
            regulatory_functions = "list",
            targets_ensembl_id = "list"
        ),
        prototype = list(
            regulatory_functions = list(),
            targets_ensembl_id = list()
        ))


#' @title Validity Check for RegulatoryNonCodingGene Class
#' @description The function checks that each ID in targets_ensembl_id slot
#' is a valid ENSEMBL ID.
#' @param object Object of class `RegulatoryNonCodingGene`
#' (and derived classes).
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \link{NonCodingGene-validity}
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name RegulatoryNonCodingGene-validity

setValidity("RegulatoryNonCodingGene", function(object) {
    ### Check the validity of all IDs in 'targets_ensembl_id' slot.
    if (length(object@targets_ensembl_id)) {
        for (i in seq_len(length(object@targets_ensembl_id))) {
            id_to_check <- object@targets_ensembl_id[[i]]
            prefixes <- c("ENSG", "ENST", "ENSP")
            counter <- 0
            for (j in seq_len(3)) {
                prefix <- prefixes[j]
                result <- try(check_ensembl_id(id_to_check, prefix, ""),
                                silent = TRUE)
                if (inherits(result, 'try-error')) {
                    counter <- counter + 1 }
            }
            if (counter == 3) {
                err_msg <- paste0("The ID in slot 'targets_ensembl_id'",
                            " is not a valid gene/transcript/",
                            "protein ENSEMBL ID. The ID must start with ",
                            "ENSG/ENST/ENSP followed by exactly 11 digits.")
                stop(err_msg)  }
        }
    }
    return(TRUE) })


# GETTERs for RegulatoryNonCodingGene

#' @title Get Regulatory Functions
#' @description Getter for regulatory_functions attribute of
#' `RegulatoryNonCodingGene` class (and derived classes).
#' @param object Object of class `RegulatoryNonCodingGene`
#' (and derived classes).
#' @return List of regulatory functions.
#' @seealso
#' \link{regulatory_functions,RegulatoryNonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' mir1 <- miRNAGene(regulatory_functions = list("Trascription modulation"))
#' regulatory_functions(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("regulatory_functions",
            function(object) standardGeneric("regulatory_functions"))


#' @title Get Regulatory Functions
#' @description Getter for regulatory_functions attribute of
#' `RegulatoryNonCodingGene` class (and derived classes).
#' @param object Object of class `RegulatoryNonCodingGene`
#' (and derived classes).
#' @return List of regulatory functions.
#' @examples
#' # Since RegulatoryNonCodingGene is a virtual class, the example is made
#' # using a derived class
#' mir1 <- miRNAGene(regulatory_functions = list("Transcription modulation"))
#' regulatory_functions(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("regulatory_functions", "RegulatoryNonCodingGene", function(object){
    object@regulatory_functions })


#' @title Get Targets
#' @description Getter for targets_ensembl_id attribute of
#' `RegulatoryNonCodingGene` class (and derived classes).
#' @param object Object of class `RegulatoryNonCodingGene`
#' (or derived classes).
#' @return List of targets ENSEMBL IDs.
#' @seealso
#' \link{targets_id,RegulatoryNonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' mir1 <- miRNAGene(targets_ensembl_id = list("ENST12341234090"))
#' targets_id(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("targets_id", function(object) standardGeneric("targets_id"))


#' @title Get Targets
#' @description Getter for targets_ensembl_id attribute of
#' `RegulatoryNonCodingGene` class (and derived classes).
#' @param object Object of class `RegulatoryNonCodingGene`
#' (and derived classes).
#' @return List of targets ENSEMBL IDs.
#' @examples
#' # Since RegulatoryNonCodingGene is a virtual class, the example is made
#' # using a derived class
#' mir1 <- miRNAGene(targets_ensembl_id = list("ENST12341234090"))
#' targets_id(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("targets_id", "RegulatoryNonCodingGene", function(object) {
    object@targets_ensembl_id })


# SETTERs for RegulatoryNonCodingGene

#' @title Add/Remove Regulatory Functions
#' @description Setter for regulatory_functions attribute of
#' `RegulatoryNonCodingGene` class (and derived classes). By choosing an
#' action, which can be 'add' or 'remove', it's possible to add and remove,
#' respectively, regulatory function from regulatory_functions attribute.
#' @param object Object of class `RegulatoryNonCodingGene`
#' (and derived classes).
#' @param action Action to be done: can be 'add' or 'remove'.
#' @param value Regulatory function to add/remove.
#' @return The modified object.
#' @seealso
#' \link{regulatory_functions<-,RegulatoryNonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' mir1 <- miRNAGene()
#' regulatory_functions(mir1)
#' regulatory_functions(mir1, action = 'add') <- 'transcription factor'
#' regulatory_functions(mir1)
#' regulatory_functions(mir1, action = 'remove') <- 'transcription factor'
#' regulatory_functions(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("regulatory_functions<-",
            function(object, action = 'add', value)
                standardGeneric("regulatory_functions<-"))


#' @title Add/Remove Regulatory Functions
#' @usage
#' \S4method{regulatory_functions}{RegulatoryNonCodingGene}(object,
#' action) <- value
#' @description Setter for regulatory_functions attribute of
#' RegulatoryNonCodingGene class (and derived classes). By choosing an
#' action, which can be 'add' or 'remove', it's possible to add and remove,
#' respectively, regulatory function from regulatory_functions attribute.
#' @param object Object of class `RegulatoryNonCodingGene`
#' (and derived classes).
#' @param action Action to be done: can be 'add' or 'remove'.
#' @param value Regulatory function to add/remove.
#' @return The modified object.
#' @examples
#' # Since RegulatoryNonCodingGene is a virtual class, the example is made
#' # using a derived class
#' mir1 <- miRNAGene()
#' regulatory_functions(mir1)
#' regulatory_functions(mir1, action = 'add') <- 'transcription factor'
#' regulatory_functions(mir1)
#' regulatory_functions(mir1, action = 'remove') <- 'transcription factor'
#' regulatory_functions(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("regulatory_functions<-", "RegulatoryNonCodingGene",
            function(object, action, value) {
                msg <- paste0("e.g. 'regulatory_functions(gene1, action =",
                            "  'remove') <- 'chromatin remodeling' '.")
                object@regulatory_functions <- modify_list(object = object,
                                                slot = "regulatory_functions",
                                                action = action,
                                                value = value,
                                                msg = msg)
                validObject(object)
                return(object) })


#' @title Add/Remove Targets
#' @description Setter for targets_ensembl_id attribute of
#' `RegulatoryNonCodingGene` class (and derived classes). By choosing an
#' action, which can be 'add' or 'remove', it's possible to add and remove,
#' respectively, targets ENSEMBL IDs from targets_ensembl_id attribute.
#' @param object Object of class `RegulatoryNonCodingGene`
#' (and derived classes).
#' @param action Action to be done: can be 'add' or 'remove'.
#' @param value Target ENSEMBL ID to add/remove.
#' @return The modified object.
#' @seealso
#' \link{targets_id<-,RegulatoryNonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' mir1 <- miRNAGene()
#' targets_id(mir1)
#' targets_id(mir1, action = 'add') <- 'ENST09090909091'
#' targets_id(mir1)
#' targets_id(mir1, action = 'remove') <- 'ENST09090909091'
#' targets_id(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("targets_id<-",
            function(object, action = 'add', value)
                standardGeneric("targets_id<-"))


#' @title Add/Remove Targets
#' @description Setter for targets_ensembl_id attribute of
#' `RegulatoryNonCodingGene` class (and derived classes). By choosing an
#' action, which can be 'add' or 'remove', it's possible to add and remove,
#' respectively, targets ENSEMBL IDs from targets_ensembl_id attribute.
#' @param object Object of class `RegulatoryNonCodingGene`
#' (and derived classes).
#' @param action Action to be done: can be 'add' or 'remove'.
#' @param value Target ENSEMBL ID to add/remove.
#' @return The modified object.
#' @examples
#' # Since RegulatoryNonCodingGene is a virtual class, the example is made
#' # using a derived class
#' mir1 <- miRNAGene()
#' targets_id(mir1)
#' targets_id(mir1, action = 'add') <- 'ENST09090909091'
#' targets_id(mir1)
#' targets_id(mir1, action = 'remove') <- 'ENST09090909091'
#' targets_id(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("targets_id<-", "RegulatoryNonCodingGene",
        function(object, action, value) {
            msg <- paste0("e.g. 'targets_id(gene1, action = 'remove') <- ",
                        "'chromatin remodeling' '.")
            object@targets_ensembl_id <- modify_list(object = object,
                                                slot = "targets_ensembl_id",
                                                action = action,
                                                value = value,
                                                msg = msg)
            validObject(object)
            return(object) })
