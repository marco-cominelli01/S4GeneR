#' @title siRNAGene Class
#' @description Class to represent genes encoding for short interfering RNAs.
#' @slot siRNA_mature_sequence `RNAString`. siRNA mature sequence.
#' @slot siRNA_off_targets `list`. siRNA off-targets.
#' @return An object of class `siRNAGene`.
#' @seealso
#' \linkS4class{RegulatoryNonCodingGene}
#' @examples
#' # An empty siRNAGene object
#' sir1 <- new("siRNAGene")
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("siRNAGene",
        contains = "SmallNonCodingGene",
        slots = c(
            siRNA_mature_sequence = "RNAString",
            siRNA_off_targets = "list"
        ),
        prototype = list(
            siRNA_mature_sequence = RNAString(),
            siRNA_off_targets = list()
        ))


#' @title Validity Check for siRNAGene Class
#' @description The function checks that mature siRNA sequence is consistent
#' with the gene sequence and coordinates defined.
#' @param object Object of class `siRNAGene`.
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \link{RegulatoryNonCodingGene-validity}
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name siRNAGene-validity

setValidity("siRNAGene", function(object) {

    # Check on length of siRNA mature sequence w.r.t.
    # gene/exons/pre-ncRNA coordinates
    check_mature_length(object = object,
                        mature_sequence='siRNA_mature_sequence')

    return(TRUE) })


# GETTERs for siRNA Gene

#' @title Get Mature siRNA sequence
#' @description Getter for siRNA_mature_sequence attribute of `siRNAGene`
#' class.
#' @param object Object of class `siRNAGene`.
#' @return Mature siRNA sequence.
#' @examples
#' sir1 <- siRNAGene(siRNA_mature_sequence = 'AUGCA')
#' mature_sequence(sir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence", "siRNAGene", function(object) {
    object@siRNA_mature_sequence })


#' @title Get Off-Targets
#' @description Getter for siRNA_off_targets attribute of `siRNAGene` class.
#' @param object Object of class `siRNAGene`.
#' @return List of off-targets.
#' @seealso
#' \link{off_targets,siRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' sir1 <- siRNAGene(siRNA_off_targets = list("BRCA1","BRCA2"))
#' off_targets(sir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("off_targets",
            function(object) standardGeneric("off_targets"))


#' @title Get Off-Targets
#' @description Getter for siRNA_off_targets attribute of `siRNAGene` class.
#' @param object Object of class `siRNAGene`.
#' @return List of off-targets.
#' @examples
#' sir1 <- siRNAGene(siRNA_off_targets = list("BRCA1","BRCA2"))
#' off_targets(sir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("off_targets", "siRNAGene", function(object) {
    object@siRNA_off_targets })

# SETTERs for siRNA Gene

#' @title Set Mature siRNA Sequence
#' @description Setter for siRNA_mature_sequence attribute of `siRNAGene`
#' class.
#' @param object Object of class `siRNAGene`.
#' @param value New value for siRNA_mature_sequence.
#' @return The modified object.
#' @examples
#' sir1 <- siRNAGene(siRNA_mature_sequence = 'AUUA')
#' mature_sequence(sir1)
#' mature_sequence(sir1) <- 'AAAGCA'
#' mature_sequence(sir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence<-", "siRNAGene", function(object, value) {
    object@siRNA_mature_sequence <- RNAString(as.character(value))
    validObject(object)
    return(object) })


#' @title Add/Remove Off-Targets
#' @description Setter for siRNA_off_targets attribute of
#' `siRNAGene` class. By choosing an action, which can be
#' 'add' or 'remove', it's possible to add and remove, respectively,
#' off-targets from siRNA_off_targets attribute.
#' @param object Object of class `siRNAGene`.
#' @param action Action to be done: can be 'add' or 'remove'.
#' @param value Off-target to add/remove.
#' @return The modified object.
#' @seealso
#' \link{off_targets<-,siRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' sir1 <- siRNAGene()
#' off_targets(sir1)
#' off_targets(sir1, action = 'add') <- 'IRF2'
#' off_targets(sir1)
#' off_targets(sir1, action = 'remove') <- 'IRF2'
#' off_targets(sir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("off_targets<-",
            function(object, action = 'add', value)
                standardGeneric("off_targets<-"))


#' @title Add/Remove Off-Targets
#' @description Setter for siRNA_off_targets attribute of
#' `siRNAGene` class. By choosing an action, which can be
#' 'add' or 'remove', it's possible to add and remove, respectively,
#' off-targets from siRNA_off_targets attribute.
#' @param object Object of class `siRNAGene`.
#' @param action Action to be done: can be 'add' or 'remove'.
#' @param value Off-target to add/remove.
#' @return The modified object.
#' @examples
#' sir1 <- siRNAGene()
#' off_targets(sir1)
#' off_targets(sir1, action = 'add') <- 'IRF2'
#' off_targets(sir1)
#' off_targets(sir1, action = 'remove') <- 'IRF2'
#' off_targets(sir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("off_targets<-", "siRNAGene",
        function(object, action, value) {
            msg <- paste0("e.g. 'siRNA_off_targets(gene1, action =",
                            " 'remove') <- 'IRF2' '.")
            object@siRNA_off_targets <- modify_list(object = object,
                                                slot = "siRNA_off_targets",
                                                action = action,
                                                value = value,
                                                msg = msg)
            validObject(object)
            return(object) })


#' @title Create an Object of Class siRNAGene
#' @usage
#' siRNAGene(
#'     gene_ensembl_id = "ENSG00000000000",
#'     hugo_symbol = NA_character_,
#'     gene_complete_name = NA_character_,
#'     gene_description = NA_character_,
#'     strand = NA,
#'     chromosome = NA,
#'     gene_start = NA,
#'     gene_end = NA,
#'     gene_sequence = DNAString(),
#'     exons_starts = NA,
#'     exons_ends = NA,
#'     pre_ncRNA_ensembl_id = "ENST00000000000",
#'     pre_ncRNA_start = NA,
#'     pre_ncRNA_end = NA,
#'     regulatory_functions = list(),
#'     targets_ensembl_id = list(),
#'     siRNA_mature_sequence = RNAString(),
#'     siRNA_off_targets = list() )
#' @description Constructor function for `siRNAGene` class.
#' @param gene_ensembl_id Gene ENSEMBL ID.
#' @param hugo_symbol Gene Hugo symbol.
#' @param gene_complete_name Gene complete name.
#' @param gene_description Gene description.
#' @param strand Strand on which the gene is annotated.
#' @param chromosome Chromosome on which the gene is annotated.
#' @param gene_start Gene start coordinate.
#' @param gene_end Gene end coordinate.
#' @param gene_sequence Gene sequence.
#' @param exons_starts Exons starts coordinates.
#' @param exons_ends Exons ends coordinates.
#' @param pre_ncRNA_ensembl_id pre-ncRNA ENSEMBL ID.
#' @param pre_ncRNA_start pre-ncRNA start coordinate.
#' @param pre_ncRNA_end pre-ncRNA end coordinate.
#' @param regulatory_functions Gene regulatory functions.
#' @param targets_ensembl_id ENSEMBL IDs of the gene targets.
#' @param siRNA_mature_sequence Mature siRNA sequence.
#' @param siRNA_off_targets siRNA off-targets.
#' @return New siRNAGene object.
#' @examples
#' # An empty siRNAGene object
#' sir1 <- siRNAGene()
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

siRNAGene <- function(gene_ensembl_id = "ENSG00000000000",
                    hugo_symbol = NA_character_,
                    gene_complete_name = NA_character_,
                    gene_description = NA_character_,
                    strand = NA, chromosome = NA,
                    gene_start = NA, gene_end = NA,
                    gene_sequence = DNAString(),
                    exons_starts = NA, exons_ends = NA,
                    pre_ncRNA_ensembl_id = "ENST00000000000",
                    pre_ncRNA_start = NA, pre_ncRNA_end = NA,
                    regulatory_functions = list(),
                    targets_ensembl_id = list(),
                    siRNA_mature_sequence = RNAString(),
                    siRNA_off_targets = list()) {

    coordinates <- constructor_helper(chrom = chromosome, str = strand,
                                    g_st = gene_start, g_en = gene_end,
                                    ex_st = exons_starts, ex_en=exons_ends,
                                    pre_st = pre_ncRNA_start,
                                    pre_en = pre_ncRNA_end,
                                    coding = FALSE)
    new("siRNAGene",
        gene_ensembl_id = gene_ensembl_id,
        hugo_symbol = hugo_symbol,
        gene_complete_name = gene_complete_name,
        gene_description = gene_description,
        gene_coordinates = coordinates[[1]],
        gene_sequence = DNAString(gene_sequence),
        exons_coordinates = coordinates[[2]],
        pre_ncRNA_ensembl_id = pre_ncRNA_ensembl_id,
        pre_ncRNA_coordinates = coordinates[[3]],
        regulatory_functions = regulatory_functions,
        targets_ensembl_id = targets_ensembl_id,
        siRNA_mature_sequence = RNAString(siRNA_mature_sequence),
        siRNA_off_targets = siRNA_off_targets) }
