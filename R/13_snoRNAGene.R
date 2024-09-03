#' @title snoRNAGene Class
#' @description Class to represent genes encoding for small nucleolar RNAs.
#' @slot snoRNA_mature_sequence `RNAString`. snoRNA mature sequence.
#' @slot snoRNA_box_type `character`. snoRNA's box type.
#' @return An object of class `snoRNAGene`.
#' @seealso
#' \linkS4class{RegulatoryNonCodingGene}
#' @examples
#' # An empty snoRNAGene object
#' sno1 <- new("snoRNAGene")
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("snoRNAGene",
        contains = "SmallNonCodingGene",
        slots = c(
            snoRNA_mature_sequence = "RNAString",
            snoRNA_box_type = "character"
        ),
        prototype = list(
            snoRNA_mature_sequence = RNAString(),
            snoRNA_box_type = NA_character_
        ))


#' @title Validity Check for snoRNAGene Class
#' @description The function checks that: each slots has a length of 1 (or
#' empty), the inserted snoRNA box type is valid and the mature sequence
#' is consistent with the gene sequence and the defined coordinates.
#' @param object Object of class `snoRNAGene`.
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \link{RegulatoryNonCodingGene-validity}
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name snoRNAGene-validity

setValidity("snoRNAGene", function(object) {
    check_slots_length(obj = object, class_obj = "snoRNAGene",
                        mother_class_obj = "SmallNonCodingGene",
                        special_slots = "snoRNA_mature_sequence")

    check_mature_length(object = object,
                        mature_sequence = 'snoRNA_mature_sequence')

    snoRNA_box_types <- c("C/D","H/ACA")
    if (!object@snoRNA_box_type %in% c(snoRNA_box_types, NA)) {
        err_msg <- paste0("Invalid snoRNA box type. Please choose among",
                        " {", paste(snoRNA_box_types, collapse = ', '), '}.')
        stop(err_msg) }

    return(TRUE) })


# GETTERs for snoRNA Gene

#' @title Get Mature snoRNA Sequence
#' @description Getter for snoRNA_mature_sequence attribute of `snoRNAGene`
#' class.
#' @param object Object of class `snoRNAGene`.
#' @return Mature snoRNA sequence.
#' @examples
#' sno1 <- snoRNAGene(snoRNA_mature_sequence = 'AUGCA')
#' mature_sequence(sno1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence", "snoRNAGene", function(object) {
    object@snoRNA_mature_sequence })


#' @title Get snoRNA Box Type
#' @description Getter for snoRNA_box_type attribute of `snoRNAGene` class.
#' @param object Object of class `snoRNAGene`.
#' @return snoRNA box type.
#' @seealso
#' \link{box_type,snoRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' sno1 <- snoRNAGene(snoRNA_box_type = 'C/D')
#' box_type(sno1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("box_type",
            function(object) standardGeneric("box_type"))


#' @title Get snoRNA Box Type
#' @description Getter for snoRNA_box_type attribute of `snoRNAGene` class.
#' @param object Object of class `snoRNAGene`.
#' @return snoRNA box type.
#' @examples
#' sno1 <- snoRNAGene(snoRNA_box_type = 'C/D')
#' box_type(sno1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("box_type", "snoRNAGene", function(object) {
    object@snoRNA_box_type })


# SETTERs for snoRNA Gene

#' @title Set Mature snoRNA Sequence
#' @description Setter for snoRNA_mature_sequence attribute of
#' `snoRNAGene` class.
#' @param object Object of class `snoRNAGene`.
#' @param value New value for snoRNA_mature_sequence.
#' @return The modified object.
#' @examples
#' sno1 <- snoRNAGene(snoRNA_mature_sequence = 'AUUA')
#' mature_sequence(sno1)
#' mature_sequence(sno1) <- 'AAAGCA'
#' mature_sequence(sno1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence<-", "snoRNAGene", function(object, value) {
    object@snoRNA_mature_sequence <- RNAString(as.character(value))
    validObject(object)
    return(object) })


#' @title Set snoRNA Box Type
#' @description Setter for snoRNA_box_type attribute of `snoRNAGene` class.
#' @param object Object of class `snoRNAGene`.
#' @param value New value for snoRNA_box_type.
#' @return The modified object.
#' @seealso
#' \link{box_type<-,snoRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' sno1 <- snoRNAGene(snoRNA_box_type = 'C/D')
#' box_type(sno1)
#' box_type(sno1) <- 'H/ACA'
#' box_type(sno1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("box_type<-",
            function(object, value) standardGeneric("box_type<-"))

#' @title Set snoRNA Box Type
#' @description Setter for snoRNA_box_type attribute of `snoRNAGene` class.
#' @param object Object of class `snoRNAGene`.
#' @param value New value for snoRNA_box_type.
#' @return The modified object.
#' @examples
#' sno1 <- snoRNAGene(snoRNA_box_type = 'C/D')
#' box_type(sno1)
#' box_type(sno1) <- 'H/ACA'
#' box_type(sno1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("box_type<-", "snoRNAGene", function(object, value) {
    object@snoRNA_box_type <- value
    validObject(object)
    return(object) })

#' @title Create an Object of Class snoRNAGene
#' @usage
#' snoRNAGene(
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
#'     snoRNA_mature_sequence = RNAString(),
#'     snoRNA_box_type = NA_character_ )
#' @description Constructor function for `snoRNAGene` class.
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
#' @param snoRNA_mature_sequence Mature snoRNA sequence.
#' @param snoRNA_box_type snoRNA box type.
#' @return New snoRNAGene object.
#' @examples
#' # An empty snoRNAGene object
#' sno1 <- snoRNAGene()
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

snoRNAGene <- function(gene_ensembl_id = "ENSG00000000000",
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
                        snoRNA_mature_sequence = RNAString(),
                        snoRNA_box_type = NA_character_) {

    coordinates <- constructor_helper(chrom = chromosome, str = strand,
                                    g_st = gene_start, g_en = gene_end,
                                    ex_st = exons_starts, ex_en=exons_ends,
                                    pre_st = pre_ncRNA_start,
                                    pre_en = pre_ncRNA_end,
                                    coding = FALSE)

    new("snoRNAGene",
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
        snoRNA_mature_sequence = RNAString(snoRNA_mature_sequence),
        snoRNA_box_type = snoRNA_box_type) }
