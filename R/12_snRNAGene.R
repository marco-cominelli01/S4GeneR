#' @title snRNAGene Class
#' @description Class to represent genes encoding for small nuclear RNAs.
#' @slot snRNA_mature_sequence `RNAString`. snRNA mature sequence.
#' @slot spliceosome_complex `character`.
#' Spliceosome complex to which snRNA is associated.
#' @return An object of class `snRNAGene`.
#' @seealso
#' \linkS4class{RegulatoryNonCodingGene}
#' @examples
#' # An empty snRNAGene object
#' sn1 <- new("snRNAGene")
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("snRNAGene",
        contains = "SmallNonCodingGene",
        slots = c(
            snRNA_mature_sequence = "RNAString",
            spliceosome_complex = "character"
        ),
        prototype = list(
            snRNA_mature_sequence = RNAString(),
            spliceosome_complex = NA_character_
        ))


#' @title Validity Check for snRNAGene Class
#' @description The function checks that: each slots has a length of 1 (or
#' empty), the inserted spliceosome complex is valid and the mature sequence
#' is consistent with the gene sequence and the defined coordinates.
#' @param object Object of class `snRNAGene`.
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \link{RegulatoryNonCodingGene-validity}
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name snRNAGene-validity

setValidity("snRNAGene", function(object) {
    check_slots_length(obj = object, class_obj = "snRNAGene",
                        mother_class_obj = "SmallNonCodingGene",
                        special_slots = "snRNA_mature_sequence")

    spliceosome_complex_types <- c("U1","U2","U4","U5","U6",
                                    "not associated to spliceosome")
    if (!is.na(object@spliceosome_complex) &&
        !object@spliceosome_complex %in% spliceosome_complex_types) {
            err_msg <- paste0("Please insert a valid type for the",
                            " spliceosome complex. Choose among {",
                            paste(spliceosome_complex_types, collapse = ', '),
                            "}")
            stop(err_msg) }

    check_mature_length(object=object,mature_sequence='snRNA_mature_sequence')

    return(TRUE) })


# GETTERs for snRNA Gene

#' @title Get Mature snRNA Sequence
#' @description Getter for snRNA_mature_sequence attribute of `snRNAGene`
#' class.
#' @param object Object of class `snRNAGene`.
#' @return Mature snRNA sequence.
#' @examples
#' sn1 <- snRNAGene(snRNA_mature_sequence = 'AUGCA')
#' mature_sequence(sn1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence", "snRNAGene", function(object) {
    object@snRNA_mature_sequence })

#' @title Get snRNA Spliceosome Complex
#' @description Getter for spliceosome_complex attribute of `snRNAGene` class.
#' @param object Object of class `snRNAGene`.
#' @return Mature snRNA sequence.
#' @seealso
#' \link{spliceosome_complex,snRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' sn1 <- snRNAGene(spliceosome_complex = 'U4')
#' spliceosome_complex(sn1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("spliceosome_complex",
            function(object) standardGeneric("spliceosome_complex"))

#' @title Get snRNA Spliceosome Complex
#' @description Getter for spliceosome_complex attribute of `snRNAGene` class.
#' @param object Object of class `snRNAGene`.
#' @return Mature snRNA sequence.
#' @examples
#' sn1 <- snRNAGene(spliceosome_complex = 'U4')
#' spliceosome_complex(sn1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("spliceosome_complex", "snRNAGene", function(object) {
    object@spliceosome_complex })


# SETTERs for snRNA Gene

#' @title Set Mature snRNA Sequence
#' @description Setter for snRNA_mature_sequence attribute of `snRNAGene`
#' class.
#' @param object Object of class `snRNAGene`.
#' @param value New value for snRNA_mature_sequence.
#' @return The modified object.
#' @examples
#' sn1 <- snRNAGene(snRNA_mature_sequence = 'AUUA')
#' mature_sequence(sn1)
#' mature_sequence(sn1) <- 'AAAGCA'
#' mature_sequence(sn1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence<-", "snRNAGene", function(object, value) {
    object@snRNA_mature_sequence <- RNAString(as.character(value))
    validObject(object)
    return(object) })


#' @title Set snRNA Spliceosome Complex
#' @description Setter for spliceosome_complex attribute of `snRNAGene` class.
#' @param object Object of class `snRNAGene`.
#' @param value New value for spliceosome_complex.
#' @return The modified object.
#' @seealso
#' \link{spliceosome_complex<-,snRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' sn1 <- snRNAGene(spliceosome_complex = 'U4')
#' spliceosome_complex(sn1)
#' spliceosome_complex(sn1) <- 'U6'
#' spliceosome_complex(sn1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("spliceosome_complex<-",
            function(object, value) standardGeneric("spliceosome_complex<-"))


#' @title Set snRNA Spliceosome Complex
#' @description Setter for spliceosome_complex attribute of `snRNAGene` class.
#' @param object Object of class `snRNAGene`.
#' @param value New value for spliceosome_complex.
#' @return The modified object.
#' @examples
#' sn1 <- snRNAGene(spliceosome_complex = 'U4')
#' spliceosome_complex(sn1)
#' spliceosome_complex(sn1) <- 'U6'
#' spliceosome_complex(sn1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("spliceosome_complex<-", "snRNAGene", function(object, value) {
    object@spliceosome_complex <- value
    validObject(object)
    return(object) })


#' @title Create an Object of Class snRNAGene
#' @usage
#' snRNAGene(
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
#'     snRNA_mature_sequence = RNAString(),
#'     spliceosome_complex = NA_character_ )
#' @description Constructor function for `snRNAGene` class.
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
#' @param snRNA_mature_sequence Mature snRNA sequence.
#' @param spliceosome_complex Spliceosome complex to which snRNA is
#' associated.
#' @return New snRNAGene object.
#' @examples
#' # An empty snRNAGene object
#' sn1 <- snRNAGene()
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

snRNAGene <- function(gene_ensembl_id = "ENSG00000000000",
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
                    snRNA_mature_sequence = RNAString(),
                    spliceosome_complex = NA_character_) {

    coordinates <- constructor_helper(chrom = chromosome, str = strand,
                                    g_st = gene_start, g_en = gene_end,
                                    ex_st = exons_starts, ex_en=exons_ends,
                                    pre_st = pre_ncRNA_start,
                                    pre_en = pre_ncRNA_end,
                                    coding = FALSE)
    new("snRNAGene",
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
        snRNA_mature_sequence = RNAString(snRNA_mature_sequence),
        spliceosome_complex = spliceosome_complex) }
