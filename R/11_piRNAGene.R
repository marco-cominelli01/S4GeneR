#' @title piRNAGene Class
#' @description Class to represent genes encoding for piwi-interacting RNAs.
#' @slot piRNA_mature_sequence `RNAString`. piRNA mature sequence.
#' @slot piRNA_associated_PIWI_proteins `list`.
#' PIWI proteins associated to piRNA.
#' @return An object of class `piRNAGene`.
#' @seealso
#' \linkS4class{RegulatoryNonCodingGene}
#' @examples
#' # An empty piRNAGene object
#' pir1 <- new("piRNAGene")
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("piRNAGene",
        contains = "SmallNonCodingGene",
        slots = c(
            piRNA_mature_sequence = "RNAString",
            piRNA_associated_PIWI_proteins = "list"
        ),
        prototype = list(
            piRNA_mature_sequence = RNAString(),
            piRNA_associated_PIWI_proteins = list()
        ))


#' @title Validity Check for piRNAGene Class
#' @description The function checks that mature piRNA sequence is consistent
#' with the gene sequence and coordinates defined.
#' @param object Object of class `piRNAGene`.
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \link{RegulatoryNonCodingGene-validity}
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name piRNAGene-validity

setValidity("piRNAGene", function(object) {
    check_mature_length(object=object,mature_sequence='piRNA_mature_sequence')
    return(TRUE) })


# GETTERs for piRNA Gene

#' @title Get Mature piRNA Sequence
#' @description Getter for piRNA_mature_sequence attribute of `piRNAGene`
#' class.
#' @param object Object of class `piRNAGene`.
#' @return Mature piRNA sequence.
#' @examples
#' pir1 <- piRNAGene(piRNA_mature_sequence = 'AUGCA')
#' mature_sequence(pir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence", "piRNAGene", function(object) {
    object@piRNA_mature_sequence })


#' @title Get Associated PIWI Proteins
#' @description Getter for piRNA_associated_PIWI_proteins attribute of
#' `piRNAGene` class.
#' @param object Object of class `piRNAGene`.
#' @return List of associated PIWI proteins.
#' @seealso
#' \link{associated_PIWI_proteins,piRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' pir1 <- piRNAGene(piRNA_associated_PIWI_proteins=list("PIWIL1","PIWIL2"))
#' associated_PIWI_proteins(pir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("associated_PIWI_proteins",
            function(object) standardGeneric("associated_PIWI_proteins"))

#' @title Get Associated PIWI Proteins
#' @description Getter for piRNA_associated_PIWI_proteins attribute of
#' `piRNAGene` class.
#' @param object Object of class `piRNAGene`.
#' @return List of associated PIWI proteins.
#' @examples
#' pir1 <- piRNAGene(piRNA_associated_PIWI_proteins=list("PIWIL1","PIWIL2"))
#' associated_PIWI_proteins(pir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("associated_PIWI_proteins", "piRNAGene", function(object) {
    object@piRNA_associated_PIWI_proteins })


# SETTERs for piRNA Gene

#' @title Set Mature piRNA Sequence
#' @description Setter for piRNA_mature_sequence attribute of `piRNAGene`
#' class.
#' @param object Object of class `piRNAGene`.
#' @param value New value for piRNA_mature_sequence.
#' @return The modified object.
#' @examples
#' pir1 <- piRNAGene(piRNA_mature_sequence = 'AUUA')
#' mature_sequence(pir1)
#' mature_sequence(pir1) <- 'AAAGCA'
#' mature_sequence(pir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence<-", "piRNAGene", function(object, value) {
    object@piRNA_mature_sequence <- RNAString(as.character(value))
    validObject(object)
    return(object) })


#' @title Add/Remove Associated PIWI Proteins
#' @description Setter for piRNA_associated_PIWI_proteins attribute of
#' `piRNAGene` class. By choosing an action, which can be
#' 'add' or 'remove', it's possible to add and remove, respectively,
#' PIWI proteins from piRNA_associated_PIWI_proteins attribute.
#' @param object Object of class `piRNAGene`.
#' @param action Action to be done: can be 'add' or 'remove'.
#' @param value PIWI protein to add/remove.
#' @return The modified object.
#' @seealso
#' \link{associated_PIWI_proteins<-,piRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' pir1 <- piRNAGene()
#' associated_PIWI_proteins(pir1)
#' associated_PIWI_proteins(pir1, action = 'add') <- 'PIWIL2'
#' associated_PIWI_proteins(pir1)
#' associated_PIWI_proteins(pir1, action = 'remove') <- 'PIWIL2'
#' associated_PIWI_proteins(pir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("associated_PIWI_proteins<-",
            function(object, action = 'add', value)
                standardGeneric("associated_PIWI_proteins<-"))


#' @title Add/Remove Associated PIWI Proteins
#' @description Setter for piRNA_associated_PIWI_proteins attribute of
#' `piRNAGene` class. By choosing an action, which can be
#' 'add' or 'remove', it's possible to add and remove, respectively,
#' PIWI proteins from piRNA_associated_PIWI_proteins attribute.
#' @param object Object of class `piRNAGene`.
#' @param action Action to be done: can be 'add' or 'remove'.
#' @param value PIWI protein to add/remove.
#' @return The modified object.
#' @examples
#' pir1 <- piRNAGene()
#' associated_PIWI_proteins(pir1)
#' associated_PIWI_proteins(pir1, action = 'add') <- 'PIWIL2'
#' associated_PIWI_proteins(pir1)
#' associated_PIWI_proteins(pir1, action = 'remove') <- 'PIWIL2'
#' associated_PIWI_proteins(pir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("associated_PIWI_proteins<-", "piRNAGene",
        function(object, action, value) {
            msg <- paste0("e.g. 'associated_PIWI_protein(gene1,",
                        " action = 'remove') <- 'PIWIL1' '.")
            object@piRNA_associated_PIWI_proteins <- modify_list(
                object = object,
                slot = "piRNA_associated_PIWI_proteins",
                action = action,
                value = value,
                msg = msg)
            validObject(object)
            return(object) })


#' @title Create an Object of Class piRNAGene
#' @usage
#' piRNAGene(
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
#'     piRNA_mature_sequence = RNAString(),
#'     piRNA_associated_PIWI_proteins = list() )
#' @description Constructor function for `piRNAGene` class.
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
#' @param piRNA_mature_sequence Mature piRNA sequence.
#' @param piRNA_associated_PIWI_proteins piRNA associated PIWI proteins.
#' @return New piRNAGene object.
#' @examples
#' # An empty piRNAGene object
#' pir1 <- piRNAGene()
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

piRNAGene <- function(gene_ensembl_id = "ENSG00000000000",
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
                    piRNA_mature_sequence = RNAString(),
                    piRNA_associated_PIWI_proteins = list()) {

    coordinates <- constructor_helper(chrom = chromosome, str = strand,
                                    g_st = gene_start, g_en = gene_end,
                                    ex_st = exons_starts, ex_en=exons_ends,
                                    pre_st = pre_ncRNA_start,
                                    pre_en = pre_ncRNA_end,
                                    coding = FALSE)
    new("piRNAGene",
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
        piRNA_mature_sequence = RNAString(piRNA_mature_sequence),
        piRNA_associated_PIWI_proteins = piRNA_associated_PIWI_proteins) }
