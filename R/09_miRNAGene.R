#' @title miRNAGene Class
#' @description Class to represent genes encoding for micro RNAs.
#' @slot miRNA_mature_sequence `RNAString`. miRNA mature sequence.
#' @slot seed_sequence `RNAString`. miRNA seed sequence.
#' @return An object of class `miRNAGene`.
#' @seealso
#' \linkS4class{RegulatoryNonCodingGene}
#' @examples
#' # An empty miRNAGene object
#' mir1 <- new("miRNAGene")
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("miRNAGene",
        contains = "SmallNonCodingGene",
        slots = c(
            miRNA_mature_sequence = "RNAString",
            seed_sequence = "RNAString"
        ),
        prototype = list(
            miRNA_mature_sequence = RNAString(),
            seed_sequence = RNAString()
        ))


#' @title Validity Check for miRNAGene Class
#' @description The function checks that seed sequence is actually contained
#' in the mature miRNA sequence and that mature miRNA sequence is consistent
#' with the gene sequence and coordinates defined.
#' @param object Object of class `miRNAGene`.
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \link{RegulatoryNonCodingGene-validity}
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name miRNAGene-validity

setValidity("miRNAGene", function(object) {
    # Check that seed sequence is actually contained in the mature sequence.
    if (length(object@miRNA_mature_sequence)
        && length(object@seed_sequence)) {
        if (length(object@seed_sequence) >
            length(object@miRNA_mature_sequence)) {
            err_msg <- paste("miRNA seed sequence can't be longer",
                            "than the mature miRNA sequence")
            stop(err_msg)  }

        seed_in_transcript <- Biostrings::matchPattern(object@seed_sequence,
                                            object@miRNA_mature_sequence)
        if (length(seed_in_transcript) == 0) {
            err_msg <- paste0("The provided seed sequence is not",
                        " compatible with the provided mature miRNA ",
                        "sequence. The former should be found in the latter.")
            stop(err_msg) } }

    # This is done if the miRNA mature sequence is not defined.
    if (length(object@gene_sequence) > 0 &&
        (length(object@gene_sequence) < length(object@seed_sequence))) {
            err_msg <- paste0("Gene sequence and miRNA seed sequence",
                        " are incompatible. The latter can't be longer than",
                        " the former.")
            stop(err_msg) }

    # Check on length of miRNA mature sequence and seed sequence w.r.t.
    # gene/exons/pre-ncRNA coordinates
    check_mature_length(object=object,mature_sequence='miRNA_mature_sequence')
    check_mature_length(object=object, mature_sequence = 'seed_sequence')

    return(TRUE) })


# GETTERs for miRNA Gene

#' @title Get Mature miRNA Sequence
#' @description Getter for miRNA_mature_sequence attribute of `miRNAGene`
#' class.
#' @param object Object of class `miRNAGene`.
#' @return Mature miRNA sequence.
#' @examples
#' mir1 <- miRNAGene(miRNA_mature_sequence = 'AUGCA')
#' mature_sequence(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence", "miRNAGene", function(object) {
    object@miRNA_mature_sequence })


#' @title Get miRNA Seed Sequence
#' @description Getter for seed_sequence attribute of `miRNAGene` class.
#' @param object Object of class `miRNAGene`.
#' @return Mature miRNA sequence.
#' @seealso
#' \link{seed_sequence,miRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' mir1 <- miRNAGene(seed_sequence = 'AUGCA')
#' seed_sequence(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("seed_sequence",
            function(object) standardGeneric("seed_sequence"))


#' @title Get miRNA Seed Sequence
#' @description Getter for seed_sequence attribute of `miRNAGene` class.
#' @param object Object of class `miRNAGene`.
#' @return Mature miRNA sequence.
#' @examples
#' mir1 <- miRNAGene(seed_sequence = 'AUGCA')
#' seed_sequence(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("seed_sequence", "miRNAGene", function(object) {
    object@seed_sequence })


# SETTERs for miRNA Gene

#' @title Set Mature miRNA Sequence
#' @description Setter for miRNA_mature_sequence attribute of `miRNAGene`
#' class.
#' @param object Object of class `miRNAGene`.
#' @param value New value for miRNA_mature_sequence.
#' @return The modified object.
#' @examples
#' mir1 <- miRNAGene(miRNA_mature_sequence = 'AUUA')
#' mature_sequence(mir1)
#' mature_sequence(mir1) <- 'AAAGCA'
#' mature_sequence(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence<-", "miRNAGene", function(object, value) {
    object@miRNA_mature_sequence <- RNAString(as.character(value))
    validObject(object)
    return(object) })


#' @title Set miRNA Seed Sequence
#' @description Setter for seed_sequence attribute of `miRNAGene` class.
#' @param object Object of class `miRNAGene`.
#' @param value New value for seed_sequence.
#' @return The modified object.
#' @seealso
#' \link{seed_sequence<-,miRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' mir1 <- miRNAGene(seed_sequence = 'AUUA')
#' seed_sequence(mir1)
#' seed_sequence(mir1) <- 'AAAGCA'
#' seed_sequence(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("seed_sequence<-",
            function(object, value) standardGeneric("seed_sequence<-"))


#' @title Set miRNA Seed Sequence
#' @description Setter for seed_sequence attribute of `miRNAGene` class.
#' @param object Object of class `miRNAGene`.
#' @param value New value for seed_sequence.
#' @return The modified object.
#' @examples
#' mir1 <- miRNAGene(seed_sequence = 'AUUA')
#' seed_sequence(mir1)
#' seed_sequence(mir1) <- 'AAAGCA'
#' seed_sequence(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("seed_sequence<-", "miRNAGene", function(object, value) {
    object@seed_sequence <- RNAString(as.character(value))
    validObject(object)
    return(object) })

#' @title Create an Object of Class miRNAGene
#' @usage
#' miRNAGene(
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
#'     miRNA_mature_sequence = RNAString(),
#'     seed_sequence = RNAString() )
#' @description Constructor function for `miRNAGene` class.
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
#' @param miRNA_mature_sequence Mature miRNA sequence.
#' @param seed_sequence miRNA seed sequence.
#' @return New miRNAGene object.
#' @examples
#' # An empty miRNAGene object
#' mir1 <- miRNAGene()
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

miRNAGene <- function(gene_ensembl_id = "ENSG00000000000",
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
                    miRNA_mature_sequence = RNAString(),
                    seed_sequence = RNAString()) {

    coordinates <- constructor_helper(chrom = chromosome, str = strand,
                                    g_st = gene_start, g_en = gene_end,
                                    ex_st = exons_starts, ex_en=exons_ends,
                                    pre_st = pre_ncRNA_start,
                                    pre_en = pre_ncRNA_end,
                                    coding = FALSE)
    new("miRNAGene",
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
        miRNA_mature_sequence = RNAString(miRNA_mature_sequence),
        seed_sequence = RNAString(seed_sequence)) }

