#' @title tRNAGene Class
#' @description Class to represent genes encoding for transfer RNAs.
#' @slot tRNA_mature_sequence `RNAString`. tRNA mature sequence.
#' @slot anticodon `RNAString`. tRNA anticodon.
#' @slot amino_acid `AAString`. tRNA aminoacid.
#' @return An object of class `tRNAGene`.
#' @seealso
#' \linkS4class{HousekeepingNonCodingGene}
#' @examples
#' # An empty tRNAGene object
#' t1 <- new("tRNAGene")
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("tRNAGene",
        contains = "HousekeepingNonCodingGene",
        slots = c(
            tRNA_mature_sequence = "RNAString",
            anticodon = "RNAString",
            amino_acid = "AAString"
        ),
        prototype = list(
            tRNA_mature_sequence = RNAString(),
            anticodon = RNAString(),
            amino_acid = AAString()
        ) )


#' @title Validity Check for tRNAGene Class
#' @description The function checks that: each slot, except the special ones,
#' is of length 1 (or empty); anticodon and aminoacid are consistent one with
#' each other and that maure tRNA sequence is compatible with both the
#' defined coordinates and the gene sequence.
#' @param object Object of class `tRNAGene`.
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \link{HousekeepingNonCodingGene-validity}
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name tRNAGene-validity

setValidity("tRNAGene", function(object) {
    ### 1. Check on slots length
    special_slots <- c("tRNA_mature_sequence", "anticodon")
    check_slots_length(object, "tRNAGene", "HousekeepingNonCodingGene",
                        special_slots = special_slots)
    ### 2. Check anticodon length.
    if (!length(object@anticodon) %in% c(0,3)) {
        stop("Anticodon, if inserted, must be of length 3.") }
    ### 3. Check that anticodon is inside the mature tRNA sequence.
    if (length(object@anticodon) && length(object@tRNA_mature_sequence)) {
        anticodon_in_tRNA <- Biostrings::matchPattern(object@anticodon,
                                            object@tRNA_mature_sequence)
        if (!length(anticodon_in_tRNA)) {
            err_msg <- paste0("Anticodon and mature tRNA sequence",
                        " are not consistent. The former should be found in",
                        " the latter.")
            stop(err_msg) } }
    ### 4. Check consistency between anticodon and aminoacid.
    ###    There are no tRNAs with an anticodon which corresponds to a stop
    ###    codon, so this might be a typo of the user.
    if (length(object@anticodon)) {
        # Considering both codon and anticodon from 5' to 3', then the codon
        # is the reverse complement of the anticodon.
        codon <- reverseComplement(object@anticodon)
        if (toString(codon) %in% c("UAG","UGA","UAA")) {
            err_msg <- paste0("The inserted anticodon corresponds to",
                        " a stop codon! Please enter a valid tRNA anticodon.")
            stop(err_msg)}
        if (length(object@amino_acid)) {
            correct_aminoacid <- Biostrings::translate(codon)
            if (object@amino_acid != correct_aminoacid) {
                err_msg <- paste0("Anticodon and aminoacid are not",
                                    " consistent one with each other.")
                stop(err_msg) } } }
    ### 5. Check length of mature tRNA sequence w.r.t gene/pre-ncRNA lengths.
    check_mature_length(object, 'tRNA_mature_sequence')
    ### 6. Check tRNA mature sequence length if no coordinates are
    ###    defined but gene sequence is.
    if (length(object@gene_sequence)) {
        if (length(object@gene_sequence)<length(object@tRNA_mature_sequence)){
            err_msg <- paste0("The sequence of the gene can't be",
                            " shorter than its mature tRNA transcript.")
            stop(err_msg) } }
    if (length(object@tRNA_mature_sequence)
        && length(object@tRNA_mature_sequence) < 3) {
        err_msg <- paste0("tRNA mature sequence should at least contain",
                        " the anticodon.")
        stop(err_msg) }
    return(TRUE) })


# GETTERs for tRNAGene

#' @title Get Mature Sequence
#' @description Getter for transcript mature sequence.
#' The specific attribute depends on the method implementation.
#' @param object Object of a class derived from `NonCodingGene`.
#' The specific type of object depends on the method implementation.
#' @return Transcript mature sequence.
#' @seealso
#' \link{mature_sequence,tRNAGene-method},
#' \link{mature_sequence,rRNAGene-method},
#' \link{mature_sequence,miRNAGene-method},
#' \link{mature_sequence,siRNAGene-method},
#' \link{mature_sequence,piRNAGene-method},
#' \link{mature_sequence,snRNAGene-method},
#' \link{mature_sequence,snoRNAGene-method},
#' \link{mature_sequence,LongNonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' # Example on a tRNAGene object
#' t1 <- tRNAGene(tRNA_mature_sequence = 'AUGCA')
#' mature_sequence(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("mature_sequence",
            function(object) standardGeneric("mature_sequence"))

#' @title Get Mature tRNA Sequence
#' @description Getter for tRNA_mature_sequence attribute of `tRNAGene` class.
#' @param object Object of class `tRNAGene`.
#' @return Mature tRNA sequence.
#' @examples
#' t1 <- tRNAGene(tRNA_mature_sequence = 'AUGCA')
#' mature_sequence(t1)
#' @export

setMethod("mature_sequence", "tRNAGene", function(object) {
    object@tRNA_mature_sequence })


#' @title Get Anticodon Sequence
#' @description Getter for anticodon attribute of `tRNAGene` class.
#' @param object Object of class `tRNAGene`.
#' @return Anticodon sequence.
#' @seealso
#' \link{anticodon,tRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' t1 <- tRNAGene(anticodon = 'GGG')
#' anticodon(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("anticodon", function(object) standardGeneric("anticodon"))


#' @title Get Anticodon Sequence
#' @description Getter for anticodon attribute of `tRNAGene` class.
#' @param object Object of class `tRNAGene`.
#' @return Anticodon sequence.
#' @examples
#' t1 <- tRNAGene(anticodon = 'GGG')
#' anticodon(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("anticodon", "tRNAGene", function(object) {
    object@anticodon })


#' @title Get tRNA Aminoacid
#' @description Getter for amino_acid attribute of `tRNAGene` class.
#' @param object Object of class `tRNAGene`.
#' @return Aminoacid sequence (one letter, or empty).
#' @seealso
#' \link{amino_acid,tRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' t1 <- tRNAGene(amino_acid = 'K')
#' amino_acid(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("amino_acid", function(object) standardGeneric("amino_acid"))

#' @title Get tRNA Aminoacid
#' @description Getter for amino_acid attribute of `tRNAGene` class.
#' @param object Object of class `tRNAGene`.
#' @return Aminoacid sequence (one letter, or empty).
#' @examples
#' t1 <- tRNAGene(amino_acid = 'K')
#' amino_acid(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("amino_acid", "tRNAGene", function(object) {
    object@amino_acid })


# SETTERs for tRNAGene

#' @title Set Mature Sequence
#' @description Setter for transcript mature sequence.
#' The specific attribute depends on the method implementation.
#' @param object Object of a class derived from `NonCodingGene`.
#' The specific type of object depends on the method implementation.
#' @param value New value for mature sequence.
#' @return The modified object.
#' @seealso
#' \link{mature_sequence<-,tRNAGene-method},
#' \link{mature_sequence<-,rRNAGene-method},
#' \link{mature_sequence<-,miRNAGene-method},
#' \link{mature_sequence<-,siRNAGene-method},
#' \link{mature_sequence<-,piRNAGene-method},
#' \link{mature_sequence<-,snRNAGene-method},
#' \link{mature_sequence<-,snoRNAGene-method},
#' \link{mature_sequence<-,LongNonCodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' # Example on a tRNAGene object
#' t1 <- tRNAGene(tRNA_mature_sequence = 'AUUA')
#' mature_sequence(t1)
#' mature_sequence(t1) <- 'AAAGCA'
#' mature_sequence(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("mature_sequence<-",
            function(object, value) standardGeneric("mature_sequence<-"))

#' @title Set Mature tRNA Sequence
#' @description Setter for tRNA_mature_sequence attribute of `tRNAGene` class.
#' @param object Object of class `tRNAGene`.
#' @param value New value for tRNA_mature_sequence.
#' @return The modified object.
#' @examples
#' t1 <- tRNAGene(tRNA_mature_sequence = 'AUUA')
#' mature_sequence(t1)
#' mature_sequence(t1) <- 'AAAGCA'
#' mature_sequence(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mature_sequence<-", "tRNAGene", function(object, value) {
    object@tRNA_mature_sequence <- RNAString(as.character(value))
    validObject(object)
    return(object) })


#' @title Set tRNA Anticodon
#' @description Setter for anticodon attribute of `tRNAGene` class.
#' @param object Object of class `tRNAGene`.
#' @param value New value for anticodon.
#' @return The modified object.
#' @seealso
#' \link{anticodon<-,tRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' t1 <- tRNAGene(anticodon = 'GGG')
#' anticodon(t1)
#' anticodon(t1) <- 'UUU'
#' anticodon(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("anticodon<-",
            function(object, value) standardGeneric("anticodon<-"))


#' @title Set tRNA Anticodon
#' @description Setter for anticodon attribute of `tRNAGene` class.
#' @param object Object of class `tRNAGene`.
#' @param value New value for anticodon.
#' @return The modified object.
#' @examples
#' t1 <- tRNAGene(anticodon = 'GGG')
#' anticodon(t1)
#' anticodon(t1) <- 'UUU'
#' anticodon(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("anticodon<-", "tRNAGene", function(object, value) {
    object@anticodon <- RNAString(as.character(value))
    validObject(object)
    return(object) })


#' @title Set tRNA Aminoacid
#' @description Setter for amino_acid attribute of `tRNAGene` class.
#' @param object Object of class `tRNAGene`.
#' @param value New value for amino_acid.
#' @return The modified object.
#' @seealso
#' \link{amino_acid<-,tRNAGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' t1 <- tRNAGene(amino_acid = 'K')
#' amino_acid(t1)
#' amino_acid(t1) <- 'T'
#' amino_acid(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("amino_acid<-",
            function(object, value) standardGeneric("amino_acid<-"))


#' @title Set tRNA Aminoacid
#' @description Setter for amino_acid attribute of `tRNAGene` class.
#' @param object Object of class `tRNAGene`.
#' @param value New value for amino_acid.
#' @return The modified object.
#' @examples
#' t1 <- tRNAGene(amino_acid = 'K')
#' amino_acid(t1)
#' amino_acid(t1) <- 'T'
#' amino_acid(t1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("amino_acid<-", "tRNAGene", function(object, value) {
    object@amino_acid <- AAString(as.character(value))
    validObject(object)
    return(object) })


#' @title Create an Object of Class tRNAGene
#' @usage
#' tRNAGene(
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
#'     essentiality_score = NA_real_,
#'     ubiquitous_expression = NA,
#'     tRNA_mature_sequence = RNAString(),
#'     anticodon = RNAString(),
#'     amino_acid = AAString() )
#' @description Constructor function for `tRNAGene` class.
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
#' @param essentiality_score Gene essentiality score.
#' @param ubiquitous_expression `logical` for ubiquitous expression.
#' @param tRNA_mature_sequence tRNA mature sequence.
#' @param anticodon tRNA anticodon.
#' @param amino_acid tRNA aminoacid.
#' @return New tRNAGene object.
#' @examples
#' # An empty tRNAGene object
#' t1 <- tRNAGene()
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

tRNAGene <- function(gene_ensembl_id = "ENSG00000000000",
                    hugo_symbol = NA_character_,
                    gene_complete_name = NA_character_,
                    gene_description = NA_character_,
                    strand = NA, chromosome = NA,
                    gene_start = NA, gene_end = NA,
                    gene_sequence = DNAString(),
                    exons_starts = NA, exons_ends = NA,
                    pre_ncRNA_ensembl_id = "ENST00000000000",
                    pre_ncRNA_start = NA, pre_ncRNA_end = NA,
                    essentiality_score = NA_real_,
                    ubiquitous_expression = NA,
                    tRNA_mature_sequence = RNAString(),
                    anticodon = RNAString(),
                    amino_acid = AAString()) {

    coordinates <- constructor_helper(chrom = chromosome, str = strand,
                                    g_st = gene_start, g_en = gene_end,
                                    ex_st = exons_starts, ex_en=exons_ends,
                                    pre_st = pre_ncRNA_start,
                                    pre_en = pre_ncRNA_end,
                                    coding = FALSE)
    new("tRNAGene",
        gene_ensembl_id = gene_ensembl_id,
        hugo_symbol = hugo_symbol,
        gene_complete_name = gene_complete_name,
        gene_description = gene_description,
        gene_coordinates = coordinates[[1]],
        gene_sequence = DNAString(gene_sequence),
        exons_coordinates = coordinates[[2]],
        pre_ncRNA_ensembl_id = pre_ncRNA_ensembl_id,
        pre_ncRNA_coordinates = coordinates[[3]],
        essentiality_score = essentiality_score,
        ubiquitous_expression = ubiquitous_expression,
        tRNA_mature_sequence = RNAString(tRNA_mature_sequence),
        anticodon = RNAString(anticodon),
        amino_acid = AAString(amino_acid)) }
