#' @title CodingGene Class
#' @description This class represents coding genes and derives from `Gene`
#' class.
#' @slot mRNA_ensembl_id `character`. mRNA ENSEMBL ID.
#' @slot mRNA_coordinates `GRanges`. mRNA genomic coordinates.
#' @slot cds_coordinates `GRanges`. CDS genomic coordinates.
#' @slot protein_ensembl_id `character`. Protein ENSEMBL ID.
#' @slot protein_sequence `AAString`. Protein sequence.
#' @slot protein_description `character`. Protein description.
#' @return An object of class `CodingGene`.
#' @seealso
#' \linkS4class{Gene}
#' @examples
#' # An empty CodingGene object.
#' cg1 <- new("CodingGene")
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setClass("CodingGene",
        contains = "Gene",
        slots = c(
            mRNA_ensembl_id = "character",
            mRNA_coordinates = "GRanges",
            cds_coordinates = "GRanges",
            protein_ensembl_id = "character",
            protein_sequence = "AAString",
            protein_description = "character"
        ),
        prototype = list(
            mRNA_ensembl_id = "ENST00000000000",
            mRNA_coordinates = GRanges(),
            cds_coordinates = GRanges(),
            protein_ensembl_id = "ENSP00000000000",
            protein_sequence = AAString(),
            protein_description = NA_character_
        ))


#' @title Validity Check for CodingGene class
#' @description The function checks that: each slot, except the special ones,
#' is of length 1 (or empty); each ID is a valid ENSEMBL ID; all coordinates
#' are consistent w.r.t. each other and that gene and protein sequences are
#' consistent with the provided coordinates.
#' @param object Object of class `CodingGene`.
#' @return `TRUE` if the object is valid, otherwise the function
#' throws an error via `stop()`.
#' @seealso
#' \link{Gene-validity}
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @name CodingGene-validity

setValidity("CodingGene", function(object){

    ### 1. Check on slots length
    special_slots <- c("protein_sequence")
    check_slots_length(obj = object, class_obj = "CodingGene",
                        mother_class_obj = "Gene",special_slots=special_slots)

    ### 2. Check the validity of mRNA and protein ENSEMBL IDs
    mRNA_id <- object@mRNA_ensembl_id
    check_ensembl_id(id = mRNA_id, prefix = 'ENST', type = 'mRNA transcript')

    protein_id <- object@protein_ensembl_id
    check_ensembl_id(id = protein_id, prefix = 'ENSP', type = 'protein')

    ### 3. Check that gene, exons, mRNA and exons coordinates (if defined)
    ###    share both chromosome and strand.
    check_coordinates_consistency(granges_list = list(
        gene_coordinates = object@gene_coordinates,
        exons_coordinates = object@exons_coordinates,
        mRNA_coordinates = object@mRNA_coordinates,
        cds_coordinates = object@cds_coordinates ))

    ### 4. Check on mRNA coordinates w.r.t gene and exons coordinates
    msg <- "mRNA coordinates can't exceed gene coordinates."
    check_boundaries(grange_in = object@mRNA_coordinates,
                    grange_out = object@gene_coordinates, msg = msg)

    check_in_exons(range_to_check = object@mRNA_coordinates,
                    exons = object@exons_coordinates, type = "mRNA")

    ### 5. Check on CDS coordinates w.r.t mRNA and exons coordinates
    msg <- "CDS coordinates can't exceed mRNA coordinates."
    check_boundaries(grange_in = object@cds_coordinates,
                    grange_out = object@mRNA_coordinates, msg = msg)
    # The additional check on gene is done because mRNA coordinates
    # could be not defined.
    msg <- "CDS coordinates can't exceed gene coordinates."
    check_boundaries(grange_in = object@cds_coordinates,
                    grange_out = object@gene_coordinates, msg = msg)

    check_in_exons(range_to_check = object@cds_coordinates,
                    exons = object@exons_coordinates, type = "CDS")

    ### 6. Check on gene sequence length w.r.t. mRNA and CDS coordinates
    coordinates <- list(object@mRNA_coordinates,
                        object@cds_coordinates)
    # If at least one of the two is not empty.
    if (!isEmpty(coordinates)) {
        mRNA_width <- width(coordinates[[1]])
        cds_width <- width(coordinates[[2]])
        min_gene_length <- max(mRNA_width, cds_width)
        if (length(object@gene_sequence)
            && length(object@gene_sequence) < min_gene_length) {
            err_msg <- paste("Length of gene sequence and mRNA/CDS",
                    "coordinates are not coherent: the sequence length",
                    "can't be less than the span of these",
                    "coordinates ranges.")
            stop(err_msg)  }
    }

    ### 7. Check on protein sequence length w.r.t. all defined coordinates
    ###    and gene sequence.
    protein_length_nt <- length(object@protein_sequence)*3
    if (protein_length_nt) {
        if (length(object@gene_sequence)
            && length(object@gene_sequence) < (protein_length_nt+3)) {
            # +3 because of the stop codon.
            err_msg <- paste("Any protein coding gene should be long",
                            "at least '(encoded_protein_length)*3 + 3'.")
            stop(err_msg) }

        coordinates <- list(gene_coordinates = object@gene_coordinates,
                            exons_coordinates = object@exons_coordinates,
                            mRNA_coordinates = object@mRNA_coordinates,
                            cds_coordinates = object@cds_coordinates)
        empty_check <- lapply(coordinates, isEmpty)
        not_empty_coordinates <- coordinates[!unlist( empty_check)]

        # If at least one coordinate is defined
        if (length(not_empty_coordinates)) {
            final_range <- Reduce(GenomicRanges::intersect,
                                not_empty_coordinates)
            if ( all(c('cds_coordinates', 'exons_coordinates') %in%
                    names(not_empty_coordinates)) ) {
                if (sum(width(final_range)) < protein_length_nt) {
                    err_msg <- paste("The length of the protein can't be",
                                "longer than intersection between CDS",
                                "and exons, divided by 3.")
                    stop(err_msg)} }
            else {
                if (!isEmpty(coordinates[[4]])) {
                    protein_limit <- protein_length_nt
                }
                # In this case I need to take into account also the stop codon
                else {protein_limit <- protein_length_nt + 3}

                if (sum(width(final_range)) < protein_limit ) {
                    err_msg <- paste("The protein is too long w.r.t.",
                                "the inserted gene/exons/mRNA/CDS",
                                "coordinates.")
                    stop(err_msg)}  }
        }
    }

    # The check on protein sequence length could be done in a more precise way
    # but in that case, alternative splicing sites and cassette exons should be
    # considered. For example, if a CDS includes 3 exons, this implementation
    # assumes that the internal exon will be entirely included, but this is not
    # always the case (because of alternative splicing sites or cassette exons)
    # so the constraint on protein sequence length is less strict.

    return(TRUE) })


# GETTERs for NonCodingGene

#' @title Get mRNA ENSEMBL ID
#' @description Getter for mRNA_ensembl_id attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @return Value of mRNA_ensembl_id.
#' @seealso
#' \link{mRNA_id,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(mRNA_ensembl_id = 'ENST12341234101')
#' mRNA_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("mRNA_id",
            function(object) standardGeneric("mRNA_id"))


#' @title Get mRNA ENSEMBL ID
#' @description Getter for mRN A_ensembl_id attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @return Value of mRNA_ensembl_id.
#' @examples
#' cg1 <- CodingGene(mRNA_ensembl_id = 'ENST12341234101')
#' mRNA_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mRNA_id", "CodingGene", function(object) {
    object@mRNA_ensembl_id })


#' @title Get mRNA Genomic Coordinates
#' @description Getter for mRNA_coordinates attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @return `GRanges` object representing mRNA coordinates.
#' @seealso
#' \link{mRNA_coordinates,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', mRNA_start = 20,
#' mRNA_end = 30)
#' mRNA_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("mRNA_coordinates",
            function(object) standardGeneric("mRNA_coordinates"))

#' @title Get mRNA Genomic Coordinates
#' @description Getter for mRNA_coordinates attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @return `GRanges` object representing mRNA coordinates.
#' @examples
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', mRNA_start = 20,
#' mRNA_end = 30)
#' mRNA_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mRNA_coordinates", "CodingGene", function(object) {
    object@mRNA_coordinates })


#' @title Get CDS Genomic Coordinates
#' @description Getter for cds_coordinates attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @return `GRanges` object representing CDS coordinates.
#' @seealso
#' \link{cds_coordinates,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', cds_start = 50,
#' cds_end = 90)
#' cds_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("cds_coordinates",
            function(object) standardGeneric("cds_coordinates"))

#' @title Get CDS Genomic Coordinates
#' @description Getter for cds_coordinates attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @return `GRanges` object representing CDS coordinates.
#' @examples
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', cds_start = 50,
#' cds_end = 90)
#' cds_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("cds_coordinates", "CodingGene", function(object) {
    object@cds_coordinates })


#' @title Get Protein ENSEMBL ID
#' @description Getter for protein_ensembl_id attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @return Protein ENSEMBL ID.
#' @seealso
#' \link{protein_id,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(protein_ensembl_id = 'ENSP12312312300')
#' protein_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("protein_id",
            function(object) standardGeneric("protein_id"))


#' @title Get Protein ENSEMBL ID
#' @description Getter for protein_ensembl_id attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @return Protein ENSEMBL ID.
#' @examples
#' cg1 <- CodingGene(protein_ensembl_id = 'ENSP12312312300')
#' protein_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("protein_id", "CodingGene", function(object) {
    object@protein_ensembl_id })


#' @title Get Protein Sequence
#' @description Getter for protein_sequence attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @return Protein sequence.
#' @seealso
#' \link{protein_sequence,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(protein_sequence = 'KRT')
#' protein_sequence(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("protein_sequence",
            function(object) standardGeneric("protein_sequence"))


#' @title Get Protein Sequence
#' @description Getter for protein_sequence attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @return Protein sequence.
#' @examples
#' cg1 <- CodingGene(protein_sequence = 'KRT')
#' protein_sequence(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("protein_sequence", "CodingGene", function(object) {
    object@protein_sequence })


#' @title Get Protein Description
#' @description Getter for protein_description attribute of `CodingGene`
#' class.
#' @param object Object of class `CodingGene`.
#' @return Protein description
#' @seealso
#' \link{protein_description,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(protein_description = "Transcription factor")
#' protein_description(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("protein_description",
            function(object) standardGeneric("protein_description"))


#' @title Get Protein Description
#' @description Getter for protein_description attribute of `CodingGene`
#' class.
#' @param object Object of class `CodingGene`.
#' @return Protein description.
#' @examples
#' cg1 <- CodingGene(protein_description = "Transcription factor")
#' protein_description(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("protein_description", "CodingGene", function(object) {
    object@protein_description })

# SETTERs for NonCodingGene

#' @title Set mRNA ENSEMBL ID
#' @description Setter for mRNA_ensembl_id attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @param value New value for mRNA_ensembl_id.
#' @return The modified object.
#' @seealso
#' \link{mRNA_id<-,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(mRNA_ensembl_id = 'ENST12345678944')
#' mRNA_id(cg1)
#' mRNA_id(cg1) <- 'ENST11111000004'
#' mRNA_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("mRNA_id<-",
            function(object, value) standardGeneric("mRNA_id<-"))


#' @title Set mRNA ENSEMBL ID
#' @description Setter for mRNA_ensembl_id attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @param value New value for mRNA_ensembl_id.
#' @return The modified object.
#' @examples
#' cg1 <- CodingGene(mRNA_ensembl_id = 'ENST12345678944')
#' mRNA_id(cg1)
#' mRNA_id(cg1) <- 'ENST11111000004'
#' mRNA_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mRNA_id<-", "CodingGene", function(object, value) {
    object@mRNA_ensembl_id <- value
    validObject(object)
    return(object) })


#' @title Set mRNA Genomic Coordinates
#' @description Setter for mRNA_coordinates attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @param value Named list with chromosome, strand, start and end values
#' to build the new mRNA coordinates.
#' @return The modified object.
#' @seealso
#' \link{mRNA_coordinates<-,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', mRNA_start = 40,
#' mRNA_end = 80)
#' mRNA_coordinates(cg1)
#' mRNA_coordinates(cg1) <- list(chromosome = 'chr2', strand = '+',
#' start = 30, end = 50)
#' mRNA_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("mRNA_coordinates<-",
            function(object, value) standardGeneric("mRNA_coordinates<-"))


#' @title Set mRNA Genomic Coordinates
#' @description Setter for mRNA_coordinates attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @param value Named list with chromosome, strand, start and end values
#' to build the new mRNA coordinates.
#' @return The modified object.
#' @examples
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', mRNA_start = 40,
#' mRNA_end = 80)
#' mRNA_coordinates(cg1)
#' mRNA_coordinates(cg1) <- list(chromosome = 'chr2', strand = '+',
#' start = 30, end = 50)
#' mRNA_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("mRNA_coordinates<-", "CodingGene", function(object, value) {
    list_names <- c("chromosome", "strand", "start", "end")
    msg <- paste0(" e.g. mRNA_coordinates(gene1) <- list(chromosome =",
                    " chr1, strand = '+', start = 10, end = 50).")
    object@mRNA_coordinates <- create_new_coordinates(value = value,
                                                    list_names = list_names,
                                                    msg = msg)
    validObject(object)
    return(object) })


#' @title Set CDS Genomic Coordinates
#' @description Setter for cds_coordinates attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @param value Named list with chromosome, strand, start and end values
#' to build the new CDS coordinates.
#' @return The modified object.
#' @seealso
#' \link{cds_coordinates<-,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', cds_start = 40,
#' cds_end = 80)
#' cds_coordinates(cg1)
#' cds_coordinates(cg1) <- list(chromosome = 'chr2', strand = '+', start=30,
#' end = 50)
#' cds_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("cds_coordinates<-",
            function(object, value) standardGeneric("cds_coordinates<-"))


#' @title Set CDS Genomic Coordinates
#' @description Setter for cds_coordinates attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @param value Named list with chromosome, strand, start and end values
#' to build the new CDS coordinates.
#' @return The modified object.
#' @examples
#' cg1 <- CodingGene(chromosome = 'chr1', strand = '+', cds_start = 40,
#' cds_end = 80)
#' cds_coordinates(cg1)
#' cds_coordinates(cg1) <- list(chromosome = 'chr2', strand = '+', start=30,
#' end = 50)
#' cds_coordinates(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("cds_coordinates<-", "CodingGene", function(object, value) {
    list_names <- c("chromosome", "strand", "start", "end")
    msg <- paste(" e.g. mRNA_coordinates(gene1) <- list(chromosome =",
                " chr1, strand = '+', start = 10, end = 50).")

    object@cds_coordinates <- create_new_coordinates(value = value,
                                                    list_names = list_names,
                                                    msg = msg)
    validObject(object)
    return(object) })


#' @title Set Protein ENSEMBL ID
#' @description Setter for protein_ensembl_id attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @param value New value for protein_ensembl_id.
#' @return The modified object.
#' @seealso
#' \link{protein_id<-,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(protein_ensembl_id = 'ENSP12345678944')
#' protein_id(cg1)
#' protein_id(cg1) <- 'ENSP11111000004'
#' protein_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("protein_id<-",
            function(object, value) standardGeneric("protein_id<-"))


#' @title Set Protein ENSEMBL ID
#' @description Setter for protein_ensembl_id attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @param value New value for protein_ensembl_id.
#' @return The modified object.
#' @examples
#' cg1 <- CodingGene(protein_ensembl_id = 'ENSP12345678944')
#' protein_id(cg1)
#' protein_id(cg1) <- 'ENSP11111000004'
#' protein_id(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("protein_id<-", "CodingGene", function(object, value) {
    object@protein_ensembl_id <- value
    validObject(object)
    return(object) })


#' @title Set Protein Sequence
#' @description Setter for protein_sequence attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @param value New value for protein_sequence.
#' @return The modified object.
#' @seealso
#' \link{protein_sequence<-,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(protein_sequence = 'KRT')
#' protein_sequence(cg1)
#' protein_sequence(cg1) <- 'KKK'
#' protein_sequence(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("protein_sequence<-",
            function(object, value) standardGeneric("protein_sequence<-"))


#' @title Set Protein Sequence
#' @description Setter for protein_sequence attribute of `CodingGene` class.
#' @param object Object of class `CodingGene`.
#' @param value New value for protein_sequence.
#' @return The modified object.
#' @examples
#' cg1 <- CodingGene(protein_sequence = 'KRT')
#' protein_sequence(cg1)
#' protein_sequence(cg1) <- 'KKK'
#' protein_sequence(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("protein_sequence<-", "CodingGene", function(object, value) {
    object@protein_sequence <-AAString(as.character(value))
    validObject(object)
    return(object) })


#' @title Set Protein Description
#' @description Setter for protein_description attribute of `CodingGene`
#' class.
#' @param object Object of class `CodingGene`.
#' @param value New value for protein_description.
#' @return The modified object.
#' @seealso
#' \link{protein_description<-,CodingGene-method}
#' @examples
#' # This is only the generic function, look in 'See Also' section
#' # for the method implementation.
#' cg1 <- CodingGene(protein_description = 'Transcription factor')
#' protein_description(cg1)
#' protein_description(cg1) <- 'Chromatin remodeling factor'
#' protein_description(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setGeneric("protein_description<-",
            function(object, value) standardGeneric("protein_description<-"))

#' @title Set Protein Description
#' @description Setter for protein_description attribute of `CodingGene`
#' class.
#' @param object Object of class `CodingGene`.
#' @param value New value for protein_description.
#' @return The modified object.
#' @examples
#' cg1 <- CodingGene(protein_description = 'Transcription factor')
#' protein_description(cg1)
#' protein_description(cg1) <- 'Chromatin remodeling factor'
#' protein_description(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("protein_description<-", "CodingGene", function(object, value) {
    object@protein_description <- value
    validObject(object)
    return(object) })


#' @title Create an Object of Class CodingGene
#' @usage
#' CodingGene(
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
#'     mRNA_ensembl_id = "ENST00000000000",
#'     mRNA_start = NA,
#'     mRNA_end = NA,
#'     cds_start = NA,
#'     cds_end = NA,
#'     protein_ensembl_id = "ENSP00000000000",
#'     protein_sequence = AAString(),
#'     protein_description = NA_character_ )
#' @description Constructor function for `CodingGene` class.
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
#' @param mRNA_ensembl_id mRNA ENSEMBL ID.
#' @param mRNA_start mRNA start coordinate.
#' @param mRNA_end mRNA end coordinate.
#' @param cds_start CDS start coordinate.
#' @param cds_end CDS end coordinate.
#' @param protein_ensembl_id Protein ENSEMBL ID.
#' @param protein_sequence Protein sequence.
#' @param protein_description Protein description.
#' @return New CodingGene object.
#' @examples
#' # An empty CodingGene object
#' cg1 <- CodingGene()
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

CodingGene <- function( gene_ensembl_id = "ENSG00000000000",
                        hugo_symbol = NA_character_,
                        gene_complete_name = NA_character_,
                        gene_description = NA_character_,
                        strand = NA, chromosome = NA,
                        gene_start = NA, gene_end = NA,
                        gene_sequence = DNAString(),
                        exons_starts = NA, exons_ends = NA,
                        mRNA_ensembl_id = "ENST00000000000",
                        mRNA_start = NA, mRNA_end = NA,
                        cds_start = NA, cds_end = NA,
                        protein_ensembl_id = "ENSP00000000000",
                        protein_sequence = AAString(),
                        protein_description = NA_character_ ) {

    coordinates <- constructor_helper(chrom = chromosome, str = strand,
                                    g_st = gene_start, g_en = gene_end,
                                    ex_st = exons_starts, ex_en=exons_ends,
                                    m_st = mRNA_start, m_en = mRNA_end,
                                    c_st = cds_start, c_en = cds_end,
                                    coding = TRUE)

    new("CodingGene",
        gene_ensembl_id = gene_ensembl_id, hugo_symbol = hugo_symbol,
        gene_complete_name = gene_complete_name,
        gene_description = gene_description,
        gene_coordinates = coordinates[[1]],
        gene_sequence = DNAString(gene_sequence),
        exons_coordinates = coordinates[[2]],
        mRNA_ensembl_id = mRNA_ensembl_id, mRNA_coordinates=coordinates[[3]],
        cds_coordinates = coordinates[[4]],
        protein_ensembl_id = protein_ensembl_id,
        protein_sequence = AAString(protein_sequence),
        protein_description = protein_description) }

