#' @title Check Length of Each Slot
#' @description The function checks that each slot of `obj`, except
#' `special_slots`, has a length of 1 (or 0).
#' @param obj Object to check.
#' @param class_obj Class of `obj`.
#' @param mother_class_obj Class from which `class_obj` derives.
#' @param special_slots Names of the slots NOT to check.
#' @return `TRUE` if the slots have a valid length, otherwise the function
#' throws an error via `stop()`.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @importFrom Hmisc all.is.numeric
#' @importFrom GenomicRanges GRanges start end seqnames strand intersect width
#' setdiff
#' @importFrom Biostrings AAString DNAString RNAString translate
#' reverseComplement
#' @importFrom methods new validObject slotNames slot
#' @importFrom stringr str_replace_all
#' @importFrom IRanges IRanges findOverlaps reduce gaps
#' @importFrom S4Vectors isEmpty
#' @importFrom dplyr filter

check_slots_length <- function(obj, class_obj,
                                mother_class_obj, special_slots = "") {
    if (is.na(mother_class_obj))
    {
        slots <- as.list(slotNames(class_obj))
    }
    else {
        all_slots <- slotNames(class_obj)
        slots <- all_slots[!all_slots %in% slotNames(mother_class_obj)]
        slots <- as.list(slots) }

    slots_to_check <- slots[!slots %in% special_slots]
    slots_content <- lapply(slots_to_check, function(X) {slot(obj, X)})
    length_check <- lapply(slots_content, length) > 1

    if (any(length_check)) {
        wrong_slots <- paste(slots_to_check[length_check], collapse = ', ')
        err_msg <- paste0("Slot(s) <", wrong_slots,
                            "> can contain at most 1 element!")
        stop(err_msg) }

    return(TRUE)
    }


#' @title Check Validity of ENSEMBL ID
#' @description The function checks that the provided `id` is a valid
#' ENSEMBL ID.
#' @param id ID to check.
#' @param prefix Expected prefix of `id`.
#' @param type Part of the message to print in case of error.
#' @param add_msg Additional part of the message to print in case of error.
#' @return `TRUE` if `id` is a valid ENSEMBL ID, otherwise the function
#' throws an error via `stop()`.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}

check_ensembl_id <- function(id, prefix, type, add_msg = "") {
    if (is.na(id)) {
        return(TRUE) }
    if (substr(id, 1, 4) != prefix || nchar(id) != 15) {
        err_msg <- paste0("Invalid ENSEMBL ", type, " ID", add_msg, "!",
                        " It should start with '", prefix, "' followed by",
                        " exactly 11 digits: e.g. '", prefix, "XXXXXXXXXXX.'")
        stop(err_msg) }

    # The check on '.' is necessary because for example
    # 'all.is.numeric(42.42)' returns TRUE.
    if ( grepl('.', substr(id, 5, 15), fixed = TRUE) ||
        !all.is.numeric(substr(id, 5, 15))) {
        err_mgs <- paste0("Invalid ENSEMBL ", type, " ID", add_msg, "!",
                        " After '", prefix, "' only numbers are accepted (",
                        "exactly 11 digits).")
        stop(err_msg)  }
    return(TRUE)
    }


#' @title Check Coordinates Consinstency
#' @usage
#' check_coordinates_consistency(
#'     granges_list,
#'     return_chr = FALSE,
#'     return_str = FALSE )
#' @description The function checks that the provided `GRanges` objects share
#' both the same chromosome (seqnames) and strand and that chromosome is a
#' valid human chromosome.
#' @param granges_list `list` with `GRanges` objects to check.
#' @param return_chr `logical` to return the chromosome or not.
#' @param return_str `logical` to return the strand or not.
#' @return Chromosome if return_chr is set to TRUE and coordinates are
#' consistent, strand if return_str is set to TRUE and coordinates are
#' consistent, TRUE if the coordinates are consistent and the two previous
#' parameters are set to FALSE, otherwise the function throws an error
#' via `stop()`.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}

check_coordinates_consistency <- function(granges_list, return_chr = FALSE,
                                            return_str = FALSE) {
    granges_names <- names(granges_list)
    chromosome <- unique(unlist(lapply(granges_list,
                                        function(grange)
                                        {as.character(seqnames(grange))})))
    strand <- unique(unlist(lapply(granges_list,
                                        function(grange)
                                            {as.character(strand(grange))})))
    if (length(chromosome) > 1) {
        granges_slots_msg <- paste(granges_names, collapse = ', ')
        err_msg <- paste0("The following coordinates must share the same",
                        " CHROMOSOME < ", granges_slots_msg, " >.")
        stop(err_msg) }
    if (length(strand) > 1) {
        granges_slots_msg <- paste(granges_names, collapse = ', ')
        err_msg <- paste0("The following coordinates must share the same",
                        " STRAND < ", granges_slots_msg, " >.")
        stop(err_msg) }
    human_chromosomes <- paste0("chr", c(seq_len(22), "X", "Y", "M"))
    if (length(chromosome) && !chromosome %in% c(human_chromosomes)) {
        err_msg <- paste0( "Please provide a valid human, or",
                            " mitochondrial chromosome. Choose among {",
                            paste0(human_chromosomes, collapse = ", "), "}.")
        stop(err_msg) }
    if (return_chr) {
        return(chromosome) }
    if (return_str) {
        return(strand) }
    return(TRUE)
    }


#' @title Check Exons Coordinates Validity
#' @description The function checks that the exons coordinates don't overlap
#' (exons represent the base pair of the genome covered by at least one
#' transcript of the gene), because if two exons overlap, they are just part
#' of the same exon. If also gene coordinates are defined, the function
#' checks that the gene starts and ends with an exon (which are first and last
#' exons of the gene), because of the definition of gene.
#' @param exons `GRanges` object representing exons coordinates.
#' @param gene `GRanges` object representing gene coordinates.
#' @return `TRUE` if `exons` is valid and consistent with `gene`, otherwise
#' the function throws an error via `stop()`.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}

check_exons_coordinates <- function(exons, gene) {
    exons <- sort(exons)
    if (length(exons)) {
        exons_width <- sum(width(exons))
        overall_width <- sum(width(reduce(exons)))
        # If they are different, it means there's some
        # overlapping of the exons.
        if (exons_width != overall_width) {
            err_msg <- paste("Exons can't overlap!",
                            "Please insert valid exons coordinates.")
            stop(err_msg)  }
        # A gene ALWAYS starts and ends with an exon.
        if (length(gene)) {
            if (start(exons[1]) != start(gene) ||
                end(exons[length(exons)]) != end(gene)) {
                err_msg <- paste("START of the FIRST exon and END of the",
                                "LAST exon must correspond to START and END",
                                "of the gene, respectively.")
                stop(err_msg)  }
        }
    }
    return(TRUE)
    }


#' @title Create a `GRanges` Object
#' @description The function creates a `GRanges` object given in input
#' chromosome, strand, start and end (contained in `value`).
#' @param value `list` with the necessary components to build
#' the `GRanges` object.
#' @param list_names The names that should be found inside `value`.
#' @param msg Part of the message to print in case of error.
#' @return The new `GRanges` object if names(`value`) are consistent with
#' names in `list_names`, otherwise the function throws an error via
#' `stop()`.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}

create_new_coordinates <- function(value, list_names, msg) {
    if (!all(names(value) %in% list_names) || length(value) != 4) {
        err_msg <- paste0("The list element names, passed as arguments,",
                    " are < ", paste(list_names, collapse = ', '), " > and",
                    " they all need to be specified. ", msg)
        stop(err_msg)}
    new_coordinates <- GRanges(value$chromosome,
                            IRanges(value$start, value$end),
                            strand = value$strand )
    return(new_coordinates)
    }


#' @title Check Coordinates Respect Boundaries
#' @description The function checks that `grange_in` coordinates are inside
#' `grange_out` coordinates.
#' @param grange_in `GRanges` object.
#' @param grange_out `GRanges` object.
#' @param msg Part of the message to print in case of error.
#' @return `TRUE` if coordinates are valid, otherwise the function throws an
#' error via `stop()`.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}

check_boundaries <- function(grange_in, grange_out, msg = "") {
    if (length(grange_in) && length(grange_out)) {
        start_in <- start(grange_in)
        end_in <- end(grange_in)
        start_out <- start(grange_out)
        end_out <- end(grange_out)

        if (start_in < start_out || start_in > end_out || end_in < start_in
            || end_in > end_out) {
            err_msg <- paste0("Out of boundaries coordinates! ", msg)
            stop(err_msg) }
    }
    return(TRUE)
    }


#' @title Check Coordinates are in Exons
#' @description The function checks that `range_to_check` start and end are
#' inside `exons` coordinates range(s).
#' @param range_to_check `GRanges` object representing coordinates to check.
#' @param exons `GRanges` object representing exons coordinates.
#' @param type Part of the message to print in case of error.
#' @return `TRUE` if coordinates are valid, otherwise the function throws an
#' error via `stop()`.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}

check_in_exons <- function(range_to_check, exons, type) {
    if (length(range_to_check) && length(exons)) {
        chrom <- as.character(seqnames(range_to_check))
        strand <- as.character(strand(range_to_check))
        start <- start(range_to_check)
        end <- end(range_to_check)
        start_to_check <- GRanges(chrom, IRanges(start), strand)
        end_to_check <- GRanges(chrom, IRanges(end), strand)
        overlap_start <- findOverlaps(start_to_check, exons)
        overlap_end <- findOverlaps(end_to_check, exons)
        if (length(overlap_start) == 0 || length(overlap_end) == 0) {
            err_msg <- paste0("The ", type,
                            " start and end must be inside the gene exons!")
            stop(err_msg) }
    }
    return(TRUE)
    }


#' @title Check Mature Sequence Validity
#' @description The function checks that the mature sequence length is
#' compatible with the defined genomic coordinates and gene sequence of
#' object. Depending on the value of short, the function also checks that
#' the mature sequence is shorter or longer than 200 nt.
#' @param object Object of class `NonCodingGene`.
#' @param mature_sequence `character` string representing the slot in which
#' the mature sequence to check is contained.
#' @param short `logical` to indicate if `object` is of class
#' SmallNonCodingGene (TRUE) or LongNonCodingGene (FALSE).
#' @return `TRUE` if mature sequence is valid, otherwise the function
#' throws an error via `stop()`.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}

check_mature_length <- function(object, mature_sequence, short = TRUE) {
    coordinates <- list(object@gene_coordinates,
                        object@exons_coordinates,
                        object@pre_ncRNA_coordinates)
    names(coordinates) <- c("gene coordinates", "exons coordinates",
                            "precursor ncRNA coordinates")
    check_coordinates <- unlist(lapply(coordinates, isEmpty))
    non_empty_ranges <- coordinates[!check_coordinates]
    final_range <- Reduce(GenomicRanges::intersect, non_empty_ranges)
    msg <- str_replace_all(mature_sequence,'_',' ')
    if (!isEmpty(final_range)) {
        if (!short && sum(width(final_range)) < 200) {
            err_msg <- paste0("Invalid <",
                        paste(names(non_empty_ranges), collapse = ', '),
                        ">. To encode a lncRNA, at least 200 nt",
                        " are required.")
            stop(err_msg) }
        if (length(slot(object, mature_sequence)) > sum(width(final_range))){
            err_msg <- paste("The length of", mature_sequence,
                        "is not compatible with",
                        "gene/exons/pre-ncRNA coordinates.")
            stop(err_msg) } }
    # Check that mature sequence and gene sequence are of compatible
    # lengths.
    if (length(object@gene_sequence)) {
        if (!short && length(object@gene_sequence) < 200) {
            err_msg <- paste0("A gene encoding for a lncRNA can't be",
                                " shorter than 200 nt.")
            stop(err_msg) }
        if (length(object@gene_sequence) < length(slot(object,
                                                        mature_sequence))) {
            err_msg <- paste0("Gene sequence and ", msg,
                            " are incompatible. The latter can't be longer",
                            " than the former.")
            stop(err_msg) } }
    if (short) {
        if (length(slot(object, mature_sequence)) >= 200) {
            err_msg <- paste0("Invalid ", msg, ". It should be",
                        " shorter than 200 nucleotides.")
            stop(err_msg)} }
    else {
        if (length(slot(object, mature_sequence)) &&
            length(slot(object, mature_sequence)) < 200) {
            err_msg <- paste0("Invalid ", msg, ". It should be",
                        " longer than 200 nucleotides.")
            stop(err_msg)} }
    return(TRUE) }


#' @title Class Constructor Helper Function
#' @usage
#' constructor_helper(
#'     chrom,
#'     str,
#'     g_st,
#'     g_en,
#'     ex_st,
#'     ex_en,
#'     m_st = NA,
#'     m_en = NA,
#'     c_st = NA,
#'     c_en = NA,
#'     pre_st = NA,
#'     pre_en = NA,
#'     coding = TRUE )
#' @description The function, depending on the value of `coding`, generates
#' a `list` with the necessary `GRanges` objects.
#' @param chrom Chromosome.
#' @param str Strand.
#' @param g_st Gene start.
#' @param g_en Gene end.
#' @param ex_st Exons starts.
#' @param ex_en Exons ends.
#' @param m_st mRNA start.
#' @param m_en mRNA end.
#' @param c_st CDS start.
#' @param c_en CDS end.
#' @param pre_st pre-ncRNA start.
#' @param pre_en pre-ncRNA end.
#' @param coding `logical` to indicate if the function will be used in the
#' constructor function of an object of class `CodingGene` or
#' `NonCodingGene`.
#' @return list containing the new GRanges objects if ex_st and ex_en
#' are of the same length, otherwise the function throws an error
#' via `stop()`.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}

constructor_helper <- function(chrom, str, g_st, g_en, ex_st, ex_en,
                                m_st = NA, m_en = NA, c_st = NA, c_en = NA,
                                pre_st = NA, pre_en = NA, coding = TRUE) {
    if (length(ex_st) != length(ex_en)) {
        stop("Each exon START must be matched by an exon END.") }

    if (!any(is.na(c(chrom, str)))) {
        if (!any(is.na(c(g_st, g_en)))) {
            gene_coordinates <- GRanges(chrom, IRanges(g_st,g_en), str) }
        else {gene_coordinates <- GRanges()}

        if (!any(is.na(c(ex_st, ex_en)))) {
            exons_coordinates <- GRanges(chrom, IRanges(ex_st,ex_en), str) }
        else {exons_coordinates <- GRanges()}

        if (coding) {
            if (!any(is.na(c(m_st, m_en)))) {
                mRNA_coordinates <- GRanges(chrom, IRanges(m_st,m_en), str) }
            else {mRNA_coordinates <- GRanges()}
            if (!any(is.na(c(c_st, c_en)))) {
                cds_coordinates <- GRanges(chrom, IRanges(c_st,c_en), str) }
            else {cds_coordinates <- GRanges()}
            output <- list(gene_coordinates, exons_coordinates,
                            mRNA_coordinates, cds_coordinates)
        }
        else {
            if (!any(is.na(c(pre_st, pre_en)))) {
                pre_ncRNA_coordinates <- GRanges(chrom,
                                                IRanges(pre_st,pre_en), str)}
            else {pre_ncRNA_coordinates <- GRanges()}
            output <- list(gene_coordinates, exons_coordinates,
                            pre_ncRNA_coordinates)
        }
    }
    else {
        if (coding) {
            output <- list(GRanges(),GRanges(),GRanges(),GRanges())
        }
        else {
            output <- list(GRanges(),GRanges(),GRanges())
        }
    }
    return(output)
    }


#' @title Compute Product Length
#' @description The function computes the length of attribute `slot` in
#' `object`.
#' @param object An object of class `NonCodingGene`.
#' @param slot Slot of which the length is computed.
#' @param name Name assigned to the output.
#' @return Length of attribute `slot` in `object`.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}

nonCodingLength <- function(object, slot, name) {
    if (length(slot(object, slot))) {
        output <- length(slot(object, slot)) }
    else {
        output <- NA }
    names(output) <- name
    return(output)
    }


#' @title Modify Input List
#' @description The function appends/removes the element specified by `value`
#' from the `list` attribute of `object`.
#' @param object An object of class `RegulatoryNonCodingGene`.
#' @param slot Slot in which there's the list to modify.
#' @param action The action to perform on the list of object in slot.
#' @param value The element to append/remove to/from the list of the object
#' in slot.
#' @param msg Part of the message to print in case of error.
#' @return The modified `list`.
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}

modify_list <- function(object, slot, action, value, msg) {
    if ( !action %in% c("add", "remove") || length(action) != 1) {
        err_msg <- paste0("Action argument accepts only 'add' and 'remove'",
                        " as input values, ", msg)
        stop(err_msg) }
    lis <- slot(object, slot)
    if (action == 'add') {
        lis <- append(lis, value) }
    else {
        lis <- lis[lis != value] }
    lis <- unique(lis)
    return(lis)
    }
