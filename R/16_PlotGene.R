utils::globalVariables(c("label","x1","x2","y","y1","y2"))

#' @title Plot Gene Exons Structure
#' @description The function, relying on ggplot2, plots the gene's exons and
#' prepares the plot for further additions (mRNA or pre-nc RNA).
#' @param x An object of class `Gene` (and derived classes).
#' @param transcr_label Name of the transcript to plot. It's used to name the
#' respective y-tick (e.g. "mRNA").
#' @param transcr_id Slot containing the transcript ENSEMBL ID.
#' @return The gene plot (without x-ticks: this function is not meant to be
#' used alone).
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @importFrom ggarchery geom_arrowsegment
#' @importFrom ggplot2 ggplot geom_rect scale_y_continuous scale_x_continuous
#' labs theme_bw geom_text theme aes element_text element_blank arrow unit

plotGene <- function(x, transcr_label, transcr_id) {
    if (isEmpty(x@exons_coordinates))
        stop("To plot the gene, exons coordinates are necessary.")

    # To plot the gene
    chromosome <- unique(as.character(seqnames(x@exons_coordinates)))
    strand <- unique(as.character(strand(x@exons_coordinates)))
    x_label <- paste("Coordinates relative to", chromosome)

    gene_df <- data.frame(x1 = start(x@exons_coordinates),
                        x2 = end(x@exons_coordinates),
                        y1 = 0.1, y2 = 0.3)
    breaks <- c(gene_df$x1, gene_df$x2)

    ensembl_ids <- data.frame(x = round(sum(range(breaks))/2),
                            y = c(0.4, 0.9),
                            label = c(x@gene_ensembl_id,
                                    slot(x, transcr_id)) )
    tss_label <- ifelse(strand == '+', "TSS", "TTS")
    tts_label <- ifelse(strand == '+', "TTS", "TSS")
    tss_tts <- data.frame(x = c(min(gene_df$x1), max(gene_df$x2)),
                        y = c(0.33, 0.33),
                        label = c(tss_label, tts_label) )
    # Monoexonic gene
    if (length(x@exons_coordinates)==1) {
        introns <- NA }
    else {
        introns <- gaps(x@exons_coordinates)[-1]
        if (strand == '+') {
            starts_intr <- start(introns)-1
            ends_intr <- end(introns) + 1 }
        else {
            starts_intr <- end(introns)+1
            ends_intr <- start(introns) - 1 }   }

    gene_plot <- ggplot() +
        geom_rect(data=gene_df, mapping=aes(xmin=x1, xmax=x2,
                                            ymin=y1, ymax=y2),
                alpha=1, fill = 'purple') +
        scale_y_continuous(breaks = c(0.2, 0.7),
                            limits = c(0, 1),
                            labels = c("Gene", transcr_label)) +
        labs(y = '', x = x_label) + theme_bw() +
        geom_text(data=ensembl_ids, aes( x=x, y=y, label=label),
                color="black", size=5) +
        geom_text(data=tss_tts[1,], aes( x=x, y=y, label=label),
                color="black", size=3, hjust = 0) +
        geom_text(data=tss_tts[2,], aes( x=x, y=y, label=label),
                color="black", size=3, hjust = 1) +
        theme(axis.text.x = element_text(angle=90, size = 8, hjust=0),
            panel.grid.minor.x = element_blank(),
            panel.grid.minor.y = element_blank(),
            axis.title.x = element_text(size = 16),
            axis.text.y = element_text(size = 12))

    # Gene has more than 1 exon
    if (!(any(is.na(introns)))) {
        gene_plot <- gene_plot +
            geom_arrowsegment(aes(x = starts_intr,
                                xend = ends_intr, y = 0.2, yend = 0.2),
                            arrow_positions = 0.5,
                            arrows = arrow(length = unit(0.1, "inches"))) }
    return(gene_plot) }


#' @title Plot NonCodingGene Objects
#' @description The function plots the gene's exons and, if provided,
#' the pre-nc RNA.
#' @param x Object of class `NonCodingGene` (and derived classes).
#' @param y Required parameter to extend the generic function 'plot()'.
#' Its value doesn't influence the final plot.
#' @return Gene plot.
#' @examples
#' mir1 <- miRNAGene(
#'     exons_starts = c(15, 40, 90),
#'     exons_ends = c(20, 50, 100),
#'     pre_ncRNA_start = 44, pre_ncRNA_end = 95,
#'     chromosome = 'chr1', strand = '-',
#'     pre_ncRNA_ensembl_id = "ENST09090911190",
#'     gene_ensembl_id = "ENSG00033300033")
#'
#' plot(mir1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("plot", "NonCodingGene", function(x, y=0) {
    gene_plot <- plotGene(x, "pre-nc RNA", "pre_ncRNA_ensembl_id")
    gene_df <- data.frame(x1 = start(x@exons_coordinates),
                        x2 = end(x@exons_coordinates),
                        y1 = 0.1, y2 = 0.3)
    breaks <- c(gene_df$x1, gene_df$x2)
    if (isEmpty(x@pre_ncRNA_coordinates)) {
        gene_plot <- gene_plot +
            scale_x_continuous(breaks = breaks,
                                limits = c(max(range(breaks)[1]-30,0),
                                            range(breaks)[2]+30),
                                expand = c(0,0))
        return(gene_plot) }
    # Additional plot of the transcript
    else {
        if (!isEmpty(x@pre_ncRNA_coordinates)) {
            intersection_ncRNA_exons <- GenomicRanges::intersect(
                x@exons_coordinates, x@pre_ncRNA_coordinates)
            ncRNA_df <- data.frame(x1 = start(intersection_ncRNA_exons),
                                x2 = end(intersection_ncRNA_exons),
                                y1 = 0.6, y2 = 0.8)
            breaks <- c(gene_df$x1, gene_df$x2,
                        ncRNA_df$x1, ncRNA_df$x2) }

        if (length(intersection_ncRNA_exons)==1) {
            ncRNA_introns <- NA }
        else {
            strand <- unique(as.character(strand(x@exons_coordinates)))
            ncRNA_introns <- gaps(intersection_ncRNA_exons)[-1]
            if (strand == '+') {
                nc_starts_intr <- start(ncRNA_introns)-1
                nc_ends_intr <- end(ncRNA_introns) + 1 }
            else {
                nc_starts_intr <- end(ncRNA_introns)+1
                nc_ends_intr <- start(ncRNA_introns) - 1 } }
        final_plot <- gene_plot +
            geom_rect(data=ncRNA_df, mapping=aes(xmin=x1, xmax=x2,
                                                ymin=y1, ymax=y2),
                    alpha=1, fill = 'red') +
            scale_x_continuous(breaks = breaks,
                                limits = c(max(range(breaks)[1]-30,0),
                                        range(breaks)[2]+30), expand = c(0,0))
        # Transcript spans on more than 1 exon
        if (!any(is.na(ncRNA_introns))) {
            final_plot <- final_plot +
                geom_arrowsegment(aes(x = nc_starts_intr,
                                    xend = nc_ends_intr, y = 0.7, yend = 0.7),
                                arrow_positions = 0.5,
                                arrows = arrow(length = unit(0.1, "inches"))) }
        return(final_plot) } })


#' @title Plot CodingGene Objects
#' @description The function plots the gene's exons and, if provided,
#' the mRNA and CDS. To plot the CDS, mRNA coordinates are necessary,
#' otherwise only the gene's exons are plotted.
#' @param x Object of class `CodingGene`.
#' @param y Required parameter to extend the generic function 'plot()'.
#' Its value doesn't influence the final plot.
#' @return Gene plot.
#' @examples
#' cg1 <- CodingGene(chromosome = 'chrM', strand = '-',
#'                 exons_starts = c(30,60,90),
#'                 exons_ends = c(40,70,100),
#'                 mRNA_start = 35, mRNA_end = 67,
#'                 cds_start = 37, cds_end = 62,
#'                 gene_ensembl_id = "ENSG12312300000",
#'                 mRNA_ensembl_id = "ENST00011133300")
#' plot(cg1)
#' @author Marco Cominelli, \email{marco.cominelli@mail.polimi.it}
#' @export

setMethod("plot", "CodingGene", function(x, y=0) {
    if (isEmpty(x@exons_coordinates))
        stop("To plot the gene object, exons coordinates are necessary.")

    gene_plot <- plotGene(x, "mRNA", "mRNA_ensembl_id")

    gene_df <- data.frame(x1 = start(x@exons_coordinates),
                        x2 = end(x@exons_coordinates),
                        y1 = 0.1, y2 = 0.3)
    breaks <- c(gene_df$x1, gene_df$x2)
    if (isEmpty(x@mRNA_coordinates)) {
        gene_plot <- gene_plot +
            scale_x_continuous(breaks = breaks,
                                limits = c(max(range(breaks)[1]-30,0),
                                            range(breaks)[2]+30),
                                expand = c(0,0))
        return(gene_plot)
    } else {
        intersection_mRNA_exons <- GenomicRanges::intersect(
            x@exons_coordinates, x@mRNA_coordinates)
        mRNA_df <- data.frame(x1 = start(intersection_mRNA_exons),
                            x2 = end(intersection_mRNA_exons),
                            y1 = 0.6, y2 = 0.8)
        strand <- unique(as.character(strand(x@exons_coordinates)))
        if (length(intersection_mRNA_exons)==1) {
            mRNA_introns <- NA }
        else {
            mRNA_introns <- gaps(intersection_mRNA_exons)[-1]
            if (strand == '+') {
                mRNA_starts_intr <- start(mRNA_introns)-1
                mRNA_ends_intr <- end(mRNA_introns) + 1 }
            else {
                mRNA_starts_intr <- end(mRNA_introns)+1
                mRNA_ends_intr <- start(mRNA_introns) - 1 } }

        if (!any(is.na(mRNA_introns))) {
            final_plot <- gene_plot +
                geom_arrowsegment(aes(x = mRNA_starts_intr,
                                xend = mRNA_ends_intr, y = 0.7, yend = 0.7),
                                arrow_positions = 0.5,
                                arrows = arrow(length = unit(0.1, "inches")))
        } else {final_plot <- gene_plot}

        if (!isEmpty(x@cds_coordinates)) {
            intersection_cds_exons <- GenomicRanges::intersect(
                x@exons_coordinates, x@cds_coordinates-1)
            utr <- setdiff(intersection_mRNA_exons,
                            intersection_cds_exons)
            utr_df <- data.frame(x1 = start(utr),
                                x2 = end(utr),
                                y1 = 0.67,
                                y2 = 0.73)
            intersection_cds_exons <- GenomicRanges::intersect(
                x@exons_coordinates, x@cds_coordinates)
            cds_df <- data.frame(x1 = start(intersection_cds_exons),
                                x2 = end(intersection_cds_exons),
                                y1 = 0.6, y2 = 0.8)
            breaks <- c(gene_df$x1, gene_df$x2, utr_df$x1, utr_df$x2,
                        cds_df$x1, cds_df$x2)
            final_plot <- final_plot +
                geom_rect(data=cds_df, mapping=aes(xmin=x1, xmax=x2,
                                                    ymin=y1, ymax=y2),
                        alpha=1, fill = 'red') +
                geom_rect(data=utr_df, mapping=aes(xmin=x1, xmax=x2,
                                                    ymin=y1, ymax=y2),
                        alpha=1, fill = 'red') +
                scale_x_continuous(breaks = breaks,
                                    limits = c(max(range(breaks)[1]-30,0),
                                            range(breaks)[2]+30),
                                    expand = c(0,0))
        }
        # CDS is not defined, but mRNA is.
        else {
            breaks <- c(gene_df$x1, gene_df$x2, mRNA_df$x1, mRNA_df$x2)

            final_plot <- final_plot +
                geom_rect(data=mRNA_df, mapping=aes(xmin=x1, xmax=x2,
                                                    ymin=y1, ymax=y2),
                        alpha=1, fill = 'red') +
                scale_x_continuous(breaks = breaks,
                                    limits = c(max(range(breaks)[1]-30,0),
                                            range(breaks)[2]+30),
                                    expand = c(0,0))
        }
    }
    return(final_plot)
    })

