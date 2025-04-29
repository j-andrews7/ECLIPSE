#' Plot enhancer ranking curve
#'
#' Plots the enhancer ranking curve based on the signal values in the regions,
#' with visual indicators for the signal cutoff and the number of super-enhancers identified.
#'
#' @param regions A GRanges object containing `rankby_signal` and `super` columns, as well as the `threshold` value in its metadata.
#' @param factor.label Optional character string for prefixing the y-axis label with a specific factor,
#'   e.g. "(H3K27ac)".
#'   Default is `NULL`.
#' @param point.color Character string specifying the color of the points.
#'   Default is "grey50".
#' @param point.size Numeric value specifying the size of the points.
#'   Default is 1.
#' @param se.point.color Character string specifying the color of the super enhancer points.
#'   Default is "red".
#' @param se.point.size Numeric value specifying the size of the super enhancer points.
#'   Default is 1.5.
#' @param plot.threshold.x Logical indicating whether to plot the x-axis threshold line.
#'   Default is `TRUE`.
#' @param plot.threshold.y Logical indicating whether to plot the y-axis threshold line.
#'   Default is `TRUE`.
#' @param threshold.color Character string specifying the color of the threshold line(s).
#'   Default is "grey30".
#' @param threshold.size Numeric value specifying the width of the threshold line(s).
#'   Default is 0.75.
#' @param show.legend Logical indicating whether to show the legend.
#'   Default is `TRUE`.
#' @param do.raster Logical indicating whether to render the plot as a raster image.
#'   Default is `FALSE`.
#' @param raster.dpi Numeric value specifying the DPI of the raster image if `do.raster` is `TRUE`.
#' @param return.plotly Logical indicating whether to return a plotly object.
#'   Default is `FALSE`.
#'
#' @return A ggplot or plotly object.
#'
#' @author Jared Andrews, based on code from ROSE.
#'
#' @export
#' @importFrom S4Vectors metadata
#' @import ggplot2
#' @importFrom plotly ggplotly
#'
#' @examples
#' library(GenomicRanges)
#'
#' set.seed(7)
#' num_elements <- 1000
#' regions <- GRanges(
#'     seqnames = rep("chr1", num_elements),
#'     ranges = IRanges(
#'         start = seq(1, by = 1000, length.out = num_elements),
#'         end = seq(1000, by = 1000, length.out = num_elements)
#'     )
#' )
#'
#' # Add random signal values and classify as super enhancers
#' regions$rankby_signal <- c(
#'     rnorm(num_elements - 50, mean = 150, sd = 10),
#'     rnorm(50, mean = 250, sd = 50)
#' )
#' regions$super <- regions$rankby_signal > 170
#' metadata(regions)$threshold <- 170
#'
#' # Plot the enhancer curve
#' plot_enhancer_curve(regions)
plot_enhancer_curve <- function(regions,
                                factor.label = NULL,
                                point.color = "black",
                                point.size = 1,
                                se.point.color = "red",
                                se.point.size = 1.5,
                                plot.threshold.x = TRUE,
                                plot.threshold.y = TRUE,
                                threshold.color = "grey30",
                                threshold.size = 0.75,
                                show.legend = TRUE,
                                do.raster = FALSE,
                                raster.dpi = 300,
                                return.plotly = FALSE) {
    if (do.raster && !requireNamespace("ggrastr", quietly = TRUE)) {
        stop("Install the 'ggrastr' package (install.packages('ggrastr')) to use the 'do.raster' option.")
    }

    ylabel <- ifelse(metadata(regions)$control_subtracted, "Sample Signal - Control Signal", "Sample Signal")

    if (!is.null(factor.label)) {
        ylabel <- paste(factor.label, ylabel)
    }

    rev_rank <- rev(seq_len(length(regions$rankby_signal)))
    ranked_sig <- regions$rankby_signal[order(regions$rankby_signal, decreasing = TRUE)]

    df <- data.frame(
        rev_rank = rev_rank,
        ranked_sig = ranked_sig,
        super = regions$super[order(regions$rankby_signal, decreasing = TRUE)]
    )

    p <- ggplot(df, aes(x = rev_rank, y = ranked_sig)) +
        geom_point(aes(color = super, size = super)) +
        scale_color_manual(values = c("FALSE" = point.color, "TRUE" = se.point.color)) +
        scale_size_manual(values = c("FALSE" = point.size, "TRUE" = se.point.size)) +
        labs(x = "Putative Enhancers", y = ylabel) +
        theme_classic()

    if (plot.threshold.x) {
        p <- p + geom_vline(
            xintercept = length(regions$rankby_signal) - sum(regions$super),
            color = threshold.color, size = threshold.size, linetype = "dashed"
        )
    }

    if (plot.threshold.y) {
        p <- p + geom_hline(
            yintercept = metadata(regions)$threshold, color = threshold.color,
            size = threshold.size, linetype = "dashed"
        )
    }

    if (show.legend) {
        p <- p + guides(color = guide_legend(title = "Super Enhancer"))
    } else {
        p <- p + theme(legend.position = "none")
    }

    if (do.raster) {
        p <- ggrastr::rasterize(p, dpi = raster.dpi)
    }

    p <- p + annotate("text", x = 0.2 * max(rev_rank), y = 0.8 * max(ranked_sig), label = paste(
        " Cutoff used: ", metadata(regions)$threshold, "\n",
        "Super-Enhancers identified: ", sum(regions$super)
    ), hjust = 0)

    if (return.plotly) {
        p <- ggplotly(p)
    }

    p
}


#' Generate QQ Plot for a Specified Distribution
#'
#' @param data A data.frame containing the values to plot.
#'   Alternatively, a GRanges object, in which case the values are extracted from the specified column of `mcols(data)`.
#' @param column A string specifying the column name to use for values.
#' @param dist A string specifying the distribution to fit. Supported distributions include:
#'   "beta", "cauchy", "chi-squared", "exponential", "gamma", "geometric", "log-normal", "lognormal",
#'   "logistic", "negative binomial", "normal", "Poisson", "t", and "weibull".
#' @param point.size Numeric value specifying the size of points in the QQ plot.
#'   Default is 1.
#' @param point.color A string specifying the color of points in the QQ plot.
#'   Default is "blue".
#' @param line.width Numeric value specifying the width of the reference line.
#'   Default is 0.8.
#' @param line.color A string specifying the color of the reference line.
#'   Default is "red".
#' @param do.raster Logical indicating whether to render the plot as a raster image.
#'   Default is `FALSE`.
#' @param raster.dpi Numeric value specifying the DPI of the raster image if `do.raster` is `TRUE`.
#' @param return.plotly Logical indicating whether to return a plotly object.
#'   Default is `FALSE`.
#' @param ... Additional arguments to pass to `MASS::fitdistr`.
#'
#' @return A ggplot2 object.
#'
#' @importFrom ggplot2 geom_point geom_abline labs theme_bw
#' @importFrom MASS fitdistr
#' @importFrom stats ppoints qnorm qexp qgamma qpois qweibull qbeta qcauchy qchisq qgeom qlnorm qlogis qnbinom qt
#' @importFrom S4Vectors mcols
#'
#' @export
#'
#' @author Jared Andrews
#'
#' @examples
#' gr <- GRanges(
#'     seqnames = rep("chr1", 1000),
#'     ranges = IRanges(start = seq(1, by = 1000, length.out = 1000), end = seq(1000, by = 1000, length.out = 1000))
#' )
#' mcols(gr) <- DataFrame(values = rnorm(1000))
#' plot_qq(data = gr, column = "values", dist = "normal")
#' plot_qq(data = gr, column = "values", dist = "gamma")
plot_qq <- function(
    data, column, dist,
    point.size = 1,
    point.color = "blue",
    line.width = 0.8,
    line.color = "black",
    do.raster = FALSE,
    raster.dpi = 300,
    return.plotly = FALSE,
    ...) {
    if (do.raster && !requireNamespace("ggrastr", quietly = TRUE)) {
        stop("Install the 'ggrastr' package (install.packages('ggrastr')) to use the 'do.raster' option.")
    }

    if (return.plotly && !requireNamespace("plotly", quietly = TRUE)) {
        stop("Install the 'plotly' package (install.packages('plotly')) to use the 'return.plotly' option.")
    }

    # Check if the data is a GRanges object
    if (class(data) == "GRanges") {
        if (!column %in% colnames(mcols(data))) {
            stop("The specified column does not exist in the object.")
        }

        values <- mcols(data)[[column]]
    } else if (class(data) == "data.frame" || class(data) == "DFrame") {
        if (!column %in% colnames(data)) {
            stop("The specified column does not exist in the object.")
        }

        values <- data[[column]]
    } else {
        stop("The data must be a GRanges or data.frame object.")
    }

    if (!is.numeric(values)) {
        stop("The specified column must contain numeric values.")
    }

    empirical_quantiles <- sort(values)

    # Round the empirical quantiles for discrete distributions
    if (dist %in% c("poisson", "geometric", "negative binomial")) {
        empirical_quantiles <- round(empirical_quantiles)
    }

    # Fit the specified distribution to the data
    fit <- tryCatch(
        MASS::fitdistr(empirical_quantiles, dist, ...),
        error = function(e) stop("Error in fitting the distribution: ", e$message)
    )

    # Generate theoretical quantiles and empirical quantiles
    n <- length(values)
    theoretical_quantiles <- switch(dist,
        "normal" = qnorm(ppoints(n), mean = fit$estimate["mean"], sd = fit$estimate["sd"]),
        "exponential" = qexp(ppoints(n), rate = fit$estimate["rate"]),
        "gamma" = qgamma(ppoints(n), shape = fit$estimate["shape"], rate = fit$estimate["rate"]),
        "poisson" = qpois(ppoints(n), lambda = fit$estimate["lambda"]),
        "weibull" = qweibull(ppoints(n), shape = fit$estimate["shape"], scale = fit$estimate["scale"]),
        "beta" = qbeta(ppoints(n), shape1 = fit$estimate["shape1"], shape2 = fit$estimate["shape2"]),
        "cauchy" = qcauchy(ppoints(n), location = fit$estimate["location"], scale = fit$estimate["scale"]),
        "chi-squared" = qchisq(ppoints(n), df = fit$estimate["df"]),
        "geometric" = qgeom(ppoints(n), prob = fit$estimate["prob"]),
        "log-normal" = qlnorm(ppoints(n), meanlog = fit$estimate["meanlog"], sdlog = fit$estimate["sdlog"]),
        "lognormal" = qlnorm(ppoints(n), meanlog = fit$estimate["meanlog"], sdlog = fit$estimate["sdlog"]),
        "logistic" = qlogis(ppoints(n), location = fit$estimate["location"], scale = fit$estimate["scale"]),
        "negative binomial" = qnbinom(ppoints(n), size = fit$estimate["size"], mu = fit$estimate["mu"]),
        "t" = qt(ppoints(n), df = fit$estimate["df"]),
        stop("Unsupported distribution: ", dist)
    )

    qq_df <- data.frame(
        Theoretical = theoretical_quantiles,
        Empirical = empirical_quantiles
    )

    p <- ggplot(qq_df, aes(x = Theoretical, y = Empirical)) +
        geom_point(size = point.size, color = point.color) +
        geom_abline(intercept = 0, slope = 1, linetype = "solid", color = line.color, linewidth = line.width) +
        labs(
            title = paste("QQ Plot for", dist, "distribution"),
            x = "Theoretical Quantiles",
            y = "Empirical Quantiles"
        ) +
        theme_bw()

    if (do.raster) {
        p <- ggrastr::rasterize(p, dpi = raster.dpi)
    }

    if (return.plotly) {
        p <- plotly::ggplotly(p)
    }

    p
}


#' Plot Locus
#'
#' This function generates a locus visualization plot using genomic regions,
#' signal tracks, and fold-change data.
#' It supports highlighting significant
#' regions and customizing the appearance of tracks.
#'
#' @details
#' Coverage tracks are generated from BAM files and normalized by library size to reads per million.
#' While not a perfect normalization, as it doesn't account for composition effects, it is often
#' sufficient for visualization of robust differences between groups with sufficient numbers of replicates.
#'
#' @param region.gr [GenomicRanges::GRanges-class] object containing a single region to plot.
#' @param window.gr [GenomicRanges::GRanges-class] object representing merged windows with
#'   associated statistics, as returned by [find_differential()].
#' @param window.rse [SummarizedExperiment::RangedSummarizedExperiment-class] object
#'   containing read count data and library sizes, as returned by [find_differential()].
#' @param g1.bam.files Named list of BAM file paths for group 1.
#' @param g2.bam.files Named list of BAM file paths for group 2.
#' @param genome String specifying the genome assembly.
#'   This is generally expected to match the UCSC conventions for genome assemblies.
#'   Default is "hg38".
#' @param region.track.name String for the name of the region track.
#'   Default is "SEs".
#' @param region.fill String specifying the fill color for the region track.
#'   Default is "orange".
#' @param padding Numeric value specifying the padding around regions.
#'   This is used to expand the region for additional context by expanding
#'   the region by a factor of (1 + padding).
#'   Default is 0.3, e.g. 30% of the region size.
#' @param fdr.thresh Numeric value specifying the FDR threshold for significance.
#'   Default is 0.05.
#' @param fc.track.position String specifying the position of the fold-change track
#'   relative to signal tracks. Options are "between", "above", or "below".
#'   Default is "between".
#' @param g1.bw.files Optional list of BigWig file paths for group 1.
#'   Not implemented currently.
#' @param g2.bw.files Optional list of BigWig file paths for group 2.
#'   Not implemented currently.
#' @param g1.fill String specifying the fill color for group 1 tracks.
#'   Default is "red".
#' @param g2.fill String specifying the fill color for group 2 tracks.
#'   Default is "blue".
#' @param param [csaw::readParam-class] object specifying read extraction parameters.
#'   Default is `readParam(minq = 20, dedup = FALSE)`.
#' @param top.tracks Optional list of additional tracks to display above the main tracks.
#'   Default is NULL.
#' @param bottom.tracks Optional list of additional tracks to display below the main tracks.
#'   Default is NULL.
#' @param main String specifying the main title of the plot.
#'   Default is "".
#' @param region.lfc.color String specifying the color of the log2 fold-change track.
#'   Default is "black".
#' @param highlight.sig Logical indicating whether to highlight significant regions.
#'   Default is TRUE.
#' @param higlight.fill String specifying the fill color for highlighted regions.
#'   Default is "#d3ff8c".
#' @param highlight.color String specifying the border color for highlighted regions.
#'   Default is "#8ff7df".
#' @param ... Additional arguments passed to the `plotTracks` function.
#'
#' @return A list of Gviz tracks.
#'
#' @examples
#' # Example usage:
#' plot_locus(
#'   region.gr = my_regions,
#'   window.gr = my_windows,
#'   window.rse = my_rse,
#'   g1.bam.files = list("sample1" = "path/to/sample1.bam"),
#'   g2.bam.files = list("sample2" = "path/to/sample2.bam"),
#'   genome = "hg38",
#'   region.track.name = "My Regions",
#'   padding = 0.2,
#'   fdr.thresh = 0.01
#' )
#'
#' @importFrom GenomicRanges start end seqnames width resize
#' @importFrom Gviz plotTracks DataTrack AnnotationTrack HighlightTrack
#' @importFrom csaw extractReads readParam
#' 
#' @author Jared Andrews
#' 
#' @export
plot_locus <- function(
    region.gr,
    window.gr,
    window.rse,
    g1.bam.files,
    g2.bam.files,
    genome = "hg38",
    region.track.name = "SEs",
    region.fill = "orange",
    padding = 0.3,
    fdr.thresh = 0.05,
    fc.track.position = c("between", "above", "below"),
    g1.bw.files = NULL,
    g2.bw.files = NULL,
    g1.fill = "red",
    g2.fill = "blue",
    param = readParam(minq = 20, dedup = FALSE),
    top.tracks = NULL,
    bottom.tracks = NULL,
    main = "",
    region.lfc.color = "black",
    highlight.sig = TRUE,
    higlight.fill = "#d3ff8c",
    highlight.color = "#8ff7df",
    ...) {
    fc.track.position <- match.arg(fc.track.position)

    # Get significant regions
    top_se_reg_sig <- window.gr[window.gr$FDR < fdr.thresh]

    # Region track
    padded_gr <- resize(region.gr, width = width(region.gr) * (1 + padding), fix = "center")
    reg_track <- AnnotationTrack(region.gr,
        genome = genome, name = region.track.name,
        background.title = "white", col.axis = "black",
        col.title = "black", fill = region.fill
    )

    # Get significant regions
    top_se_reg_sig <- window.gr[window.gr$FDR < fdr.thresh]

    # Region foldchanges
    fctrack <- DataTrack(
        range = top_se_reg, data = top_se_reg$rep.logFC, name = "Region FCs",
        genome = genome, ylim = c(
            -max(abs(top_se_reg$rep.logFC)),
            max(abs(top_se_reg$rep.logFC))
        ),
        type = "b", col = region.lfc.color,
        background.title = "white", baseline = 0, col.axis = "black",
        col.title = "black"
    )


    # Signal Tracks
    # naive1_track <- DataTrack(range = naiveB1_treat_bw_path, type = "histogram",
    #                           name = "naiveB_r1", genome = "hg38",
    #                         fill = "blue", background.title = "blue", baseline = 0)
    # naive2_track <- DataTrack(range = naiveB2_treat_bw_path, type = "histogram",
    #                           name = "naiveB_r2", genome = "hg38",
    #                         fill = "blue", background.title = "blue", baseline = 0)
    # act1_track <- DataTrack(range = actB1_treat_bw_path, type = "histogram",
    #                         name = "activatedB_r1", genome = "hg38",
    #                         fill = "red", background.title = "red", baseline = 0)
    # act2_track <- DataTrack(range = actB2_treat_bw_path, type = "histogram",
    #                         name = "activatedB_r2", genome = "hg38",
    #                         fill = "red", background.title = "red", baseline = 0)

    # Signal Tracks from BAMs
    g1_collected <- list()
    g1_covs <- list()
    g2_collected <- list()
    g2_covs <- list()
    max_cov <- 1
    g1_lib_sizes <- window.rse$totals[colnames(window.rse) %in% names(g1.bam.files)] / 1e6
    g2_lib_sizes <- window.rse$totals[colnames(window.rse) %in% names(g2.bam.files)] / 1e6

    for (i in seq_along(g1.bam.files)) {
        reads <- extractReads(bam.file = g1.bam.files[[i]], padded_gr, param = param)
        cov <- as(coverage(reads) / g1_lib_sizes[i], "GRanges")
        g1_covs[[i]] <- cov
        if (max(cov$score) > max_cov) {
            max_cov <- max(cov$score)
        }
    }

    for (i in seq_along(g2.bam.files)) {
        reads <- extractReads(bam.file = g2.bam.files[[i]], padded_gr, param = param)
        cov <- as(coverage(reads) / g2_lib_sizes[i], "GRanges")
        g2_covs[[i]] <- cov
        if (max(cov$score) > max_cov) {
            max_cov <- max(cov$score)
        }
    }

    for (i in seq_along(g1_covs)) {
        g1_collected[[i]] <- DataTrack(g1_covs[[i]],
            type = "histogram", lwd = 0, ylim = c(0, max_cov),
            name = names(g1.bam.files)[i], col.axis = "black",
            col.title = "black", genome = genome,
            col.histogram = g1.fill, background.title = "white", baseline = 0
        )
    }

    for (i in seq_along(g2_covs)) {
        g2_collected[[i]] <- DataTrack(g2_covs[[i]],
            type = "histogram", lwd = 0, ylim = c(0, max_cov),
            name = names(g2.bam.files)[i], col.axis = "black",
            col.title = "black", genome = genome,
            col.histogram = g2.fill, background.title = "white", baseline = 0
        )
    }

    tl <- c(reg_track, g1_collected, fctrack, g2_collected)

    if (fc.track.position == "above") {
        tl <- c(reg_track, fctrack, g1_collected, g2_collected)
    } else if (fc.track.position == "below") {
        tl <- c(reg_track, g1_collected, g2_collected, fctrack)
    }

    # Highlight sig bins
    if (highlight.sig) {
        ht <- HighlightTrack(
            trackList = tl,
            start = start(top_se_reg_sig),
            end = end(top_se_reg_sig),
            chromosome = as.character(seqnames(top_se_reg_sig)),
            fill = "#d3ff8c",
            col = "#8ff7df"
        )
    } else {
        ht <- tl
    }

    tracklist <- unlist(list(top.tracks, ht, bottom.tracks))

    browser()

    plotTracks(tracklist,
        from = start(padded_gr),
        to = end(padded_gr),
        chromosome = as.character(seqnames(padded_gr)),
        transcriptAnnotation = "symbol",
        collapseTranscripts = "longest",
        main = main,
        ...
    )
}
