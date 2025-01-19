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
