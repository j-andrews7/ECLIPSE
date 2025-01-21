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
