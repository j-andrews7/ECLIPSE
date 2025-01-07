# X11 License
# Copyright (C) 2013 by Whitehead Institute for Biomedical Research
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE WHITEHEAD INSTITUTE FOR BIOMEDICAL RESEARCH BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
# USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# Except as contained in this notice, the name of the Whitehead Institute for Biomedical Research shall not be used in
# advertising or otherwise to promote the sale, use or other dealings in this Software without prior written authorization
# from the Whitehead Institute for Biomedical Research.

### Below functions from original ROSE_callSuper.R ###

#' Calculate cutoff by sliding a diagonal line
#'
#' Calculates a cutoff value by sliding a diagonal line over the sorted input vector and
#' finding the point where it is tangential (or as close as possible) to the data.
#'
#' @param inputVector Numeric vector of data points.
#' @param drawPlot Logical indicating whether to draw a plot. Defaults to `TRUE``.
#' @param ... Additional arguments passed to the `plot`` function.
#'
#' @return A list containing:
#' - `absolute`: The cutoff value (y-value at the tangential point).
#' - `overMedian`: Fold over the median (cutoff divided by the median of `inputVector`).
#' - `overMean`: Fold over the mean (cutoff divided by the mean of `inputVector`).
#'
#' @author Young Lab, with documentation by Jared Andrews
#' 
#' @importFrom stats optimize median
#' @importFrom graphics plot points abline title axis
#' @keywords internal
#'
#' @rdname INTERNAL_calculate_cutoff
calculate_cutoff <- function(inputVector, drawPlot = TRUE, slope = NULL, ...) {
    inputVector <- sort(inputVector)
    inputVector[inputVector < 0] <- 0 # Set negative values to zero
    if (is.null(slope)) {
        slope <- (max(inputVector) - min(inputVector)) / length(inputVector)
    }
    xPt <- floor(optimize(numPts_below_line,
        lower = 1, upper = length(inputVector),
        myVector = inputVector, slope = slope
    )$minimum)
    y_cutoff <- inputVector[xPt] # The cutoff value

    if (drawPlot) {
        plot(1:length(inputVector), inputVector, type = "l", ...)
        b <- y_cutoff - (slope * xPt)
        abline(v = xPt, h = y_cutoff, lty = 2, col = 8)
        points(xPt, y_cutoff, pch = 16, cex = 0.9, col = 2)
        abline(coef = c(b, slope), col = 2)
        title(paste("x=", xPt,
            "\ny=", signif(y_cutoff, 3),
            "\nFold over Median=", signif(y_cutoff / median(inputVector), 3), "x",
            "\nFold over Mean=", signif(y_cutoff / mean(inputVector), 3), "x",
            sep = ""
        ))
        axis(1, sum(inputVector == 0), sum(inputVector == 0), col.axis = "pink", col = "pink")
    }
    return(list(
        absolute = y_cutoff,
        overMedian = y_cutoff / median(inputVector),
        overMean = y_cutoff / mean(inputVector)
    ))
}


#' Determine number of points below a diagonal line
#'
#' Determines the number of points below a diagonal line passing through a given point in the data.
#'
#' @param myVector Numeric vector of sorted data points.
#' @param slope Numeric value of the slope of the diagonal line.
#' @param x Numeric value representing the x-coordinate of the point through which the line passes.
#'
#' @return An integer representing the number of points in `myVector` that lie below the diagonal line.
#'
#' @keywords internal
#' @author Young Lab, with documentation by Jared Andrews
#' @rdname INTERNAL_numPts_below_line
numPts_below_line <- function(myVector, slope, x) {
    yPt <- myVector[x]
    b <- yPt - (slope * x)
    xPts <- 1:length(myVector)
    return(sum(myVector <= (xPts * slope + b)))
}
