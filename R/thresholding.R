#' Multiple Inflection Points via Second Difference
#'
#' **Find multiple inflection points** by examining the second difference of
#' a numeric series. We look for _local maxima_ in the absolute second difference
#' (optionally above a threshold) as candidate breakpoints.
#'
#' @param x Numeric vector of x-values (ideally sorted).
#' @param y Numeric vector of y-values of the same length as `x`.
#' @param use.threshold Logical indicating whether to only consider second-difference
#'   peaks above `threshold`.
#' @param threshold Numeric value for the minimum magnitude of the second difference.
#'
#' @return A list with:
#' * `inflection_indices`: The indices where \eqn{| \Delta^2 y |} is a local max.
#' * `inflection_xs`: The corresponding `x` values.
#' * `inflection_ys`: The corresponding `y` values.
#' * `second_differences`: The vector of second differences (length = `length(x) - 2`).
#'
#' @details
#' 1. Compute first differences, \eqn{\Delta y[i] = y[i+1] - y[i]}.
#' 2. Compute second differences, \eqn{\Delta^2 y[i] = \Delta y[i+1] - \Delta y[i]}.
#' 3. Identify local maxima in \eqn{| \Delta^2 y |}, optionally only those above `threshold`.
#'
#' @export
#' @author Jared Andrews
#'
#' @examples
#' x_data <- seq(500)
#' # Exponential data
#' y_data <- 1 * exp(0.02 * (x_data - 1))
#'
#' # Sigmoidal
#' y_data2 <- 1 / (1 + exp(-0.02 * (x_data - 250)))
#'
#' res_sd <- get_second_diff(x_data, y_data)
#' res_sd$inflection_xs
#' res_sd$inflection_ys
#'
#' res_sd <- get_second_diff(x_data, y_data2)
#' res_sd$inflection_xs
#' res_sd$inflection_ys
get_second_diff <- function(x, y,
                            use.threshold = TRUE,
                            threshold = 0) {
    n <- length(x)
    if (n < 3) {
        stop("Need at least 3 points for second-difference approach.")
    }
    dy <- diff(y) # first differences
    ddy <- diff(dy) # second differences

    candidate_idxs <- c()

    for (i in seq_len(length(ddy))) {
        left <- if (i == 1) 0 else abs(ddy[i - 1])
        mid <- abs(ddy[i])
        right <- if (i == length(ddy)) 0 else abs(ddy[i + 1])

        is_local_max <- (mid > left) && (mid >= right)
        meets_threshold <- (!use.threshold) || (mid > threshold)

        if (is_local_max && meets_threshold) {
            # i in [1, n-2], so the corresponding x,y index is i+1
            candidate_idxs <- c(candidate_idxs, i + 1)
        }
    }

    list(
        inflection_indices = candidate_idxs,
        inflection_xs      = x[candidate_idxs],
        inflection_ys      = y[candidate_idxs],
        second_differences = ddy
    )
}


#' Multiple Breakpoints via Piecewise Linear Regression
#'
#' **Fit multiple linear segments** using the \pkg{segmented} package. The user
#' specifies the number of breakpoints. The model estimates their _x_-coordinates,
#' and we compute corresponding _y_ predictions at those points.
#'
#' Additionally, we find the **closest actual data point** (in the original
#' `x,y` vectors) to each estimated breakpoint and return those, too.
#'
#' @param x Numeric vector of x-values.
#' @param y Numeric vector of y-values of the same length as `x`.
#' @param n.breakpoints Integer number of breakpoints to identify (>=1).
#'
#' @return A list with:
#' * `model_break_xs`: The **model-estimated** x-coordinates of each breakpoint.
#' * `model_break_ys`: The predicted y-values at those breakpoints (via the segmented model).
#' * `closest_xs`: The x-values from the *original data* nearest each breakpoint.
#' * `closest_ys`: The corresponding y-values from the original data.
#' * `segmented_model`: The resulting `segmented` model object.
#'
#' @details
#' 1. Fit `lm(y ~ x)`.
#' 2. Pass this to `segmented::segmented()`, requesting `n.breakpoints`.
#' 3. Extract the x-locations from `segfit$psi`, then predict y-values at those points.
#' 4. For each estimated breakpoint, find the **nearest actual x** in the data
#'    (via `which.min(abs(x - bps_x[i]))`), and return the corresponding
#'    `(x, y)` as `closest_xs` and `closest_ys`.
#'
#' Note that these breakpoints may not coincide with _actual_ data points, but
#' rather are estimated for an optimal piecewise fit. The additional `closest_xs`
#' and `closest_ys` let you see which real data points are nearest those break estimates.
#'
#' @importFrom stats lm predict
#'
#' @export
#' @author Jared Andrews
#'
#' @examples
#' x_data <- seq(500)
#' # Exponential data
#' y_data <- 1 * exp(0.02 * (x_data - 1))
#'
#' # Sigmoidal
#' y_data2 <- 1 / (1 + exp(-0.02 * (x_data - 250)))
#'
#' res_seg <- get_segmented(x_data, y_data, n.breakpoints = 1)
#' res_seg$model_break_xs # Model-estimated x
#' res_seg$model_break_ys # Predicted y from model at that x
#' res_seg$closest_xs # Actual data x nearest the model breakpoint
#' res_seg$closest_ys # Actual data y for that nearest x
#'
#' res_seg2 <- get_segmented(x_data, y_data2, n.breakpoints = 2)
#' res_seg2$model_break_xs
#' res_seg2$closest_xs
get_segmented <- function(x, y, n.breakpoints = 2) {
    if (n.breakpoints < 1) {
        stop("Number of breakpoints must be >= 1.")
    }

    df <- data.frame(x = x, y = y)
    lmfit <- lm(y ~ x, data = df)

    segfit <- segmented::segmented(lmfit, seg.Z = ~x, npsi = n.breakpoints)

    # Extract model-estimated breakpoints in x and predicted y at those x
    model_bps_x <- segfit$psi[, 2] # second column are the estimates
    model_bps_y <- predict(segfit, newdata = data.frame(x = model_bps_x))

    # For each model breakpoint, find the closest x in the actual data and get its corresponding y
    closest_xs <- numeric(length(model_bps_x))
    closest_ys <- numeric(length(model_bps_x))

    for (i in seq_along(model_bps_x)) {
        idx_closest <- which.min(abs(x - model_bps_x[i]))
        closest_xs[i] <- x[idx_closest]
        closest_ys[i] <- y[idx_closest]
    }

    list(
        model_break_xs   = model_bps_x,
        model_break_ys   = model_bps_y,
        closest_xs       = closest_xs,
        closest_ys       = closest_ys,
        segmented_model  = segfit
    )
}


#' Multiple Elbows via Chord Distance (Local Maxima)
#'
#' **Draw a chord** from the first to the last point, compute perpendicular distances
#' for each \eqn{(x[i], y[i])}, and find _local maxima_ in these distances, optionally
#' above a threshold, as candidate elbows.
#'
#' @param x Numeric vector of x-values (sorted or strictly increasing).
#' @param y Numeric vector of y-values of the same length as `x`.
#' @param use.threshold Logical; if `TRUE`, only local maxima with distance
#'   `>= threshold` are kept.
#' @param threshold Numeric distance threshold.
#'
#' @return A list with:
#' * `inflection_indices`: Indices of local maxima in chord-distance.
#' * `inflection_xs`: The x-values at those indices.
#' * `inflection_ys`: The y-values at those indices.
#' * `distances`: The full array of perpendicular distances to the chord.
#'
#' @details
#' 1. The chord is formed between `(x[1], y[1])` and `(x[n], y[n])`.
#' 2. For each point, the perpendicular distance is computed.
#' 3. Each local maximum is returned as a potential elbow.
#'
#' @export
#' @author Jared Andrews
#'
#' @examples
#' x_data <- seq(500)
#' # Exponential data
#' y_data <- 1 * exp(0.02 * (x_data - 1))
#'
#' # Sigmoidal
#' y_data2 <- 1 / (1 + exp(-0.02 * (x_data - 250)))
#'
#' res_chord <- get_chord_distance(x_data, y_data)
#' res_chord$inflection_xs
#' res_chord$inflection_ys
#'
#' res_chord <- get_chord_distance(x_data, y_data2)
#' res_chord$inflection_xs
#' res_chord$inflection_ys
get_chord_distance <- function(x, y,
                               use.threshold = TRUE,
                               threshold = 0) {
    n <- length(x)
    if (n < 2) {
        stop("Need at least 2 points for chord method.")
    }

    x1 <- x[1]
    y1 <- y[1]
    x2 <- x[n]
    y2 <- y[n]

    dx <- x2 - x1
    dy <- y2 - y1
    chord_len_sq <- dx^2 + dy^2

    if (chord_len_sq == 0) {
        # Degenerate: first == last
        return(list(
            inflection_indices = 1,
            inflection_xs      = x[1],
            inflection_ys      = y[1],
            distances          = rep(0, n)
        ))
    }

    distances <- numeric(n)
    for (i in seq_len(n)) {
        ax <- x[i] - x1
        ay <- y[i] - y1
        cross_val <- abs(ax * dy - ay * dx)
        dist_i <- cross_val / sqrt(chord_len_sq)
        distances[i] <- dist_i
    }

    candidate_idxs <- c()
    for (i in seq(2, n - 1)) {
        left <- distances[i - 1]
        mid <- distances[i]
        right <- distances[i + 1]

        is_local_max <- (mid >= left) && (mid >= right)
        meets_thresh <- (!use.threshold) || (mid >= threshold)

        if (is_local_max && meets_thresh) {
            candidate_idxs <- c(candidate_idxs, i)
        }
    }

    list(
        inflection_indices = candidate_idxs,
        inflection_xs      = x[candidate_idxs],
        inflection_ys      = y[candidate_idxs],
        distances          = distances
    )
}


#' Find Points Closest to n*MAD from the Median
#'
#' This function identifies the point(s) in the data whose `y`-value is
#' *closest* to \eqn{\mathrm{median}(y) \pm n \times \mathrm{MAD}(y)}.
#' In other words, it picks an actual data point near the "(n) MAD boundary"
#' (upper, lower, or both).
#'
#' @param x A numeric vector of x-values (same length as `y`).
#' @param y A numeric vector of y-values (same length as `x`).
#' @param n Numeric value for the number of MADs to use.
#'   3 by default.
#' @param direction One of "upper", "lower", or "both":
#'   - **"upper"**: Find the data point nearest \eqn{\mathrm{median}(y) + n \times \mathrm{MAD}(y)}.
#'   - **"lower"**: Find the data point nearest \eqn{\mathrm{median}(y) - n \times \mathrm{MAD}(y)}.
#'   - **"both"**: Return two points: one near the upper boundary and one near the lower boundary.
#'
#' @return A list with:
#'   - `median_y`: The median of `y`.
#'   - `mad_y`: The median absolute deviation of `y`.
#'   - `upper_threshold`: \eqn{\mathrm{median}(y) + n \times \mathrm{MAD}(y)}.
#'   - `lower_threshold`: \eqn{\mathrm{median}(y) - n \times \mathrm{MAD}(y)}.
#'   - `break_xs`: The x-values closest to the chosen threshold(s).
#'   - `break_ys`: The corresponding y-values.
#'   - `indices`: The indices in `x,y` for these breakpoints.
#'
#' @details
#' 1. We compute:
#'    \deqn{
#'      m = \mathrm{median}(y), \quad
#'      a = \mathrm{MAD}(y)
#'    }
#' 2. Define:
#'    \deqn{
#'      \mathrm{upper} = m + (n)a, \quad
#'      \mathrm{lower} = m - (n)a
#'    }
#' 3. For each threshold requested (*upper*, *lower*, or *both*),
#'    we find the index \eqn{i} that minimizes
#'    \eqn{| y[i] - \mathrm{threshold} |}.
#' 4. We return those `x[i]`, `y[i]` as the "breakpoints."
#'
#' This is a simple “rule-of-thumb” method: once the data cross
#' roughly \eqn{\pm n \times \mathrm{MAD}} from the center, we call that a
#' "significant departure." It doesn't use ranking or second derivatives—just
#' a robust measure of spread to define the boundary.
#'
#' @importFrom stats mad median
#' @export
#' @author Jared Andrews
#'
#' @seealso \code{\link[stats]{mad}}
#'
#' @examples
#' x_data <- seq(500)
#' # Exponential data
#' y_data <- 1 * exp(0.02 * (x_data - 1))
#'
#' # Sigmoidal
#' y_data2 <- 1 / (1 + exp(-0.02 * (x_data - 250)))
#'
#' res_mad <- get_mad(x_data, y_data, n = 3, direction = "upper")
#' res_mad$break_xs
#' res_mad$break_ys
#'
#' res_mad <- get_mad(x_data, y_data2, n = 3, direction = "upper")
#' res_mad$break_xs
#' res_mad$break_ys
get_mad <- function(x, y,
                    n = 3,
                    direction = c("upper", "lower", "both")
                    ) {
    direction <- match.arg(direction)

    if (length(x) != length(y)) {
        stop("x and y must have the same length.")
    }

    m <- median(y)
    a <- mad(y, center = m)

    up_thresh <- m + (n * a)
    low_thresh <- m - (n * a)

    # Function to find index of the data point closest to a threshold th
    find_closest_idx <- function(th) {
        deltas <- abs(y - th)
        which.min(deltas)
    }

    breakpoint_indices <- integer(0)
    if (direction %in% c("upper", "both")) {
        i_up <- find_closest_idx(up_thresh)
        breakpoint_indices <- c(breakpoint_indices, i_up)
    }
    if (direction %in% c("lower", "both")) {
        i_low <- find_closest_idx(low_thresh)
        breakpoint_indices <- c(breakpoint_indices, i_low)
    }

    list(
        median_y        = m,
        mad_y           = a,
        upper_threshold = up_thresh,
        lower_threshold = low_thresh,
        break_xs        = x[breakpoint_indices],
        break_ys        = y[breakpoint_indices],
        indices         = breakpoint_indices
    )
}
