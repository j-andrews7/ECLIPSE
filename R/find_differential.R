#' Perform differential analysis of super enhancers between conditions
#'
#' This function identifies differential signal between two conditions
#' using a negative binomial quasi-likelihood model. It processes enhancer regions,
#' performs window-based read counting, normalization, filtering based on
#' background signals, statistical testing, and merging of windows into
#' consolidated regions.
#'
#' @details
#' The analysis proceeds through the following steps:
#' 1. Regions are segmented into sliding windows.
#' 2. Reads from BAM files are counted in each window.
#' 3. Windows with low read counts relative to background are filtered out.
#' 4. Normalization is performed to account for trended bias.
#' 5. Dispersion is modeled using a negative binomial generalized linear
#'    model with quasi-likelihood (QL).
#' 6. Differential analysis is performed using a QL F-test to
#'    assign statistical significance to each window.
#'    Group 2 will be set as the reference level.
#' 7. Windows are merged into regions based on proximity, and test statistics
#'    from individual windows are combined via [csaw::mergeResults()].
#' 8. Differential analysis statistics are added to the metadata of merged regions.
#'
#' @param regions A `GRanges` object of regions to compare between conditions.
#'   These regions will be reduced, i.e. overlapping regions will be merged to one.
#' @param g1.bam.files A character vector of BAM file paths for group 1.
#' @param g2.bam.files A character vector of BAM file paths for group 2.
#' @param g1.name A character string specifying the name of group 1.
#'   Default is "group1".
#' @param g2.name A character string specifying the name of group 2.
#'   Default is "group2".
#' @param window.size An integer specifying the window width (bp).
#'   Default is `150`.
#' @param window.spacing An integer specifying the distance (bp) between windows.
#'   Default is `50`.
#' @param param A [csaw::readParam-class] object specifying parameters for read counting.
#'   This can be used to restrict chromosomes or exclude intervals (like blacklists).
#'   Default is `readParam(minq = 20, dedup = FALSE)`.
#' @param paired A boolean specifying whether reads are paired-end.
#'   Default is `FALSE`.
#' @param fragment.length An integer specifying the fragment length.
#'   If paired-end, this argument is ignored.
#'   If `NULL`, calculated using [csaw::correlateReads()], which is rather slow.
#'   Default is `200`.
#' @param bg.bin.width An integer specifying the bin width (bp) to calculate background abundance.
#'   Default is `2000`.
#' @param bg.fc A numeric specifying the minimum fold change above background
#'   abundance required to keep windows.
#'   Default is `3`.
#' @param merged.adj.width An integer specifying the maximum distance between adjacent windows.
#'   Adjacent windows within this distance will be merged.
#'   Default is `100`.
#' @param merged.max.width An integer specifying the maximum width (bp) of merged windows to mitigate
#'   excessive daisy chaining.
#'   Default is `30000`.
#' @param BPPARAM A [BiocParallel::BiocParallelParam-class] specifying parallelization strategy.
#'   Default is [BiocParallel::SerialParam-class].
#' @param debug Boolean specifying whether to return a list of intermediate objects for key analysis steps.
#'   Default is `FALSE`.
#'
#' @return A named list containing:
#'   * **merged_regions** - A [GenomicRanges::GRanges-class] object containing regions with associated
#'     differential analysis results including log-fold changes, p-values, FDR, and other statistics.
#'   * **filtered_window_counts** - A [SummarizedExperiment::RangedSummarizedExperiment-class]
#'     object containing filtered window counts, useful for coverage calculation, viz, etc.
#'
#' @importFrom GenomicRanges slidingWindows reduce mcols
#' @importFrom csaw readParam correlateReads maximizeCcf regionCounts windowCounts
#'   filterWindowsGlobal normOffsets asDGEList mergeResults reform
#' @importFrom edgeR estimateDisp glmQLFit glmQLFTest
#' @importFrom BiocParallel SerialParam
#' @importFrom Rsamtools indexBam
#' @importFrom stats relevel model.matrix
#'
#' @export
#'
#' @author Nicolas Peterson, Jared Andrews
#'
#' @examples
#' \dontrun{
#' significant.regions <- find_differential(
#'     regions = my.granges,
#'     g1.bam.files = my.bam1,
#'     g2.bam.files = my.bam2
#' )
#' }
#'
find_differential <- function(regions, g1.bam.files, g2.bam.files,
                              g1.name = "group1",
                              g2.name = "group2",
                              window.size = 150,
                              window.spacing = 50,
                              param = readParam(minq = 20, dedup = FALSE),
                              paired = FALSE,
                              fragment.length = 200,
                              bg.bin.width = 2000,
                              bg.fc = 3,
                              merged.adj.width = 100,
                              merged.max.width = 50000,
                              BPPARAM = SerialParam(),
                              debug = FALSE) {
    # ------------ CHECKS ------------
    # Check if regions is a GRanges object
    if (!inherits(regions, "GRanges")) {
        stop("`regions` must be a GRanges object.")
    }

    bams <- c(g1.bam.files, g2.bam.files)

    # Check if bam.files is a character vector or BamFile object
    if (!is.character(bams) && !inherits(bams, "BamFile")) {
        stop("`bams` must be a character vector of file paths or list of BamFile objects.")
    }

    # Create bam indexes if they don't exist
    for (i in seq_along(bams)) {
        if (is.character(bams[i])) {
            if (!file.exists(paste0(bams[i], ".bai"))) {
                message(paste0("BAM index for ", bams[i], " does not exist, creating."))
                indexBam(bams[i])
            }
        } else if (is(bams[i], "BamFile")) {
            if (length(bams[i]$index) == 0 || !file.exists(bams[i]$index)) {
                message(paste0("BAM index for ", bams[i], " not found. Generating an index."))
                bams[i]$index <- unname(indexBam(bams[i]))
            }
        }
    }

    conditions <- c(
        rep(g1.name, length(g1.bam.files)),
        rep(g2.name, length(g2.bam.files))
    )

    if (!is.factor(conditions)) {
        conditions <- factor(conditions, levels = c(g2.name, g1.name))
    }

    # ------------ WINDOWS ------------
    message("Reducing regions...")
    regions <- reduce(regions)

    # Add a REGION_ID column to keep track of the SEs the bins/windows are derived from
    mcols(regions)$REGION_ID <- seq_along(regions)

    # set up windows for the SEs
    message(paste0("Defining ", window.size, " bp windows spaced ", window.spacing, " bp apart..."))
    region_windows <- slidingWindows(regions, width = window.size, step = window.spacing)

    # determine the number of windows created for each range
    num_region_windows <- elementNROWS(region_windows)

    # map the region IDs to the corresponding window
    region_ids <- rep(mcols(regions)$REGION_ID, num_region_windows)

    # Assign the region IDs
    region_windows <- unlist(region_windows)
    mcols(region_windows)$REGION_ID <- region_ids


    # ------------ COUNT ------------
    # Set up parameters
    my_args <- list(
        bam.files = bams,
        regions = region_windows,
        param = param,
        BPPARAM = BPPARAM
    )

    if (paired == TRUE) {
        message("Paired-end read mode.")
        param <- reform(param, pe = "both")
    } else {
        message("Single-end read mode.")

        # Define the fragment length
        if (is.null(fragment.length)) {
            message("Approximating the fragment length...")
            ccf <- correlateReads(bams, param = param, max.dist = 500, BPPARAM = BPPARAM) # slow
            fragment.length <- maximizeCcf(ccf)
        }
        message(paste0("Fragment length set at ", fragment.length, " bp."))
        my_args$ext <- fragment.length
    }

    message("Counting reads in windows...")
    time1 <- Sys.time()

    # count reads within the windows - NOTE: counts will be different depending on PAIRED == T or F
    window_counts <- do.call(regionCounts, my_args)
    time2 <- Sys.time()
    message(paste0("Time elapsed: ", (time2 - time1)))

    # ------------ BACKGROUND ------------
    message("Calculating background abundance using ", bg.bin.width, " bp bins...")
    time1 <- Sys.time()

    bg_bins <- windowCounts(bams, bin = TRUE, width = bg.bin.width, param = param, BPPARAM = BPPARAM)
    time2 <- Sys.time()
    message(paste0("Time elapsed: ", (time2 - time1)))

    time1 <- Sys.time()
    bg_stat <- filterWindowsGlobal(window_counts, background = bg_bins)
    time2 <- Sys.time()
    message(paste0("Time elapsed: ", (time2 - time1)))

    # ------------ FILTER ------------
    message("Filtering windows using a minimum fold-change over backgound of ", bg.fc, "...")

    keep <- bg_stat$filter > log2(bg.fc)
    filtered_window_counts <- window_counts[keep]

    len_filtered <- length(filtered_window_counts)
    len_original <- length(window_counts)
    pct_filtered <- round((len_filtered / len_original) * 100, 2)

    message(
        "Kept ", len_filtered, " of ", len_original, " windows ",
        "(", pct_filtered, "%)"
    )

    # ------------ NORMALIZATION ------------
    message("Performing normalization...")

    offsets <- normOffsets(filtered_window_counts) # loessfit to address trended bias

    # ------------ DIFFERENTIAL ANALYSIS ------------

    design <- model.matrix(~conditions)
    colnames(design) <- c("intercept", "Condition")

    message("Modeling dispersion...")
    y <- asDGEList(filtered_window_counts)
    y$offset <- assay(offsets, "offset")
    y <- estimateDisp(y, design)

    fit <- glmQLFit(y, design, robust = TRUE)

    message(
        "Performing quasi-likelihood F-Test using levels: ",
        paste(levels(conditions), collapse = " ")
    )

    results <- glmQLFTest(fit)

    rowData(filtered_window_counts) <- cbind(rowData(filtered_window_counts), results$table)

    # ------------ MERGE WINDOWS ------------
    message(
        "Merging adjacent windows within ", merged.adj.width,
        " bp and limiting merged regions to ",
        merged.max.width, " bp..."
    )

    merged_results <- mergeResults(filtered_window_counts,
        mcols(filtered_window_counts),
        tol = merged.adj.width,
        merge.args = list(max.width = merged.max.width)
    )

    message("Created ", length(merged_results$regions), " regions from ", len_filtered, " windows...")

    # get the combined statistics and assign the original REGION_ID
    merged_combined <- merged_results$combined
    merged_combined$REGION_ID <- filtered_window_counts@rowRanges$REGION_ID[merged_combined$rep.test]

    # get the GRanges of the combined windows
    merged_regions <- merged_results$regions
    mcols(merged_regions) <- DataFrame(merged_combined)

    # report regions that have FDR <= 0.05
    sig <- table(merged_regions$FDR <= 0.05)[2]
    total <- length(merged_regions$FDR)
    sig_pct <- round((sig / total) * 100, 2)
    message("Summary: ", sig, " of ", total, " merged regions (", sig_pct, "%) have an FDR <= 0.05.")

    if (debug == TRUE) {
        return(
            list(region_windows, window_counts, filtered_window_counts, merged_results, merged_regions)
        )
    } else {
        return(list(
            merged_regions = merged_regions,
            filtered_window_counts = filtered_window_counts
        ))
    }
}
