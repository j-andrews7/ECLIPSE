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
#'   model with quasi-likelihood (QL).
#' 6. Differential analysis is performed using a QL F-test to
#'   assign statistical significance to each window.
#' 7. Windows are merged into regions based on proximity, and test statistics
#'   from individual windows are combined.
#' 8. Differential analysis statistics are added to the metadata of merged regions.
#'
#' @param regions A `GRanges` object of regions to compare between conditions.
#'   These regions will be reduced, i.e. overlapping regions will be merged to one.
#' @param bam.files A character vector of BAM file paths that correspond to each sample.
#' @param conditions A character or factor vector indicating the experimental condition for *each* BAM file.
#' @param window.size An integer specifying the window width (bp).
#'   Default is `150`.
#' @param window.spacing An integer specifying the distance (bp) between windows.
#'   Default is `50`.
#' @param paired A boolean specifying whether reads are paired-end.
#'   Default is `FALSE`.
#' @param fragment.length An integer specifying the fragment length.
#'   If paired-end, this argument is ignored.
#'   If `NULL`, calculated using `correlateReads` from `csaw`, which is rather slow.
#'   Default is `200`.
#' @param quality An integer specifying the minimum mapping quality score to include reads.
#'   Default is `20`.
#' @param bg.bin.width An integer specifying the bin width (bp) to calculate background abundance.
#'   Default is `2000`.
#' @param bg.fc A numeric specifying the minimum fold change above background
#'   abundance required to keep windows.
#'   Default is `3`.
#' @param merged.adj.width An integer specifying the maximum distance between adjacent windows.
#'   Adjacent windows within this distance will be merged.
#'   Default is `100`.
#' @param merged.max.width An integer specifying the maximum width (bp) of merged windows.
#'   Default is `50000`
#' @param BPPARAM A \linkS4class{BiocParallelParam} specifying parallelization strategy.
#'   Default is \linkS4class{SerialParam}.
#' @param debug Boolean specifying whether to return a list of intermediate objects for key analysis steps.
#'   Default is `FALSE`.
#'
#' @return A `GRanges` object containing regions with associated
#' differential analysis results including log-fold changes, p-values, FDR, and other statistics.
#'
#' @importFrom GenomicRanges slidingWindows reduce mcols
#' @importFrom csaw readParam correlateReads maximizeCcf regionCounts windowCounts
#'   filterWindowsGlobal normOffsets asDGEList mergeResults reform
#' @importFrom edgeR estimateDisp glmQLFit glmQLFTest
#' @importFrom BiocParallel SerialParam
#'
#' @export
#'
#' @author Nicolas Peterson, Jared Andrews
#'
#' @examples
#' \dontrun{
#' significant.regions <- find_differential(regions = my.granges,
#'                                          bam.files = my.bam,
#'                                          conditions = ("ctrl", "ctrl", "drug", "drug")
#' )
#' }
#'
find_differential <- function(regions, g1.bam.files, g2.bam.files,
                              g1.name, g2.name,
                              window.size = 150,
                              window.spacing = 50,
                               # optional args to readParam for restricting chromosomes or excluding intervals
                               # restrict = NULL, discard = NULL,
                              paired = FALSE,
                              fragment.length = 200,
                              quality = 20,
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
        stop("`bam.files` must be a character vector of file paths or list of BamFile objects.")
    }

    conditions <- c(rep(g1.name, length(g1.bam.files)),
                    rep(g2.name, length(g2.bam.files)))

    if (!is.factor(conditions)) {
        conditions <- factor(conditions,
                             levels = c(g1.name, g2.name),
                             ordered = TRUE)
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
    my_param <- readParam(minq = quality, dedup = FALSE)

    my_args <- list(
        bam.files = bam.files,
        regions = region_windows,
        param = my_param,
        BPPARAM = BPPARAM
    )

    if (paired == TRUE) {
        message("Paired-end read mode.")
        my_param <- reform(my_param, pe = "both")
    } else {
        message("Single-end read mode.")

        # Define the fragment length
        if(is.null(fragment.length)) {
            message("Approximating the fragment length...")
            ccf <- correlateReads(bam.files, param = my_param, max.dist = 500, BPPARAM = BPPARAM) # slow
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

    bg_bins <- windowCounts(bam.files, bin = TRUE, width = bg.bin.width, param = my_param, BPPARAM = BPPARAM)
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

    message("Kept ", len_filtered, " of ", len_original, " windows ",
            "(", pct_filtered, "%)")

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

    message("Performing quasi-likelihood F-Test using levels: ",
            paste(levels(conditions), collapse = " "))

    results <- glmQLFTest(fit, coef = "Condition")

    rowData(filtered_window_counts) <- cbind(rowData(filtered_window_counts), results$table)

    # ------------ MERGE WINDOWS ------------
    message("Merging adjacent windows within ", merged.adj.width,
            " bp and limiting merged regions to ",
            merged.max.width, " bp...")

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
    sig_pct <- round((sig/total)*100,2)
    message("Summary: ", sig, " of ", total, " merged regions (", sig_pct, "%) have an FDR <= 0.05.")

    if (debug == TRUE) {
        return(
            list(region_windows, window_counts, filtered_window_counts, merged_results, merged_regions)
        )
    } else {
        return(merged_regions)
    }
}
