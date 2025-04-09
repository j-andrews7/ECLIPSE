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
#' 1. Provided super enhancer (`GRanges`) regions are concatenated to form a consensus set.
#' 2. Consensus regions are segmented into sliding windows.
#' 3. Reads from BAM files are counted in each window.
#' 4. Windows with low read counts relative to background are filtered out.
#' 5. Normalization is performed to account for trended bias.
#' 6. Dispersion is modeled using a negative binomial generalized linear model with quasi-likelihood (QL).
#' 7. Differential analysis is performed using a QL F-test to assign statistical significance to each window.
#' 8. Windows are merged into regions based on proximity, and test statistics from individual windows are combined.
#' 9. Differential analysis statistics are added to the metadata of merged regions.
#'
#' @param se.1 A `GRanges` object of ROSE-identified enhancers for the first condition
#' @param se.2 A `GRanges` object of ROSE-identified enhancers for the second condition
#' @param bam.files A character vector of BAM file paths that correspond to each sample.
#' @param conditions A character or factor vector indicating the experimental condition for *each* BAM file.
#' @param window.size An integer specifying the window width (bp).
#'  Default is `150`.
#' @param window.spacing An integer specifying the distance (bp) between windows.
#'  Default is `50`.
#' @param paired A boolean specifying whether reads are paired-end.
#'  Default is `NULL`.
#' @param read.length An integer specifying the read length. If paired-end, this argument is ignored.
#'  If `NULL`, calculated using `correlateReads` from `csaw`. Can also be specified by the user.
#'  Default is `NULL`.
#' @param quality An integer specifying the minimum mapping quality score to include reads.
#'  Default is `20`.
#' @param bg.bin.width An integer specifying the bin width (bp) to calculate background abundance.
#'  Default is `2000`.
#' @param bg.fc A numeric specifying the minimum fold change above background abundance required to keep windows.
#'  Default is `3`.
#' @param merged.adj.width An integer specifying the maximum distance between adjacent windows.
#'  Adjacent windows within this distance will be merged.
#'  Default is `100`.
#' @param merged.max.width An integer specifying the maximum width (bp) of merged windows.
#'  Default is `50000`
#' @param bpparam A \linkS4class{BiocParallelParam} specifying parallelization strategy.
#'  Default is \linkS4class{SerialParam}.
#' @param debug Boolean specifying whether to return a list of intermediate objects for key analysis steps.
#'  Default is `FALSE`.
#'
#' @return A `GRanges` object containing merged super enhancer regions with associated
#' differential analysis results including log-fold changes, p-values, FDR, and other statistics.
#'
#' @importFrom GenomicRanges slidingWindows reduce mcols
#' @importFrom csaw readParam correlateReads maximizeCcf regionCounts windowCounts filterWindowsGlobal normOffsets asDGEList
#' @importFrom edgeR estimateDisp glmQLFit glmQLFTest
#' @importFrom BiocParallel SerialParam
#'
#' @export
#'
#' @author Nicolas Peterson
#'
#' @examples
#
#' \dontrun{
#' significant.regions <- differentialEnhancers(se.1 = ctrl.enhancers,
#'                                              se.2 = treat.enhancers,
#'                                              bam.files = my.bam,
#'                                              conditions = ("ctrl", "ctrl", "drug", "drug")
#'                                              )
#' }
#'
differentialEnhancers <- function(se.1, se.2, bam.files, conditions,
                               window.size = 150, window.spacing = 50,
                               # optional args to readParam for restricting chromosomes or excluding intervals
                               # restrict = NULL, discard = NULL,
                               paired = NULL, read.length = NULL, quality = 20,
                               bg.bin.width = 2000, bg.fc = 3,
                               merged.adj.width = 100, merged.max.width = 50000,
                               bpparam = SerialParam(),
                               debug = FALSE) {

    # ------------ CHECKS ------------

    # Check that the GRanges objects have metadata supplied by ECLIPSE
    meta.names <- unique(c(names(mcols(se.1)), names(mcols(se.2))))

    valid <- "super" %in% meta.names && # check if super exists
        is.logical(se.1$super) && is.logical(se.2$super) # check for boolean values

    if(!valid) {
        stop("GRanges objects must have metadata reflecting their super enhancer status.")
    }

    # check conditions
    if (length(unique(conditions)) != 2) {
        stop("Only two conditions are allowed for differential analysis.")
    }
    if (length(bam.files) != length(conditions)) {
        stop("Each BAM file must have an associated condition.
             Ensure the length of `conditions` matches the number of BAM files provided.")
    }
    if (!is.factor(conditions)) {
        conditions <- factor(conditions) # factors are designated alphabetically
    }

    # check paired-end
    if (is.null(paired)) {
        stop("Please specify whether reads are paired-end or single-end.")
    }

    # ------------ WINDOWS ------------
    # concatenate super enhancer regions (super == T)
    all.SEs <- c(se.1[se.1$super], se.2[se.2$super])

    # make a GRanges object that combines overlapping SE regions
    message("Collecting super enhancer regions...")
    all.SEs <- reduce(all.SEs)

    # Add a REGION_ID column to keep track of the SEs the bins/windows are derived from
    mcols(all.SEs)$REGION_ID <- seq_along(all.SEs)

    # set up windows for the SEs
    message(paste0("Defining ", window.size, " bp windows spaced ", window.spacing, " bp apart..."))
    se.windows <- slidingWindows(all.SEs, width = window.size, step = window.spacing)

    # determine the number of windows created for each range
    se.num.windows <- elementNROWS(se.windows)

    # map the region IDs to the corresponding window
    se.region.ids <- rep(mcols(all.SEs)$REGION_ID, se.num.windows)

    # Assign the region IDs
    se.windows <- unlist(se.windows)
    mcols(se.windows)$REGION_ID <- se.region.ids


    # ------------ COUNT ------------
    # Set up parameters
    my.param <- readParam(minq = quality, dedup = FALSE)

    if (paired == TRUE) {
        message("Paired-end read mode.")
        my.param <- reform(my.param, pe = "both")
    }

    if (paired == FALSE) {
        message("Single-end read mode.")

        # Define the fragment length
        if(is.null(read.length)) {
            message("Approximating the fragment length...")
            ccf <- correlateReads(bam.files, param = my.param, max.dist = 500, BPPARAM = bpparam) # slow
            read.length <- maximizeCcf(ccf)
        }
        message(paste0("Read length set at ", read.length, " bp."))
    }

    message("Counting reads into windows...")
    time1 <- Sys.time()

    my.args <- list(
        bam.files = bam.files,
        regions = se.windows,
        param = my.param,
        BPPARAM = bpparam
    )
    if (paired == FALSE) { my.args$ext <- read.length }

    # count reads within the windows - NOTE: counts will be different depending on PAIRED == T or F
    se.window.counts <- do.call(regionCounts, my.args)
    time2 <- Sys.time()
    message(paste0("Time elapsed: ", (time2-time1)))

    # ------------ BACKGROUND ------------
    message("Calculating background abundance using ", bg.bin.width, " bp bins...")
    time1 <- Sys.time()

    bg.bins <- windowCounts(bam.files, bin = TRUE, width = bg.bin.width, param = my.param, BPPARAM = bpparam)
    time2 <- Sys.time()
    message(paste0("Time elapsed: ", (time2-time1)))

    time1 <- Sys.time()
    bg.stat <- filterWindowsGlobal(se.window.counts, background = bg.bins)
    time2 <- Sys.time()
    message(paste0("Time elapsed: ", (time2-time1)))

    # ------------ FILTER ------------
    message("Filtering windows using a minimum fold-change threshold of ", bg.fc, "...")

    # determine logical values for windows that meet the threshold
    to.keep <- bg.stat$filter > log2(bg.fc)

    # create filtered GRanges
    filtered.by.global <- se.window.counts[to.keep]

    # calculate some stats to report
    len.filtered <- length(filtered.by.global)
    len.original <- length(se.window.counts)
    pct.filtered <- round((len.filtered/len.original)*100,2)

    message("Kept ", len.filtered, " of ", len.original, " windows ",
            "(", pct.filtered, "%)")

    # ------------ NORMALIZATION ------------
    message("Performing normalization...")

    offsets <- normOffsets(filtered.by.global) # loessfit to address trended bias

    # ------------ DIFFERENTIAL ANALYSIS ------------

    design <- model.matrix(~conditions)
    colnames(design) <- c("intercept", "Condition")

    message("Modeling dispersion...")
    y <- asDGEList(filtered.by.global)
    y$offset <- assay(offsets, "offset")
    y <- estimateDisp(y, design)

    fit <- glmQLFit(y, design, robust = TRUE)

    message("Performing quasi-likelihood F-Test using levels: ", paste(levels(conditions), collapse = " "))

    results <- glmQLFTest(fit, coef = "Condition")

    rowData(filtered.by.global) <- cbind(rowData(filtered.by.global), results$table)

    # ------------ MERGE WINDOWS ------------
    message("Merging adjacent windows within ", merged.adj.width, " bp and limiting merged regions to ",
            merged.max.width, " bp...")
    merged.results <- mergeResults(filtered.by.global,
                                   mcols(filtered.by.global),
                                   tol = merged.adj.width,
                                   merge.args = list(max.width = merged.max.width)
                                   )
    message("Created ", length(merged.results$regions), " regions from ", len.filtered, " windows...")

    # get the combined statistics and assign the original REGION_ID
    merged.combined <- merged.results$combined
    merged.combined$REGION_ID <- filtered.by.global@rowRanges$REGION_ID[merged.combined$rep.test]

    # get the GRanges of the combined windows
    merged.regions <- merged.results$regions
    mcols(merged.regions) <- DataFrame(merged.combined)

    # report regions that have FDR <= 0.05
    sig <- table(merged.regions$FDR <= 0.05)[2]
    total <- length(merged.regions$FDR)
    sig.pct <- round((sig/total)*100,2)
    message("Summary: ", sig, " of ", total, " merged regions (", sig.pct, "%) have an FDR <= 0.05.")

    if (debug == TRUE) {
        return(
            list(se.windows, se.window.counts, filtered.by.global, merged.results, merged.regions))
        }
    else { 
        return(merged.regions) 
        }
}


