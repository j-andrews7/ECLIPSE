#' Extend reads directionally upstream and downstream
#'
#' Extends reads by a specified number of bases upstream or downstream
#' based on strand orientation.
#'
#' @param regions A `GRanges` object containing genomic regions.
#' @param upstream Number of bases to extend upstream. 
#'   Default is 0.
#' @param downstream Number of bases to extend downstream. 
#'   Default is 0.
#'
#' @return A `GRanges` object with extended regions.
#'
#' @export
#' @importFrom GenomicRanges strand start end 'ranges<-'
#' @importFrom IRanges IRanges
#'
#' @author Jared Andrews
#'
#' @examples
#' library(GenomicRanges)
#' regions <- GRanges(seqnames = "chr1", ranges = IRanges(100, 200), strand = "+")
#' extended <- extend_reads(regions, upstream = 50, downstream = 100)
#' extended
extend_reads <- function(regions, upstream = 0, downstream = 0) {
    pos <- strand(regions) == "+" | strand(regions) == "*"
    ext_start <- start(regions) - ifelse(pos, upstream, downstream)
    ext_end <- end(regions) + ifelse(pos, downstream, upstream)
    ranges(regions) <- IRanges(ext_start, ext_end)
    regions
}


#' Unstitch regions based on overlaps with TSS
#'
#' This function takes stitched regions and identifies overlaps with TSS regions.
#' If a given stitched region overlaps TSS from more than a specified number of unique genes,
#' the region is unstitched into its original components from the `original` regions.
#'
#' @details
#' This attempts to emulate the method used by ROSE, but the results it returns are not identical due to a presumed bug in ROSE
#' wherein it does not always properly find all overlaps.
#'
#' @param stitched A `GRanges` object representing the stitched regions.
#' @param original A `GRanges` object representing the original regions prior to stitching.
#' @param tss A `GRanges` object representing the TSS regions.
#' @param id.col A character string specifying the column name in the TSS object that contains gene IDs. 
#'   Default is "GENEID".
#' @param threshold An integer specifying the threshold for the number of unique genes before unstitching is applied. 
#'   Default is 2.
#'
#' @return A list with a `hits` data.frame with the following columns:
#' \describe{
#'   \item{qhits}{Query hits corresponding to the stitched regions.}
#'   \item{chrom}{Chromosome names of the stitched regions.}
#'   \item{start}{Start positions of the stitched regions.}
#'   \item{end}{End positions of the stitched regions.}
#'   \item{num_olap}{Number of unique genes overlapping each stitched region.}
#'   \item{ids}{Comma-separated gene IDs overlapping each stitched region.}
#' }
#'
#' The `regions` element of the list contains a GRanges object with the stitched regions
#' and the original regions resulting from unstitching.
#'
#' @importFrom GenomicRanges findOverlaps seqnames start end mcols
#' @importFrom S4Vectors subjectHits
#'
#' @export
#'
#' @author Jared Andrews
#'
#' @examples
#' library(GenomicRanges)
#' 
#' # Example stitched regions, original peaks, and TSSes
#' stitched <- GRanges(seqnames = c("chr1", "chr1", "chr2"),
#'                     ranges = IRanges(start = c(100, 300, 500),
#'                     end = c(200, 400, 600)),
#'                     strand = c("+", "-", "+"))
#' 
#' original <- GRanges(seqnames = c("chr1", "chr1", "chr2", "chr2"),
#'                     ranges = IRanges(start = c(100, 150, 300, 500),
#'                     end = c(150, 200, 350, 550)),
#'                     strand = c("+", "+", "-", "+"))
#' 
#' tss <- GRanges(seqnames = c("chr1", "chr1", "chr1", "chr1", "chr2"),
#'                ranges = IRanges(start = c(120, 140, 160, 320, 520), 
#'                end = c(130, 150, 170, 330, 530)),
#'                strand = c("+", "+", "+", "-", "+"),
#'                GENEID = c("gene1", "gene11", "gene111", "gene2", "gene3"))
#' 
#' result <- unstitch_regions(stitched, original, tss, id.col = "GENEID", threshold = 2)
#' 
#' result$regions
unstitch_regions <- function(stitched, original, tss, id.col = "GENEID", threshold = 2) {
    # Find overlaps between stitched regions and TSS regions
    overlaps <- findOverlaps(stitched, tss)

    # Extract GENEID from TSS overlaps
    overlap_df <- as.data.frame(overlaps)
    overlap_df[[id.col]] <- as.character(mcols(tss)[[id.col]][overlap_df$subjectHits])

    # Calculate the number of unique genes overlapping each stitched region
    unique_genes_per_stitched <- lapply(unique(overlap_df$queryHits), function(x) {
        y <- overlap_df[overlap_df$queryHits == x, ]
        length(unique(y[[id.col]]))
    })

    # Snag the actual symbols as well
    unique_symbols_per_stitched <- lapply(unique(overlap_df$queryHits), function(x) {
        y <- overlap_df[overlap_df$queryHits == x, ]
        unique(y[[id.col]])
    })

    names(unique_genes_per_stitched) <- unique(overlap_df$queryHits)
    names(unique_symbols_per_stitched) <- unique(overlap_df$queryHits)

    hit_df <- data.frame(
        qhits = unique(overlap_df$queryHits),
        chrom = seqnames(stitched[unique(overlap_df$queryHits)]),
        start = start(stitched[unique(overlap_df$queryHits)]),
        end = end(stitched[unique(overlap_df$queryHits)]),
        num_olap = unlist(unique_genes_per_stitched),
        ids = unlist(lapply(unique_symbols_per_stitched, paste, collapse = ","))
    )

    hit_df$unstitch <- FALSE
    hit_df$unstitch[hit_df$num_olap > threshold] <- TRUE

    # Regions that do not need to be unstitched
    need_unstitch <- hit_df$qhits[hit_df$unstitch == TRUE]
    regions_to_keep <- stitched[-need_unstitch]

    # Regions that need to be unstitched
    if (any(need_unstitch)) {
        stitched_to_unstitch <- stitched[need_unstitch]
        overlaps_unstitch <- findOverlaps(stitched_to_unstitch, original)
        regions_to_unstitch <- original[subjectHits(overlaps_unstitch)]
        # Combine the regions
        final_regions <- c(regions_to_keep, regions_to_unstitch)
    } else {
        final_regions <- regions_to_keep
    }

    out <- list(regions = final_regions, hits = hit_df)
    out
}


#' Get region signal from BAM files and add to GRanges
#'
#' Calculates the libary size-normalized coverage of reads within specified regions from a
#' signal BAM file with optional subtraction of a control BAM file.
#' Allows read extension and additional processing like setting negative coverage to zero.
#'
#' @details The defaults of this function are set to so as to match the functionality
#' of ROSE as closely as possible.
#'
#' To detail the process:
#' - The total number of reads in the signal BAM file is calculated.
#' - The signal reads are extended downstream by a specified number of bases (200 bp by default).
#' - The coverage for each basepair in the regions of interest are calculated.
#'     ROSE does this manually by calling samtools for each region, which is slow.
#' - Basepairs with coverage below a specified threshold (`floor`, 1 by default) are removed.
#' - For each region, the coverage is summed and divided by the total number of reads to get the signal.
#'
#' @param sample.bam Character or BamFile object representing the signal BAM file.
#' @param regions Character or GRanges object representing genomic regions of interest.
#' @param control.bam Character or BamFile object representing the control BAM file.
#'   Default is `NULL`.
#' @param floor Numeric value for the minimum coverage threshold to count for region.
#'   Coverage must be greater than this value to be counted in the signal.
#'   Default is 1.
#' @param read.ext Numeric value for extending reads downstream. Default is 200.
#'
#' @return A GRanges object for `regions` with additional columns for sample and control (if provided) signal.
#'   The metadata of the object will also contain the scaling factor for library size normalization for the `sample.bam`
#'   and `control.bam` (if provided) in the `sample_mmr` and `control_mmr` elements.
#'
#' @export
#'
#' @importFrom Rsamtools BamFile idxstatsBam
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges granges trim
#' @importFrom HelloRanges do_bedtools_coverage
#' @importFrom IRanges IRanges
#' @importFrom genomation readBed
#'
#' @author Jared Andrews
#'
#' @examples
#' library(GenomicRanges)
#' \dontrun{
#' sample.bam <- "path/to/sample.bam"
#' regions <- GRanges(seqnames = "chr1", ranges = IRanges(1000, 2000))
#' regions <- add_region_signal(sample.bam, regions, floor = 10)
#' }
add_region_signal <- function(sample.bam,
                              regions,
                              control.bam = NULL,
                              floor = 1,
                              read.ext = 200) {
    if (is.character(sample.bam)) {
        sample.bam <- BamFile(sample.bam)
    }

    if (is.character(control.bam)) {
        control.bam <- BamFile(control.bam)
    }

    if (is.character(regions)) {
        regions <- readBed(file = regions)
    }

    samp_stat <- idxstatsBam(sample.bam)
    samp_total <- sum(samp_stat$mapped)
    samp_mmr <- samp_total / 1000000

    samp_reads <- granges(import(sample.bam))
    samp_reads <- trim(samp_reads)

    if (read.ext > 0) {
        samp_reads <- extend_reads(samp_reads, downstream = read.ext)
    }

    samp_cov <- do_bedtools_coverage(a = regions, b = samp_reads, d = TRUE)

    samp_cov$coverage <- samp_cov$coverage[samp_cov$coverage > floor]
    samp_cov$coverage <- sum(samp_cov$coverage)

    samp_cov$signal <- samp_cov$coverage / samp_mmr

    regions$sample_signal <- samp_cov$signal
    metadata(regions)$sample_mmr <- samp_mmr

    ctrl_sig <- NULL
    if (!is.null(control.bam)) {
        ctrl_stat <- idxstatsBam(control.bam)
        ctrl_total <- sum(ctrl_stat$mapped)
        ctrl_mmr <- ctrl_total / 1000000

        ctrl_reads <- granges(import(control.bam))
        ctrl_reads <- trim(ctrl_reads)
        if (read.ext > 0) {
            ctrl_reads <- extend_reads(ctrl_reads, downstream = read.ext)
        }

        ctrl_cov <- do_bedtools_coverage(a = regions, b = ctrl_reads, d = TRUE)

        ctrl_cov$coverage <- ctrl_cov$coverage[ctrl_cov$coverage > floor]
        ctrl_cov$coverage <- sum(ctrl_cov$coverage)

        ctrl_cov$signal <- ctrl_cov$coverage / ctrl_mmr

        regions$control_signal <- ctrl_cov$signal
        metadata(regions)$control_mmr <- ctrl_mmr
    }

    regions
}


#' Get ranking signal for regions and add to GRanges object
#'
#' Computes the ranking signal for the given genomic regions by
#' optionally subtracting control signal and setting negative values to zero.
#'
#' @param regions A GRanges object containing the `sample_signal` and optionally `control_signal` in metadata columns.
#' @param negative.to.zero Logical indicating whether to set negative values in the ranking signal to zero.
#'   Default is `TRUE`, as that is what ROSE does.
#'
#' @return A GRanges object with an added `rank_signal` column containing the
#'   computed `rank_signal` column in its metadata columns, sorted by said column.
#'   Adds a `region_rank` column as well. Adds a `control_subtracted` element to the `metadata`.
#'
#' @export
#'
#' @importFrom GenomicRanges sort
#'
#' @author Jared Andrews, Jacqueline Myers
#'
#' @examples
#' library(GenomicRanges)
#' regions <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                    ranges = IRanges(start = c(100, 200, 300),
#'                    end = c(200, 300, 400)))
#' regions$sample_signal <- rnorm(length(regions))
#' regions$control_signal <- rnorm(length(regions))
#' ranked_regions <- add_signal_rank(regions)
add_signal_rank <- function(regions, negative.to.zero = TRUE) {
    if (is.null(regions$sample_signal)) {
        stop("regions must contain signal, run 'get_region_signal'")
    }

    rank_sig <- regions$sample_signal

    metadata(regions)$control_subtracted <- FALSE
    if (!is.null(regions$control_signal)) {
        rank_sig <- rank_sig - regions$control_signal
        metadata(regions)$control_subtracted <- TRUE
    }

    if (negative.to.zero) {
        rank_sig[rank_sig < 0] <- 0
    }

    regions$rank_signal <- rank_sig
    regions <- sort(regions, decreasing = TRUE, by = ~rank_signal)
    regions$region_rank <- seq_len(NROW(regions))
    regions
}


#' Classify enhancers based on signal thresholds
#'
#' Classifies enhancers as super enhancers or typical enhancers based on a ranking signal and
#' a specified thresholding method.
#' Optionally applies user-transformations to the signal before classification, which is highly recommended
#' to ameliorate the effects of outliers on the classification.
#'
#' @param regions A GRanges object containing `rank_signal` and optionally other metadata.
#' @param transformation A function to apply to the ranking signal before threshold determination.
#'   Default is `NULL`.
#' @param drop.zeros Logical indicating whether to drop regions with zero signal.
#'   Default is `FALSE`.
#' @param thresh.method Character specifying the method to determine the signal threshold.
#'   Must be one of "ROSE", "first", "curvature", or "arbitrary".
#'   Default is "ROSE".
#' @param first.threshold Numeric value for the fraction of steepest slope when using the "first" threshold method.
#'   Higher values will result in fewer SEs called.
#'   Default is 0.5.
#' @param arbitrary.threshold Numeric value for the arbitrary threshold if the "arbitrary" method is selected.
#'   Default is 0.4, which is a reasonable setting when a cumulative proportion of signal transformation is applied.
#'
#' @return A GRanges object with a new `super` logical column indicating whether the enhancer is classified as a super enhancer.
#'   A `rankby_signal` column is added as the final signal values used for ranking (post-transformation, if applied).
#'   Any transformations applied, the thresholding method used, the threshold, and the number of dropped regions if
#'   `drop.zeros = TRUE` are added to the metadata of the GRanges object.
#'
#' @export
#'
#' @importFrom S4Vectors 'metadata<-'
#' @importFrom KneeArrower findCutoff
#'
#' @author Jared Andrews
#'
#' @examples
#' library(GenomicRanges)
#' regions <- GRanges(seqnames = c("chr1", "chr2", "chr3"),
#'                    ranges = IRanges(start = c(100, 200, 300),
#'                                      end = c(200, 300, 400)))
#' regions$rank_signal <- rnorm(length(regions))
#' classified_regions <- classify_enhancers(regions, thresh.method = "ROSE")
classify_enhancers <- function(regions,
                               transformation = NULL,
                               drop.zeros = FALSE,
                               thresh.method = "ROSE",
                               first.threshold = 0.5,
                               arbitrary.threshold = 0.4) {
    if (is.null(regions$rank_signal)) {
        stop("regions must contain ranking signal, run 'add_signal_rank'")
    }

    # transformation must be a function if provided
    if (!is.null(transformation) && !is.function(transformation)) {
        stop("transformation must be a function")
    }

    if (!thresh.method %in% c("ROSE", "first", "curvature", "arbitrary")) {
        stop("thresh.method must be one of 'ROSE', 'first', 'curvature', or 'arbitrary'")
    }

    metadata(regions)$threshold_method <- thresh.method

    if (drop.zeros) {
        num_no_sig_regions <- NROW(regions[regions$rank_signal == 0])
        message(paste("Dropped", num_no_sig_regions, "regions due to no signal"))
        metadata(regions)$dropped_zero_count <- num_no_sig_regions
        regions <- regions[regions$rank_signal > 0]
        regions$region_rank <- seq_len(NROW(regions))
    }

    # Keep track of whether to use transformed signal.
    use_transformed <- FALSE
    if (!is.null(transformation)) {
        message("Applying provided transformation to signal prior to determining cutoff")
        metadata(regions)$transformation <- transformation
        regions$transformed_signal <- transformation(regions$rank_signal)
        use_transformed <- TRUE
    }

    if (use_transformed) {
        regions$rankby_signal <- regions$transformed_signal
    } else {
        regions$rankby_signal <- regions$rank_signal
    }

    # Use y-axis position for threshold, i.e. the signal value rather than rank
    if (thresh.method == "ROSE") {
        cutoff_options <- calculate_cutoff(regions$rankby_signal, drawPlot = FALSE)
        cutpoint <- cutoff_options$absolute
    } else if (thresh.method == "first") {
        cutpoint <- findCutoff(
            rank(regions$rankby_signal),
            regions$rankby_signal,
            method = "first",
            frac.of.steepest.slope = first.threshold
        )
        cutpoint <- cutpoint$y
    } else if (thresh.method == "curvature") {
        cutpoint <- findCutoff(rank(regions$rankby_signal),
            regions$rankby_signal,
            method = "curvature"
        )
        cutpoint <- cutpoint$y
    } else if (thresh.method == "arbitrary") {
        cutpoint <- arbitrary.threshold
    }

    metadata(regions)$threshold <- cutpoint
    message(paste("Using", cutpoint, "as cutoff for SE classification"))

    regions$super <- regions$rankby_signal > cutpoint
    message(paste(sum(regions$super), "super enhancers called"))

    regions
}


#' Run ROSE (Rank Ordering of Super-Enhancers)
#'
#' This function performs the ROSE for identifying super-enhancers by stitching
#' together peaks, calculating the signal in the regions, ranking the regions by signal,
#' and classifying them as super-enhancers.
#'
#' @details
#' This function allows for near identical functionality to the original ROSE software,
#' but also provides a number of additional options to minimize the impact of signal outliers
#' on the classification threshold.
#' In particular, the `transformation` argument allows for
#' the application of a function to the ranking signal before threshold determination.
#' This is highly recommended to prevent outliers from skewing the threshold.
#'
#' The `thresh.method` argument allows for the selection of the method to determine the threshold.
#' - The "ROSE" method is the default and uses a sliding diagonal line to determine the cutoff.
#' - The "first" method uses a first derivative to determine the point where the slope is a given fraction of the maximum.
#'   The `first.threshold` argument controls this fraction.
#' - The "curvature" method finds the point at which the circle tanget to the curve has the smallest radius.
#' - The "arbitrary" method allows for the user to specify a fixed threshold, useful for transformations
#'   that result in a consistent curve shape with a known maximum, like cumulative proportion of signal.
#'
#' @param treatment A character string or `BamFile` object representing the sample BAM file.
#' @param peaks A character string or `GRanges` object representing the peaks.
#' @param control A character string or `BamFile` object representing the control BAM file.
#'   Default is `NULL`.
#' @param stitch.distance Numeric value for the distance within which peaks are stitched together.
#'   Default is 12500.
#' @param tss.exclusion.distance Numeric value for distance to add to TSS exclusion range prior to stitching.
#'   Peaks *fully contained* within this exclusion range will be excluded from the stitching process.
#'   Default is 0, meaning peaks overlapping TSS will not be excluded from stitching.
#' @param txdb Transcript annotations used for peak exclusion for TSS overlap and unstitching process for
#'   regions spanning TSSes from more than `max.unique.gene.tss.overlap` genes. Also used for enhancer annotation.
#' @param org.db An `OrgDb` object containing organism database information. Used for enhancer annotation.
#'   Default is `NULL`.
#' @param drop.y Logical indicating whether to drop peaks on chromosome Y, as done by the original ROSE implementation.
#'   Default is `TRUE`.
#' @param max.unique.gene.tss.overlap Maximum number of unique genes that a region can overlap the TSS of before being unstitched.
#'   Note that multiple overlapping regions from the same gene, e.g. multiple isoforms, are counted as one.
#'   Default is `NULL`. Ignored if `txdb` is `NULL`.
#' @param tss.overlap.distance Numeric value for the distance (in bp) to add to TSS for unstitching.
#'   Default is 50. Ignored if `txdb` is `NULL`.
#' @param negative.to.zero Logical indicating whether to set negative ranking signals to zero.
#'   Default is `TRUE`.
#' @param thresh.method Character string specifying the method to determine the signal threshold.
#'   Must be one of "ROSE", "first", "curvature", or "arbitrary".
#'   Default is "ROSE".
#' @param transformation A function to apply to the ranking signal before threshold determination.
#'   Default is `NULL`.
#' @param floor Numeric value representing the minimum coverage threshold to count.
#'  Default is 1.
#' @param read.ext Numeric value for extending reads downstream.
#'   Default is 200. Ignored if inputs are `GRanges` objects.
#' @param drop.zeros Logical indicating whether to drop regions with zero signal.
#'   Default is `FALSE`.
#' @param first.threshold Numeric value for the fraction of steepest slope when using the "first" threshold method.
#'   Default is 0.5.
#' @param arbitrary.threshold Numeric value for the arbitrary threshold if the "arbitrary" method is selected.
#'   Default is 0.4.
#' @param annotate Logical indicating whether annotations should be provided.
#'   Default is `TRUE`.
#' @param annotate.dist Numeric specifying the flanking distance (in bp) around the region to annotate.
#'   Default is 50000.
#' @param promoter.dist Integer vector of length 2, where first and second values are upstream/downstream
#'   distances (in bp) from the transcription start site. Used to specify the promoter region for annotations
#'   and determining active genes.
#'   Default is `c(2000, 200)`.
#' @param active.genes Character vector of gene *symbols* to retain in the annotations.
#'   Default is `NULL`.
#' @param identify.active.genes Logical indicating whether active genes should be identified based on
#'   overlaps with `peaks`. Cannot be used simultaneously with `active.genes`.
#'   Default is `FALSE`.
#' @param omit.unkown Logical indicating whether uncharacterized genes (i.e., gene symbols starting with *LOC*)
#'   should be excluded from annotations.
#'   Default is `TRUE`.
#'
#' @return A `GRanges` object containing the classified super-enhancers and associated metadata.
#'
#' @author Jared Andrews, Nicolas Peterson
#'
#' @importFrom Rsamtools BamFile indexBam
#' @importFrom genomation readBed
#' @importFrom GenomicRanges reduce seqnames trim
#' @importFrom S4Vectors queryHits
#' @importFrom IRanges findOverlaps
#'
#' @export
#'
#' @examples
#' \dontrun{
#' sample.bam <- "path/to/sample.bam"
#' peaks <- "path/to/peaks.bed"
#' result <- run_rose(sample.bam, peaks)
#' }
run_rose <- function(
    treatment,
    peaks,
    control = NULL,
    stitch.distance = 12500,
    tss.exclusion.distance = 0,
    txdb = NULL,
    org.db = NULL,
    drop.y = TRUE,
    max.unique.gene.tss.overlap = NULL,
    tss.overlap.distance = 50,
    negative.to.zero = TRUE,
    thresh.method = "ROSE",
    transformation = NULL,
    floor = 1,
    read.ext = 200,
    drop.zeros = FALSE,
    first.threshold = 0.5,
    arbitrary.threshold = 0.4,
    annotate = TRUE,
    annotate.dist = 50000,
    promoter.dist = c(2000, 200),
    active.genes = NULL,
    identify.active.genes = FALSE,
    omit.unkown = TRUE) {

    if (is.character(treatment)) {
        treatment <- BamFile(treatment)
    }
    message(paste0("Treatment BAM file: ", sub(".*/(.*\\.bam)$", "\\1", treatment$path)))

    if (length(treatment$index) == 0 || !file.exists(treatment$index)) {
        message("Treatment BAM index not found. Generating an index.")
        treatment$index <- unname(indexBam(treatment))
    }

    if (!is.null(control)) {
        if (is.character(control)) {
            control <- BamFile(control)
        }
        message(paste0("Control BAM file: ", sub(".*/(.*\\.bam)$", "\\1", control$path)))

        if (length(control$index) == 0 || !file.exists(control$index)) {
            message("Control BAM index not found. Generating an index.")
            control$index <- unname(indexBam(control))
        }
    }

    message("Reading peaks")
    if (is.character(peaks)) {
        peaks <- readBed(peaks)
        peaks <- trim(peaks)
    }

    if (tss.exclusion.distance > 0) {
        if (is.null(txdb)) {
            stop("txdb must be provided if tss.exclusion.distance is greater than 0")
        }

        message("Excluding peaks within TSS exclusion distance of ", tss.exclusion.distance)
        tss <- promoters(txdb, upstream = tss.exclusion.distance, downstream = tss.exclusion.distance, columns = "GENEID")

        # Create susbset of stitched peaks that do not overlap TSS and remove
        overlaps <- findOverlaps(peaks, tss, type = "within")

        contained_indices <- queryHits(overlaps)

        message(length(contained_indices), " peaks fully contained within TSS exclusion window and will be excluded from stitching")
        peaks <- peaks[-unique(contained_indices)]
    }

    message("Stitching peaks with stitch distance of ", stitch.distance)
    peaks_stitched <- reduce(peaks, min.gapwidth = stitch.distance)

    # Drop chrY as ROSE does
    if (drop.y) {
        peaks.chr <- as.vector(seqnames(peaks_stitched))
        message("Dropped ", length(which(peaks.chr == "chrY")), " peaks on chrY")
        peaks_stitched <- peaks_stitched[peaks.chr != "chrY"]
    }

    if (!is.null(max.unique.gene.tss.overlap)) {
        if (is.null(txdb)) {
            stop("txdb must be provided if max.unique.gene.tss.overlap is not NULL")
        }

        tss <- promoters(txdb, upstream = tss.overlap.distance, downstream = tss.overlap.distance, columns = "GENEID")
        message("Unstitching regions overlapping TSS from more than ", max.unique.gene.tss.overlap, " unique genes")
        unstitched <- unstitch_regions(peaks_stitched, peaks, tss, threshold = max.unique.gene.tss.overlap)
        peaks_stitched <- unstitched$regions
        hits <- unstitched$hits
        message("Unstitched ", sum(hits$unstitch), " regions")
    }

    message("Calculating normalized signal for ", length(peaks_stitched)," stitched regions")
    regions <- add_region_signal(treatment, peaks_stitched, control.bam = control, floor = floor, read.ext = read.ext)

    message("Ranking regions")
    regions <- add_signal_rank(regions, negative.to.zero = negative.to.zero)

    message("Classifying enhancers")
    regions <- classify_enhancers(regions,
        transformation = transformation, drop.zeros = drop.zeros,
        thresh.method = thresh.method, first.threshold = first.threshold, arbitrary.threshold = arbitrary.threshold
    )

    if(annotate) {
        message("Annotating regions")
        regions <- annotate_enhancers(regions, peaks, tx.db = txdb, org.db = org.db, annotate.dist = annotate.dist,
                                      promoter.dist = promoter.dist, active.genes = active.genes,
                                      identify.active.genes = identify.active.genes, omit.unknown = omit.unkown)
    }

    regions
}
