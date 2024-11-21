#' Annotate genomic ranges using TxDb and OrgDB objects
#'
#' Annotates a `GRanges` object with gene information provided by a `TxDb` and `OrgDb` object.
#' It can also filter for active genes or remove uncharacterized genes.
#'
#' @details
#' Annotation is performed using `detailRanges` from the `csaw` package, which
#' uses information provided by a `TxDb` and `OrgDb` object to identify
#' genes that either overlap the specified regions or are located within
#' a user-defined distance around each region. To filter the annotations,
#' users can specify a vector of active genes (symbols) or identify them
#' automatically based on overlaps with peaks and promoters (`TxDb`). Both
#' unfiltered and filtered annotations are appended to the metadata.
#'
#' @param regions A `GRanges` object representing genomic regions to annotate.
#' @param peaks A `GRanges` object representing peaks to identify active promoters.
#' @param tx.db A `TxDb` object containing transcript database information.
#'   Default is `NULL`.
#' @param org.db An `OrgDb` object containing organism database information.
#'   Default is `NULL`.
#' @param annotate.dist Numeric specifying the flanking distance (in bp) around the region to annotate.
#'   Default is 50000.
#' @param promoter.dist Integer vector of length 2, where first and second values are upstream/downstream
#'   distances (in bp) from the transcription start site. Used to specify the promoter region for annotations
#'   and determining active genes.
#'   Default is `c(2000, 200)`.
#' @param annotate.ids String specifying the gene identifier for annotations.
#'   Default is `SYMBOL` (`ENTREZID` and `ENSEMBL` are not yet supported).
#' @param active.genes Character vector of gene *symbols* to retain in the annotations.
#'   Default is `NULL`.
#' @param identify.active.genes Logical indicating whether active genes should be identified based on
#'   overlaps with `peaks`. Cannot be used simultaneously with `active.genes`.
#'   Default is `FALSE`.
#' @param omit.unkown Logical indicating whether uncharacterized genes (i.e., gene symbols starting with *LOC*)
#'   should be excluded from annotations.
#'   Default is `TRUE`.
#'
#' @return A `GRanges` object with additional metadata columns for annotations.
#'
#' @importFrom csaw detailRanges
#' @importFrom GenomicRanges findOverlaps promoters
#' @importFrom stringi stri_split_fixed stri_replace_first_regex stri_join
#' @importFrom data.table %chin%
#'
#' @export
#'
#' @author Nicolas Peterson
#'
#' @examples
#' # Example usage
#' \dontrun{
#' annotated.regions <- annotate_enhancers(regions = my.enhancers, peaks = my.peaks, tx.db = TxDb,
#'                                         org.db = OrgDb, identify.active.genes = TRUE,
#'                                         omit.unknown = TRUE)
#' }
annotate_enhancers <- function(regions, peaks, tx.db, org.db, annotate.dist = 50000,
                               promoter.dist = c(2000, 200), annotate.ids = "SYMBOL",
                               #annotate.ids = c("SYMBOL", "ENTREZID", "ENSEMBL"),
                               active.genes = NULL, identify.active.genes = FALSE,
                               omit.unknown = TRUE) {

    # Check that TxDb and OrgDb objects are provided
    if(class(tx.db) != "TxDb" || class(org.db) != "OrgDb") {
        stop("TxDb and OrgDb objects must be provided to perform annotation.")
    }

    # Check for mutually exclusive arguments
    if (identify.active.genes && !is.null(active.genes)) {
        stop("Please specify only one of 'identify.active.genes' or 'active.genes', not both.")
    }

    # Check for appropriate vector length
    if (length(promoter.dist) != 2) {
        stop("Please ensure `promoter.dist` is an integer vector of length 2.")
    }

    # Ensure correct gene identifier choice
    #annotate.ids <- match.arg(annotate.ids)

    # Annotate regions, returns a list of character vectors
    annotations <- detailRanges(incoming = regions, txdb = tx.db, orgdb = org.db, dist = annotate.dist,
                                promoter = promoter.dist, name.field = annotate.ids)

    # Identify active genes, if required
    if(identify.active.genes && is.null(active.genes)) {
        message("Identifying active genes")

        # Retrieve promoters
        promoters <- suppressWarnings(promoters(TxDb, upstream = promoter.dist[1], downstream = promoter.dist[2])) # suppress out-of-bounds warning
        promoters <- trim(promoters) # trim out of bounds ranges

        # Identify active promoters by checking overlaps with peaks
        overlaps <- findOverlaps(query = peaks, subject = promoters)
        promoters.with.peaks <- promoters[subjectHits(overlaps)] # active promoters

        # Get unique transcript IDs from active promoters
        tx.id <- as.character(promoters.with.peaks$tx_id)
        tx.id <- unique(tx.id)

        # Get unique gene IDs
        gene.ids <- suppressMessages(
            select(TxDb, keys = tx.id, columns = "GENEID", keytype = "TXID"))
        unique.gene.ids <- unique(na.omit(gene.ids$GENEID))

        # Map gene IDs to symbols to create active.genes
        active.genes <- suppressMessages(
            mapIds(OrgDb, keys = unique.gene.ids, column = "SYMBOL", keytype = "ENTREZID"))
    }

    # initialize for concatenation
    filtered <- NULL

    if(!is.null(active.genes) || omit.unknown) {
        message("Filtering annotations")

        filtered <- lapply(annotations, function(x) { # x is a list of character vectors

            # Split each string in x to isolate gene info, Gene:Strand:Type or Gene:Strand:Dist
            split <- stri_split_fixed(x, ",")

            # Extract gene names
            genes <- lapply(split, stri_replace_first_regex, pattern = ":.*", replacement = "")

            # Set up a list of logical vectors based on user-specified args

            # -- Evaluates TRUE if gene is found in active.genes
            # -- If active.genes is not provided, all TRUE
            if (!is.null(active.genes)) {
                is.active <- lapply(genes, FUN = `%chin%`, active.genes)
            } else {
                is.active <- lapply(genes, function(x) TRUE)
            }

            # -- If omit.unknown is TRUE and gene starts with "LOC", evaluates TRUE
            # -- If omit.unknown is FALSE, all false
            if (omit.unknown) {
                is.uncharacterized <- lapply(genes, FUN = grepl, pattern = "^LOC[0-9]+")
            } else {
                is.uncharacterized <- lapply(genes, function(x) FALSE)
            }

            # Determine which genes to keep by comparing the lists of logical vectors
            to.keep <- Map(function(active, uncharacterized) { active & !uncharacterized },
                           is.active, is.uncharacterized)

            # Filter genes in `split` using to.keep
            kept <- Map(function(genes, keep) genes[keep], split, to.keep)

            # rejoin filtered gene info
            lapply(kept, stri_join, collapse = ",")
        })
        # Convert nested list to a single list of vectors (like annotations)
        filtered <- lapply(filtered, unlist)

        # Update names to distinguish the filtered annotations
        names(filtered) <- paste0("filtered_", names(filtered))
    }

    message("Adding annotations to metadata")

    # Concatenate the unfiltered and filtered annotations
    combined <- c(annotations, filtered)

    # Add each vector to metadata
    for (name in names(combined)) {
        mcols(regions)[[name]] <- combined[[name]]
    }

    return(regions)
}
