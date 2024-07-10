
#' @description
#' Parses a single pair of \code{.matrix} and \code{.bed} files.
#'
#' @param matrixPath
#' The path to the interactions matrix file.
#' @param bedPath
#' The path to the bed file.
#'
#' @return
#' A data.table of interactions.
#'
#' @keywords internal
#' @noRd
#' @importFrom InteractionSet GInteractions
#' @importFrom GenomicRanges GRanges
#' @importFrom data.table setorder
.parseOneHiCPro <- function(matrixPath, bedPath, replicate, condition) {

    message("\nParsing '", matrixPath, "' and '", bedPath, "'.")

    interactions <- data.table::fread(
        matrixPath,
        header = FALSE,
        stringsAsFactors = FALSE,
        col.names = c("startIndex", "stopIndex", "interaction"),
        data.table = TRUE
    )

    bed <- data.table::fread(
        bedPath,
        header = FALSE,
        stringsAsFactors = FALSE,
        col.names = c("chromosome", "start", "end", "index"),
        data.table = TRUE
    )

    setorder(bed, chromosome, start, end)
    # Adding 1 to follow Bioconductor GRanges recommended format
    bed[,start := start+1]
    # Keeping only intra-chromosomal interactions
    # Add 1 if BED index start with 0
    allChromosomes <- vector("character", length = max(bed$index) + 1)
    allChromosomes[bed[,index]+1] <- bed[,chromosome]
    interactions <- interactions[
        allChromosomes[startIndex + 1] == allChromosomes[stopIndex + 1]]

    order1 <- match(interactions$startIndex, bed$index)
    order2 <- match(interactions$stopIndex, bed$index)
    allRegions <- GenomicRanges::GRanges(bed)

    gi <- InteractionSet::GInteractions(
        allRegions[order1],
        allRegions[order2],
        regions = allRegions,
        mode="strict"
    )
    assay <- as.matrix(interactions$interaction, ncol=1)
    interactionSet <- .createInteractionSet(assay, gi, allRegions, condition, replicate)

    return(interactionSet)
}

#' @title Parser for HiCPro data
#' @description
#' Parses interactions in pairs of \code{.matrix} and \code{.bed} files and
#' returns an InteractionSet object.
#'
#' @param matrixPaths
#' A vector of paths to HiC-Pro matrix files.
#' @param bedPaths
#' A vector of paths to HiC-Pro bed files.
#' @param replicates
#' A vector of replicate names repeated along the conditions.
#' @param conditions
#' A vector of condition names repeated along the replicates.
#'
#' @return
#' An InteractionSet.
#'
#' @examples
#' \dontrun{
#'     # Path to each matrix file
#'     matrixPaths = c(
#'       'path/to/condition-1.replicate-1.matrix',
#'       'path/to/condition-1.replicate-2.matrix',
#'       'path/to/condition-2.replicate-1.matrix',
#'       'path/to/condition-2.replicate-2.matrix',
#'       'path/to/condition-3.replicate-1.matrix'
#'     )
#'
#'     # Path to each bed file
#'     bedPaths = c(
#'       'path/to/condition-1.replicate-1.bed',
#'       'path/to/condition-1.replicate-2.bed',
#'       'path/to/condition-2.replicate-1.bed',
#'       'path/to/condition-2.replicate-2.bed',
#'       'path/to/condition-3.replicate-1.bed'
#'     )
#'
#'     # Replicate and condition of each file. Can be names instead of numbers.
#'     replicates <- c(1, 2, 1, 2, 1)
#'     conditions <- c(1, 1, 2, 2, 3)
#'
#'     # Instantiation of data set
#'     hic.experiment <- HiCDOCDataSetFromHiCPro(
#'       matrixPaths = matrixPaths,
#'       bedPaths = bedPaths,
#'       replicates = replicates,
#'       conditions = conditions
#'     )
#' }
#'
#' @usage
#' HiCDOCDataSetFromHiCPro(matrixPaths, bedPaths, replicates, conditions)
#'
#' @importFrom pbapply pbmapply
#' @export
parseHiCPro <- function(matrixPaths, bedPaths, replicates, conditions) {

    if (is.factor(matrixPaths)) {
        matrixPaths <- as.vector(matrixPaths)
    }
    if (!is.character(matrixPaths)) {
        stop("'matrixPaths' must be a vector of characters.", call. = FALSE)
    }

    if (is.factor(bedPaths)) {
        bedPaths <- as.vector(bedPaths)
    }
    if (!is.character(bedPaths)) {
        stop("'bedPaths' must be a vector of characters.", call. = FALSE)
    }

    if (length(matrixPaths) != length(bedPaths)) {
        stop(
            "'matrixPaths' and 'bedPaths' must have the same length.",
            call. = FALSE
        )
    }

    for (path in c(matrixPaths, bedPaths)) {
        if (!file.exists(path)) {
            stop("'", path, "' does not exist.", call. = FALSE)
        }
    }

    if (is.factor(replicates)) {
        replicates <- as.vector(replicates)
    }
    if (is.null(replicates)) {
        stop("'replicates' must be a vector of replicates.", call. = FALSE)
    }

    if (is.factor(conditions))
        conditions <- as.vector(conditions)
    if (is.null(conditions)) {
        stop("'conditions' must be a vector of conditions.", call. = FALSE)
    }

    if (length(conditions) != length(replicates)) {
        stop(
            "'conditions' and 'replicates' must have the same length",
            call. = FALSE
        )
    }

    interactionSet <- pbapply::pbmapply(
        .parseOneHiCPro,
        matrixPaths,
        bedPaths,
        replicates,
        conditions
    )

    mergedinteractionSet <- Reduce(f = mergeInteractionSet, x = interactionSet)
    mergedinteractionSet <- .sortHiCData(mergedinteractionSet)

    return(mergedinteractionSet)
}
