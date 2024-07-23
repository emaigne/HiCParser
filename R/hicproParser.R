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
.parseOneHiCPro <- function(matrixPath, bedPath, condition, replicate) {
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
    bed[, start := start + 1]
    # Keeping only intra-chromosomal interactions
    # Add 1 if BED index start with 0
    allChromosomes <- vector("character", length = max(bed$index) + 1)
    allChromosomes[bed[, index] + 1] <- bed[, chromosome]
    interactions <- interactions[
        allChromosomes[startIndex + 1] == allChromosomes[stopIndex + 1]
    ]

    order1 <- match(interactions$startIndex, bed$index)
    order2 <- match(interactions$stopIndex, bed$index)
    allRegions <- GenomicRanges::GRanges(bed)

    gi <- InteractionSet::GInteractions(
        allRegions[order1],
        allRegions[order2],
        regions = allRegions,
        mode = "strict"
    )
    assay <- as.matrix(interactions$interaction, ncol = 1)
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
#' @param conditions
#' A vector of condition names repeated along the replicates.
#' @param replicates
#' A vector of replicate names repeated along the conditions.
#'
#' @return
#' An InteractionSet.
#'
#' @examples
#' # Path to each matrix file
#' matrixPaths <- c(
#'     'path/to/condition-1.replicate-1.matrix',
#'     'path/to/condition-1.replicate-2.matrix',
#'     'path/to/condition-1.replicate-3.matrix',
#'     'path/to/condition-2.replicate-1.matrix',
#'     'path/to/condition-2.replicate-2.matrix',
#'     'path/to/condition-2.replicate-3.matrix'
#' )
#'
#' # Path to each bed file
#' bedPaths <- c(
#'     'path/to/condition-1.replicate-1.bed',
#'     'path/to/condition-1.replicate-2.bed',
#'     'path/to/condition-1.replicate-3.bed',
#'     'path/to/condition-2.replicate-1.bed',
#'     'path/to/condition-2.replicate-2.bed',
#'     'path/to/condition-2.replicate-3.bed'
#' )
#'
#' # Condition and replicate of each file. Can be names instead of numbers.
#' conditions <- c(1, 1, 1, 2, 2, 2)
#' replicates <- c(1, 2, 3, 1, 2, 3)
#'
#' if(FALSE){
#'   # Instantiation of data set
#'   hic.experiment <- parseHiCPro(
#'       matrixPaths = matrixPaths,
#'       bedPaths = bedPaths,
#'       conditions = conditions,
#'       replicates = replicates
#'   )
#' }
#'
#' @usage
#' parseHiCPro(matrixPaths, bedPaths, conditions, replicates)
#'
#' @importFrom pbapply pbmapply
#' @export
parseHiCPro <- function(matrixPaths, bedPaths, conditions, replicates) {
    matrixPaths <- .checkPaths("matrixPaths"=matrixPaths)
    bedPaths <- .checkPaths("bedPaths"=bedPaths)

    if(length(bedPaths) == 1) {
        bedPaths <- rep(bedPaths, length(matrixPaths))
    }
    if (length(matrixPaths) != length(bedPaths)) {
        stop(
            "'matrixPaths' and 'bedPaths' must have the same length.",
            call. = FALSE
        )
    }

    repCond <- .checkConditionsReplicates(conditions, replicates)
    if (min(lengths(repCond)) != length(matrixPaths)) {
        stop(
            "'conditions/replicates' and 'matrixPaths/bedPaths' ",
            "must have the same length",
            call. = FALSE
        )
    }

    interactionSet <- pbapply::pbmapply(
        .parseOneHiCPro,
        matrixPaths,
        bedPaths,
        repCond[["conditions"]],
        repCond[["replicates"]]
    )

    mergedinteractionSet <- Reduce(f = mergeInteractionSet, x = interactionSet)
    mergedinteractionSet <- .sortHiCData(mergedinteractionSet)

    return(mergedinteractionSet)
}
