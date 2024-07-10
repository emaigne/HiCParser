
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

#' @description
#' Parses interactions in pairs of \code{.matrix} and \code{.bed} files and
#' fills the interactions slots of the provided \code{\link{HiCDOCDataSet}}.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.parseHiCPro <- function(object, replicates, conditions) {

    matrixPaths <- lapply(object@input, `[[`, 1)
    bedPaths <- lapply(object@input, `[[`, 2)

    interactionSet <- pbapply::pbmapply(
        .parseOneHiCPro,
        matrixPaths,
        bedPaths,
        replicates,
        conditions
    )

    mergedinteractionSet <- Reduce(f = .mergeInteractionSet, x = interactionSet)

    object <- new(
        "HiCDOCDataSet",
        mergedinteractionSet,
        input = object@input
    )

    return(object)
}
