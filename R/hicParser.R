#' @description
#' Parses a single interactions file in \code{.hic} format. Calls the C++
#' \code{parseHiCFile} parser.
#'
#' @param path
#' The path to the interactions file.
#' @param binSize
#' The resolution (span of each position in number of bases) to select within
#' the \code{.hic} file.
#'
#' @return
#' An \code{\link{InteractionSet}}.
#'
#' @keywords internal
#' @noRd
#' @importFrom data.table setalloccol
.parseOneHiC <- function(path, binSize, condition, replicate) {
    message("\nParsing '", path, "'.")
    interactions <- parseHiCFile(path, binSize)
    # Automagical stuff to transform Rcpp DataFrame to data.table
    interactions <- data.table::setalloccol(interactions)
    # Create InteractionSet object
    interactions <- .setFromTabular(interactions, condition, replicate)
    return(interactions)
}

#' @title Parser for data in hic format
#' @description
#' Parses interactions in \code{.hic} format and returns an InteractionSet
#' object.
#'
#' @param paths
#' A vector of paths to \code{.hic} files.
#' @param binSize
#' The resolution (span of each position in number of bases) to select within
#' the \code{.hic} files.
#' @param replicates
#' A vector of replicate names repeated along the conditions.
#' @param conditions
#' A vector of condition names repeated along the replicates.
#'
#' @return
#' An InteractionSet.
#'
#' @examples
#' #Path to each file
#' paths <- c(
#'     "path/to/condition-1.replicate-1.hic",
#'     "path/to/condition-1.replicate-2.hic",
#'     "path/to/condition-2.replicate-1.hic",
#'     "path/to/condition-2.replicate-2.hic",
#'     "path/to/condition-3.replicate-1.hic"
#' )
#' # Replicate and condition of each file. Can be names instead of numbers.
#' replicates <- c(1, 2, 1, 2, 1)
#' conditions <- c(1, 1, 2, 2, 3)
#' # Resolution to select
#' binSize <- 500000
#' if(FALSE){
#' # Instantiation of data set
#' hic.experiment <- parseHiC(
#'     paths,
#'     replicates = replicates,
#'     conditions = conditions,
#'     binSize = binSize
#' )
#' }
#' @usage
#' parseHiC(paths, binSize, replicates, conditions)
#'
#' @importFrom pbapply pbmapply
#' @export
parseHiC <- function(paths, binSize, replicates, conditions) {
    paths <- .checkPaths("paths"=paths)
    repCond <- .checkReplicatesConditions(replicates, conditions)
    if (min(lengths(repCond)) != length(paths)) {
        stop(
            "'conditions/replicates' and 'paths' ",
            "must have the same length",
            call. = FALSE
        )
    }

    if (!is.numeric(binSize) || length(binSize) != 1) {
        stop("'binSize' must be an integer.", call. = FALSE)
    }
    binSize <- as.integer(binSize)

    interactionSet <- pbapply::pbmapply(
        .parseOneHiC,
        path      = paths,
        binSize   = binSize,
        condition = repCond[["conditions"]],
        replicate = repCond[["replicates"]]
    )

    mergedinteractionSet <- Reduce(
        f = mergeInteractionSet,
        x = interactionSet
    )
    mergedinteractionSet <- .sortHiCData(mergedinteractionSet)

    return(mergedinteractionSet)
}
