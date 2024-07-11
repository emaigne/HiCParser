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
#' \dontrun{
#' #' # Path to each file
#' paths <- c(
#'     "path/to/condition-1.replicate-1.hic",
#'     "path/to/condition-1.replicate-2.hic",
#'     "path/to/condition-2.replicate-1.hic",
#'     "path/to/condition-2.replicate-2.hic",
#'     "path/to/condition-3.replicate-1.hic"
#' )
#'
#' # Replicate and condition of each file. Can be names instead of numbers.
#' replicates <- c(1, 2, 1, 2, 1)
#' conditions <- c(1, 1, 2, 2, 3)
#'
#' # Resolution to select
#' binSize <- 500000
#'
#' # Instantiation of data set
#' hic.experiment <- parseHiC(
#'     paths,
#'     replicates = replicates,
#'     conditions = conditions,
#'     binSize = binSize
#' )
#' }
#' @usage
#' parseHiC(paths, replicates, conditions, binSize)
#'
#' @importFrom pbapply pbmapply
#' @export
parseHiC <- function(paths, binSize, replicates, conditions) {
    if (is.factor(paths)) {
        paths <- as.vector(paths)
    }
    if (!is.character(paths)) {
        stop("'paths' must be a vector of characters.", call. = FALSE)
    }
    for (path in paths) {
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
    if (length(replicates) != length(paths)) {
        stop("'replicates' should have the same length as 'paths'")
    }

    if (is.factor(conditions)) {
        conditions <- as.vector(conditions)
    }
    if (is.null(conditions)) {
        stop("'conditions' must be a vector of conditions.", call. = FALSE)
    }
    if (length(conditions) != length(paths)) {
        stop("'conditions' should have the same length as 'paths'")
    }

    if (!is.numeric(binSize) || length(binSize) != 1) {
        stop("'binSize' must be an integer.", call. = FALSE)
    }
    binSize <- as.integer(binSize)

    interactionSet <- pbapply::pbmapply(
        .parseOneHiC,
        path      = paths,
        binSize   = binSize,
        condition = conditions,
        replicate = replicates
    )

    mergedinteractionSet <- Reduce(
        f = mergeInteractionSet,
        x = interactionSet
    )
    mergedinteractionSet <- .sortHiCData(mergedinteractionSet)

    return(mergedinteractionSet)
}
