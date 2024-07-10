#' @description
#' Parses a single interactions file in \code{.cool} or \code{.mcool} format.
#'
#' @param path
#' The path to the interactions file.
#' @param binSize
#' The resolution (span of each position in number of bases). Optionally
#' provided to select the appropriate resolution in \code{.mcool} files.
#' Defaults to NULL.
#'
#' @return
#' A data.table of interactions.
#'
#' @keywords internal
#' @noRd
#' @importFrom data.table setorder
#' @importFrom GenomicRanges GRanges
.parseOneCool <- function(path, binSize = NA, replicate, condition) {
    if(!requireNamespace('rhdf5')) stop("'rhdf5' package is required. Please install it and retry.")
    message("\nParsing '", path, "'.")
    # TODO : ici s'il y a un binSize ça marche pas sur les données d'exemples
    uri <- function(path) {
        if (!is.numeric(binSize)) return(path)
        return(
            paste(
                "resolutions",
                format(binSize, scientific = FALSE),
                path,
                sep = "/"
            )
        )
    }

    bins <- data.table::data.table(
        chromosome = factor(
            rhdf5::h5read(file = path, name = uri("bins/chrom"))
        ),
        start = rhdf5::h5read(file = path, name = uri("bins/start")),
        end = rhdf5::h5read(file = path, name = uri("bins/end"))
    )
    bins[, start := as.integer(start)]
    bins[, start := start+1]
    bins[, end := as.integer(end)]

    setorder(bins, chromosome, start, end)
    bins[, index := seq_len(nrow(bins))]

    interactions <- data.table::data.table(
        id1 = rhdf5::h5read(file = path, name = uri("pixels/bin1_id")),
        id2 = rhdf5::h5read(file = path, name = uri("pixels/bin2_id")),
        interaction = rhdf5::h5read(file = path, name = uri("pixels/count"))
    )
    interactions[, id1 := as.integer(id1) + 1]
    interactions[, id2 := as.integer(id2) + 1]
    interactions[, interaction := as.numeric(interaction)]

    order1 <- match(interactions$id1, bins$index)
    order2 <- match(interactions$id2, bins$index)
    allRegions <- GenomicRanges::GRanges(bins)

    # GInteractions part
    gi <- InteractionSet::GInteractions(
        allRegions[order1],
        allRegions[order2],
        regions = allRegions,
        mode="strict"
    )
    assay <- as.matrix(interactions$interaction, ncol = 1)

    interactionSet <- .createInteractionSet(assay, gi, allRegions, condition, replicate)
    return(interactionSet)
}

#' @title Parser for data in cool format
#' @description
#' Parses interactions in \code{.cool} or \code{.mcool} format and returns
#' an InteractionSet object.
#' @param paths
#' A vector of paths to \code{.cool} or \code{.mcool} files.
#' @param replicates
#' A vector of replicate names repeated along the conditions.
#' @param conditions
#' A vector of condition names repeated along the replicates.
#' @param binSize
#' The resolution (span of each position in number of bases). Optionally
#' provided to select the appropriate resolution in \code{.mcool} files.
#' Defaults to NULL.
#'
#' @return
#' An InteractionSet.
#' @examples
#'   \dontrun{
#'     # Path to each file
#'     paths = c(
#'       'path/to/condition-1.replicate-1.cool',
#'       'path/to/condition-1.replicate-2.cool',
#'       'path/to/condition-2.replicate-1.cool',
#'       'path/to/condition-2.replicate-2.cool',
#'       'path/to/condition-3.replicate-1.cool'
#'     )
#'     # Replicate and condition of each file. Can be names instead of numbers.
#'     replicates <- c(1, 2, 1, 2, 1)
#'     conditions <- c(1, 1, 2, 2, 3)
#'     # Resolution to select in .mcool files
#'     binSize = 500000
#'     # Instantiation of data set
#'     object <- HiCDOCDataSetFromCool(
#'       paths,
#'       replicates = replicates,
#'       conditions = conditions,
#'       binSize = binSize # Specified for .mcool files.
#'     )
#'   }
#' @importFrom pbapply pbmapply
#' @export
parseCool <- function(paths, binSize=NA, replicates, conditions) {
    if(!requireNamespace('rhdf5')) stop("'rhdf5' package is required. Please install it and retry.")
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
        conditions <- as.vector(replicates)
    }
    if (is.null(replicates)) {
        stop("'replicates' must be a vector of replicates.", call. = FALSE)
    }

    if (is.factor(conditions)) {
        conditions <- as.vector(conditions)
    }
    if (is.null(conditions)) {
        stop("'conditions' must be a vector of conditions.", call. = FALSE)
    }

    if (!is.na(binSize) && (!is.numeric(binSize) || length(binSize) != 1)) {
        stop("'binSize' must be an integer.", call. = FALSE)
    }

    interactionSetCool <- pbapply::pbmapply(
        .parseOneCool,
        path = paths,
        binSize = binSize,
        condition = conditions,
        replicate = replicates
    )

    mergedinteractionSetCool <- Reduce(
        f = mergeInteractionSet,
        x = interactionSetCool
    )
    mergedinteractionSetCool <- .sortHiCData(mergedinteractionSetCool)

    return(mergedinteractionSetCool)
}
