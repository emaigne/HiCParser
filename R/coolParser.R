#' @description
#' Parses a single interactions file in \code{.cool} or \code{.mcool} format.
#'
#' @param path
#' The path to the interactions file.
#' @param binSize
#' The resolution (span of each position in number of bases). Optionally
#' provided to select the appropriate resolution in \code{.mcool} files.
#' Defaults to NULL.
#' @details
#' To read `.cool` of `.mcool` files, the `rhdf5` package is required.
#' Please install it before running the function.
#'
#' @return
#' A data.table of interactions.
#'
#' @keywords internal
#' @noRd
#' @importFrom data.table setorder
#' @importFrom GenomicRanges GRanges
.parseOneCool <- function(path, binSize = NA, condition, replicate) {
    if (!requireNamespace("rhdf5")) {
        stop("'rhdf5' package is required. Please install it and retry.")
    }
    message("\nParsing '", path, "'.")
    uri <- function(path) {
        if (!is.numeric(binSize)) {
            return(path)
        }
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
    bins[, start := start + 1]
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
        mode = "strict"
    )
    assay <- as.matrix(interactions$interaction, ncol = 1)

    interactionSet <-
        .createInteractionSet(assay, gi, allRegions, condition, replicate)
    return(interactionSet)
}

#' @title Parser for data in cool format
#' @description
#' Parses interactions in \code{.cool} or \code{.mcool} format and returns
#' an InteractionSet object.
#' @param paths
#' A vector of paths to \code{.cool} or \code{.mcool} files.
#' @param conditions
#' A vector of condition names repeated along the replicates.
#' @param replicates
#' A vector of replicate names repeated along the conditions.
#' @param binSize
#' The resolution (span of each position in number of bases). Optionally
#' provided to select the appropriate resolution in \code{.mcool} files.
#' Defaults to NULL.
#' @details
#' To read `.cool` of `.mcool` files, the `rhdf5` package is required.
#' Please install it before running the function.
#'
#' @return
#' An InteractionSet.
#' @examples
#' # EXAMPLE FOR .cool FORMAT
#' # Path to each file
#' pathsCool <- c(
#'     "path/to/condition-1.replicate-1.cool",
#'     "path/to/condition-1.replicate-2.cool",
#'     "path/to/condition-1.replicate-3.cool",
#'     "path/to/condition-2.replicate-1.cool",
#'     "path/to/condition-2.replicate-2.cool",
#'     "path/to/condition-2.replicate-3.cool"
#' )
#' # Condition and replicate of each file. Can be names instead of numbers.
#' conditions <- c(1, 1, 1, 2, 2, 2)
#' replicates <- c(1, 2, 3, 1, 2, 3)
#' if (FALSE) {
#'     library(rhdf5)
#'     object <- parseCool(
#'         paths,
#'         conditions = conditions,
#'         replicates = replicates
#'     )
#' }
#'
#' # EXAMPLE FOR .mcool FORMAT
#' # Resolution
#' binSize <- 500000
#' # Path to each file
#' paths <- c(
#'     "path/to/condition-1.replicate-1.mcool",
#'     "path/to/condition-1.replicate-2.mcool",
#'     "path/to/condition-1.replicate-3.mcool",
#'     "path/to/condition-2.replicate-1.mcool",
#'     "path/to/condition-2.replicate-2.mcool",
#'     "path/to/condition-2.replicate-3.mcool"
#' )
#' # Condition and replicate of each file. Can be names instead of numbers.
#' conditions <- c(1, 1, 1, 2, 2, 2)
#' replicates <- c(1, 2, 3, 1, 2, 3)
#' if (FALSE) {
#'     # Instantiation of data set
#'     library(rhdf5)
#'     object <- parseCool(
#'         paths,
#'         conditions = conditions,
#'         replicates = replicates,
#'         binSize = binSize
#'     )
#' }
#' @importFrom pbapply pbmapply
#' @export
parseCool <- function(paths, binSize = NA, conditions, replicates) {
    if (!requireNamespace("rhdf5")) {
        stop("'rhdf5' package is required. Please install it and retry.")
    }
    paths <- .checkPaths("paths" = paths)

    repCond <- .checkConditionsReplicates(conditions, replicates)
    if (min(lengths(repCond)) != length(paths)) {
        stop(
            "'conditions/replicates' and 'paths' ",
            "must have the same length",
            call. = FALSE
        )
    }

    if (!is.na(binSize) && all(grepl(".cool", paths, fixed = TRUE))) {
        warning("binSize specified but all files are not .cool, ignored")
        binSize <- NA
    }
    if (is.na(binSize) && any(grepl("mcool", paths, fixed = TRUE))) {
        stop("binSize must be specified fot mcool files")
    }
    if (!is.na(binSize) && (!is.numeric(binSize) || length(binSize) != 1)) {
        stop("'binSize' must be an integer.", call. = FALSE)
    }

    interactionSetCool <- pbapply::pbmapply(
        .parseOneCool,
        path = paths,
        binSize = binSize,
        condition = repCond[["conditions"]],
        replicate = repCond[["replicates"]]
    )

    mergedinteractionSetCool <- Reduce(
        f = mergeInteractionSet,
        x = interactionSetCool
    )
    mergedinteractionSetCool <- .sortHiCData(mergedinteractionSetCool)

    return(mergedinteractionSetCool)
}
