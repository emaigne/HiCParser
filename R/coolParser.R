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
.parseOneCool <- function(path, binSize = NA, replicate, condition) {

    message("\nParsing '", path, "'.")

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

#' @description
#' Parses interactions in \code{.cool} or \code{.mcool} format and fills the
#' interactions slot of the provided \code{\link{HiCDOCDataSet}}.
#'
#' @param object
#' A \code{\link{HiCDOCDataSet}}.
#' @param binSize
#' The resolution (span of each position in number of bases). Optionally
#' provided to select the appropriate resolution in \code{.mcool} files.
#' Defaults to NULL.
#'
#' @return
#' A filled \code{\link{HiCDOCDataSet}}.
#'
#' @keywords internal
#' @noRd
.parseCool <- function(object, binSize = NA, replicates, conditions) {
    interactionSetCool <- pbapply::pbmapply(
        .parseOneCool,
        path = object@input,
        binSize = binSize,
        condition = conditions,
        replicate = replicates
    )

    mergedinteractionSetCool <- Reduce(
        f = .mergeInteractionSet,
        x = interactionSetCool
    )

    new(
        "HiCDOCDataSet",
        mergedinteractionSetCool,
        input = object@input
    )
}

