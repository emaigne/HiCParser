#' .modeVector Extract the mode of vector.
#'
#' @param x A vector
#'
#' @return The mode of a vector
#'
#' @examples
#' .modeVector(c(1, 2, 2, 2, 4))
#' @noRd
.modeVector <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
}

#' @description
#' Fills parameters and slots describing the data. Called by a
#' \code{\link{HiCParser}} constructor.
#'
#' @param object
#' A \code{InteractionSet}.
#'
#' @return
#' A sorted \code{InteractionSet}.
#'
#' @keywords internal
#' @noRd
.sortHiCData <- function(object) {
    # Reduce the levels in interaction part
    object <- InteractionSet::reduceRegions(object)
    objectRegions <- InteractionSet::regions(object)
    chromosomeNames <- unique(as.character(
        GenomeInfoDb::seqnames(objectRegions)
    ))
    chromosomeNames <- gtools::mixedsort(chromosomeNames)
    GenomeInfoDb::seqlevels(
        InteractionSet::regions(object),
        pruning.mode = "coarse"
    ) <- chromosomeNames

    # Add chromosome column for split purpose
    chromosomes <-
        GenomeInfoDb::seqnames(InteractionSet::anchors(object, "first"))
    chromosomes <- S4Vectors::Rle(factor(chromosomes, levels = chromosomeNames))
    S4Vectors::mcols(object) <- S4Vectors::DataFrame("chromosome" = chromosomes)

    # Sorting interactions and assay
    ids <- InteractionSet::anchors(object, id = TRUE)
    neworder <- order(chromosomes, ids$first, ids$second)
    object <- object[neworder, ]
    return(object)
}

#' Create the object InteractionSet from given conditions and replicates
#'
#' @param assay matrix of assays
#' @param gi GInteractions object
#' @param allRegions regions object
#' @param condition condition (for colData), vector of length 1 or more
#' @param replicate replicate (for colData), vector of length 1 or more
#'
#' @return an interactionSet object
#' @keywords internal
#' @noRd
#' @importFrom S4Vectors DataFrame
.createInteractionSet <- function(assay, gi, allRegions, condition, replicate){
    # Add 0 on the diagonal if there is only off diagonal interaction
    ids <- InteractionSet::anchors(gi, id = TRUE)
    idsDiagonals <- ids$first[ids$first == ids$second]
    notPresent <- setdiff(unique(c(ids$first, ids$second)), idsDiagonals)
    nb <- length(notPresent)
    colnames(assay) <- NULL
    if(nb>0){
        notPresentRegions <- allRegions[match(notPresent, allRegions$index)]
        gi <- c(gi,
                InteractionSet::GInteractions(
                    notPresentRegions,
                    notPresentRegions,
                    regions = allRegions,
                    mode="strict"))
        assay <- rbind(assay, as.matrix(rep(0, nb*ncol(assay)),
                                        ncol=ncol(assay),
                                        nrow=nb))
    }

    interactionSet <- InteractionSet::InteractionSet(
        assays = assay,
        interactions = gi,
        colData=S4Vectors::DataFrame(
            "condition" = condition,
            "replicate" = replicate
        )

    )

    # Keep only intra-chromosomal interactions
    interactionSet <- interactionSet[InteractionSet::intrachr(interactionSet),]
    return(interactionSet)
}

#' @description
#' Fill \code{\link{InteractionSet}} with possibly missing values
#'
#' @param interactionSet
#' An \code{\link{InteractionSet}}.
#' @param interactionSetUnion
#' The full \code{\link{InteractionSet}}.
#' @param fill
#' Fill missing values with this.
#'
#' @return
#' The full \code{\link{InteractionSet}}.
#'
#' @keywords internal
#' @noRd
#' @importFrom InteractionSet InteractionSet
#' @importFrom GenomicRanges match
#' @importFrom SummarizedExperiment assay colData
.fillInteractionSet <- function(
        interactionSet,
        interactionSetUnion,
        fill = NA
) {
    over <- GenomicRanges::match(interactionSet, interactionSetUnion)
    totalColumns <- ncol(interactionSet)
    newAssays <- matrix(
        rep(fill, length(interactionSetUnion) * totalColumns),
        ncol = totalColumns
    )
    newAssays[over, ] <- SummarizedExperiment::assay(interactionSet)
    return(
        InteractionSet::InteractionSet(
            newAssays,
            interactionSetUnion,
            colData = SummarizedExperiment::colData(interactionSet)
        )
    )
}
