#' Create the object interactionSet to use in HiCDOCDataSet
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
