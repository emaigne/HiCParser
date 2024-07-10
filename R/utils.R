#' @title Merge two InteractionSet objects
#' @description
#' Merge two different \code{\link{InteractionSet}}.
#'
#' @param interactionSet1
#' The first \code{\link{InteractionSet}}.
#' @param interactionSet2
#' The second \code{\link{InteractionSet}}.
#' @param fill
#' Fill missing values with this.
#'
#' @return
#' The merged \code{\link{InteractionSet}}.
#'
#' @importFrom GenomicRanges union
#' @importFrom InteractionSet interactions
#' @importFrom BiocGenerics cbind
#' @export
mergeInteractionSet <- function(interactionSet1, interactionSet2, fill = NA) {
    unionInteractions <- GenomicRanges::union(
        InteractionSet::interactions(interactionSet1),
        InteractionSet::interactions(interactionSet2)
    )
    # Complete InteractionSets
    interactionSet1 <- .fillInteractionSet(
        interactionSet1,
        unionInteractions,
        fill
    )
    interactionSet2 <- .fillInteractionSet(
        interactionSet2,
        unionInteractions,
        fill
    )

    # Merge
    newiset <- BiocGenerics::cbind(interactionSet1, interactionSet2)
    return(newiset)
}
