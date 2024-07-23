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
#' @examples
#' path <- system.file("extdata", "liver_18_10M_500000.tsv", package = "HiCParser")
#' object1 <- parseTabular(path, sep = "\t")
#' # Creating an object with a different condition
#' object2 <- parseTabular(path, sep = "\t")
#' object2@colData[1,"condition"] <- 2
#' objectMerged <- mergeInteractionSet(object1, object2)
#'
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


#' Check input paths for parser functions (format and existence)
#'
#' @param ... a vector of paths, named, containing input data
#'
#' @return a  caracter vector of paths
#' @noRd
#'
#' @examples
#' path <- system.file("extdata", "liver_18_10M_500000.tsv", package = "HiCParser")
#' .checkPaths("ExamplePaths" = path)
.checkPaths <- function(...){
    args <- list(...)
    path <- args[[1]]
    namePath <- names(vapply(match.call(), deparse, FUN.VALUE = "character"))[-1]
    if (is.factor(path)) {
        path <- as.vector(path)
    }
    if (anyNA(path)) {
        stop(namePath," can't contain NA values.", call. = FALSE)
    }
    if (!is.character(path)) {
        stop(namePath," must be a non empty vector of characters.", call. = FALSE)
    }
    for (p in path) {
        if (!file.exists(p)) {
            stop("'", p, "' does not exist.", call. = FALSE)
        }
    }
    return(path)
}

#' Check conditions and replicate arguments of parsers
#'
#' @param ... a vector of replicates or conditions
#'
#' @return a  caracter vector of paths
#' @noRd
#'
#' @examples
#' .checkConditionsReplicates(c(1, 1, 1, 2, 2, 2), c(1, 2, 3, 1, 2, 3))
.checkConditionsReplicates <- function(conditions, replicates){
    if (is.factor(replicates)) {
        replicates <- as.vector(replicates)
    }
    if (is.factor(conditions)) {
        conditions <- as.vector(conditions)
    }

    if (is.null(replicates)) {
        stop("'replicates' must be a vector of replicates.", call. = FALSE)
    }
    if (is.null(conditions)) {
        stop("'conditions' must be a vector of conditions.", call. = FALSE)
    }
    if (anyNA(c(replicates, conditions))) {
        stop("'replicates' and 'conditions' can't contain NA values", call. = FALSE)
    }
    if (length(conditions) != length(replicates)) {
        stop(
            "'conditions' and 'replicates' must have the same length",
            call. = FALSE
        )
    }
    return(list("conditions"=conditions, "replicates"=replicates))
}
