#' @title Initialize an InteractionSet object from tabular data
#' @description
#' Parses interactions in tabular format and fills the conditions, replicates,
#' and interactions slots of the provided \code{\link{InteractionSet}}.
#'
#' @param tabular
#' a path to a tabulated data file
#' @param input
#' The name of the input file
#' @param conditions
#' The names of the conditions
#' (infered from the table if not specified).
#' @param replicates
#' The names of the replicates
#' (infered from the table if not specified).
#'
#' @return
#' An \code{\link{InteractionSet}}.
#'
#' @keywords internal
#' @noRd
#' @importFrom data.table setorder setnames shift merge.data.table
#' @importFrom gtools mixedsort
#' @importFrom GenomicRanges GRanges
.setFromTabular <- function(tabular, conditions = NULL, replicates = NULL) {
    if (colnames(tabular)[1] != "chromosome") {
        stop(
            "First column of the input file must be named 'chromosome'.",
            call. = FALSE
        )
    }
    if (colnames(tabular)[2] == "position.1") {
        data.table::setnames(tabular, "position.1", "position 1")
    }
    if (colnames(tabular)[2] != "position 1") {
        stop(
            "Second column of the input file must be named 'position 1'.",
            call. = FALSE
        )
    }
    if (colnames(tabular)[3] == "position.2") {
        data.table::setnames(tabular, "position.2", "position 2")
    }
    if (colnames(tabular)[3] != "position 2") {
        stop(
            "Third column of the input file must be named 'position 2'.",
            call. = FALSE
        )
    }
    if (is.null(conditions) != is.null(replicates)) {
        stop(
            "Conditions and replicates should be both NULL, or none.",
            call. = FALSE
        )
    }
    tabular[, chromosome := as.character(chromosome)]
    tabular[, chromosome := factor(
        chromosome,
        levels = gtools::mixedsort(unique(chromosome))
    )]
    setorder(tabular, chromosome, `position 1`, `position 2`)
    # Assays part, fill with NA
    assays <- as.matrix(tabular[, 4:ncol(tabular), drop = FALSE])

    if (!is.null(conditions) | !is.null(replicates)) {
        if (
            (length(conditions) != ncol(assays)) |
                (length(replicates) != ncol(assays))
        ) {
            stop(
                "Number of conditions and replicates should match the number ",
                "of counts in the matrix.",
                call. = FALSE
            )
        }
    } else {
        if (!all(grepl("^.+?\\..+$", colnames(assays)))) {
            stop(
                "Fourth to last column of the input file must be named 'C.R', ",
                "with C the condition number/name and R the replicate ",
                "number/name.",
                call. = FALSE
            )
        }
    }
    assays[assays == 0] <- NA

    # GInteraction part
    tabular <- tabular[, .(chromosome, `position 1`, `position 2`)]
    data.table::setnames(tabular, "position 1", "bin.1")
    data.table::setnames(tabular, "position 2", "bin.2")

    diagonal <- (tabular$bin.1 == tabular$bin.2)
    binSize <- .modeVector(abs(
        tabular[!diagonal, ]$bin.1 - tabular[!diagonal, ]$bin.2
    ))

    tabular[, bin.1 := bin.1 / binSize]
    tabular[, bin.2 := bin.2 / binSize]

    allRegions <- data.table::melt(
        tabular[, .(chromosome, bin.1, bin.2)],
        id.vars = "chromosome",
        value.name = "indexC"
    )
    allRegions[, variable := NULL]
    allRegions <- unique(allRegions)
    setorder(allRegions, chromosome, indexC)

    # Constructing unique index for all chromosomes,
    # taking into account the difference in bins.
    allRegions[
        ,
        index := indexC - data.table::shift(indexC, fill = 0),
        by = .(chromosome)
    ]
    allRegions[index == 0, index := 1]
    allRegions[, index := cumsum(index)]
    allRegions[, end := (indexC + 1) * binSize]
    allRegions[, start := (indexC) * binSize + 1]
    data.table::setcolorder(
        allRegions,
        c("chromosome", "start", "end", "index", "indexC")
    )

    tabular <- data.table::merge.data.table(
        tabular,
        allRegions[, .(chromosome, startIndex = index, bin.1 = indexC)],
        all.x = TRUE,
        sort = FALSE,
        by = c("chromosome", "bin.1")
    )
    tabular <- data.table::merge.data.table(
        tabular,
        allRegions[, .(chromosome, stopIndex = index, bin.2 = indexC)],
        all.x = TRUE,
        sort = FALSE,
        by = c("chromosome", "bin.2")
    )
    tabular[, bin.1 := NULL]
    tabular[, bin.2 := NULL]

    allRegions[, indexC := NULL]
    order1 <- match(tabular$startIndex, allRegions$index)
    order2 <- match(tabular$stopIndex, allRegions$index)
    allRegions <- GenomicRanges::GRanges(allRegions)

    gi <- InteractionSet::GInteractions(
        allRegions[order1],
        allRegions[order2],
        regions = allRegions,
        mode = "strict"
    )

    if (is.null(conditions)) {
        conditions <- gsub("^(.+?)\\..+$", "\\1", colnames(assays))
    }
    if (is.null(replicates)) {
        replicates <- gsub("^.+?\\.(.+)$", "\\1", colnames(assays))
    }

    interactionSet <-
        .createInteractionSet(assays, gi, allRegions, conditions, replicates)
    return(interactionSet)
}

#' @title Parser for tabular data
#' @description
#' Read the file, produce an \code{InteractionSet} object.
#'
#' @param path
#' A path to a tabular file.
#' @param sep
#' The separator of the tabular file. Default to tabulation.
#'
#' @details
#' Accepts a tabular file with \code{chromosome}, \code{position 1},
#' \code{position 2}, and multiple replicate columns listing interaction counts.
#' Null interactions do not have to be listed. Values must be separated by
#' tabulations. The header must be
#' \code{chromosome position 1 position 2 x.y x.y x.y ...} with \code{x}
#' replaced by condition names and \code{y} replaced by replicate names.
#'
#' @return
#' An InteractionSet object.
#'
#' @examples
#' path <- system.file("extdata", "liver_18_10M_500000.tsv",
#'     package = "HiCParser"
#' )
#' object <- parseTabular(path, sep = "\t")
#'
#' @usage
#' parseTabular(path, sep = '\t')
#'
#' @importFrom data.table fread
#' @export
parseTabular <- function(path, sep = "\t") {
    path <- .checkPaths("path" = path)
    if (length(path) > 1) {
        stop("'path' must be 1-lengths character vector.", call. = FALSE)
    }

    message("Parsing '", path, "'.")

    interactions <- data.table::fread(
        file = path,
        sep = sep,
        header = TRUE,
        # comment.char = "#",
        check.names = FALSE,
        data.table = TRUE,
        stringsAsFactors = FALSE
    )

    interactionSet <- .setFromTabular(interactions)
    interactionSet <- .sortHiCData(interactionSet)

    return(interactionSet)
}
