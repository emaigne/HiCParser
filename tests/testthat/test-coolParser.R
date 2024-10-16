test_that("Parser for .cool data works as expected", {
    paths <-
        system.file("extdata", "liver_18_10M_500000.cool",
            package = "HiCParser"
        )

    # Replicate and condition of each file. Can be names instead of numbers.
    replicates <- "GG"
    conditions <- "AA"
    binSize <- 500000

    # Instantiation of data set
    expect_message(
        object <- parseCool(
            paths,
            replicates = replicates,
            conditions = conditions
        ),
        "liver_18_10M_500000.cool"
    )
    expect_equal(length(object), 210)

    # Interactions
    expect_true("matrix" %in% class(SummarizedExperiment::assay(object)))
    expect_s4_class(InteractionSet::regions(object), "GRanges")
    expect_s4_class(
        InteractionSet::interactions(object),
        "StrictGInteractions"
    )
    expect_s4_class(S4Vectors::mcols(object), "DataFrame")
    expect_true(is.numeric(SummarizedExperiment::assay(object)))
})

test_that("Parser for .cool data returns expected errors", {
    paths <-
        system.file("extdata", "liver_18_10M_500000.cool",
            package = "HiCParser"
        )
    expect_error(
        parseCool(
            paths,
            replicates = c(1, 2),
            conditions = c(1, 2)
        ),
        "'conditions/replicates' and 'paths' must have the same length"
    )
})
