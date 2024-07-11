test_that("Parser for .hic data works as expected", {
    paths <-
        system.file("extdata", "liver_18_10M.hic", package = "HiCParser")

    # Replicate and condition of each file. Can be names instead of numbers.
    replicates <- 1
    conditions <- 1
    binSize <- 500000

    # Instantiation of data set
    expect_message(
        object <- parseHiC(
            paths,
            binSize = binSize,
            replicates = replicates,
            conditions = conditions
        ),
        "liver_18_10M.hic"
    )
    expect_equal(length(object), 210)

    # Interactions
    expect_true("matrix" %in% class(SummarizedExperiment::assay(object)))
    expect_s4_class(InteractionSet::regions(object), "GRanges")
    expect_s4_class(InteractionSet::interactions(object), "StrictGInteractions")
    expect_s4_class(S4Vectors::mcols(object), "DataFrame")
    expect_true(is.numeric(SummarizedExperiment::assay(object)))
})
