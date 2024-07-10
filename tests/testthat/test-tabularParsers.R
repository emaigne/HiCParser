test_that("tabular parser for one file produce correct format", {
    path <- system.file(
        "extdata",
        "liver_18_10M_500000.tsv",
        package = "HiCParser"
    )
    expect_message(object <- parseTabular(path), "liver_18_10M_500000.tsv")
    # Class and slots
    expect_s4_class(object, "InteractionSet")

    # TODO
    # Here test order of interactions and assays.
    expect_equal(length(object), 210)

    # Interactions
    expect_true("matrix" %in% class(SummarizedExperiment::assay(object)))
    expect_s4_class(InteractionSet::regions(object), "GRanges")
    expect_s4_class(InteractionSet::interactions(object), "StrictGInteractions")
    expect_s4_class(S4Vectors::mcols(object), "DataFrame")
    expect_true(is.numeric(SummarizedExperiment::assay(object)))
})
