test_that("mergeInteractionSet works as expected", {
    path <- system.file("extdata", "liver_18_10M_500000.tsv", package = "HiCParser")
    object1 <- parseTabular(path, sep = "\t")
    # Reducing object1
    object1 <- object1[1:100]

    # Creating an 2nd fictive object, different from the first
    # Ginteractions
    set.seed(12345)
    allRegions <- regions(object1)
    gi2 <- InteractionSet::GInteractions(
        anchor1 = sample(allRegions$index, 300, replace=TRUE),
        anchor2 = sample(allRegions$index, 300, replace=TRUE),
        regions = allRegions,
        mode="strict"
    )
    gi2 <- unique(gi2)
    # Assays
    assays2 <- assays2[1:length(gi2),, drop=FALSE]
    # colData
    colData2[1,"condition"] <- 2
    # ISet object
    object2 <- InteractionSet::InteractionSet(assays2, gi2, colData=colData2)

    expect_error(
        objectMerged <- mergeInteractionSet(object1, object2),
        NA
    )
    expect_type(objectMerged, "S4")
    expect_true("InteractionSet" %in% class(objectMerged))
    expect_identical(objectMerged@colData,
                     DataFrame("condition"=c("1","2"),
                               "replicate"=c("R1", "R1")))
    expect_true(length(objectMerged) == 189)
    nbNA <- apply(objectMerged@assays@data[[1]], 2,function(x) sum(is.na(x)))
    expect_identical(nbNA, c(89L, 23L))
    expect_true("StrictGInteractions" %in% class(objectMerged@interactions))
    expect_identical(
        which(objectMerged@interactions@anchor1>objectMerged@interactions@anchor2),
        vector("integer", 0))
})
