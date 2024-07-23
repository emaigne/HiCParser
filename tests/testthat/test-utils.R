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
    assays2 <- matrix(sample(seq_len(2000), length(gi2)), ncol=1)
    # colData
    colData2 <- DataFrame("condition" = "2", "replicate" = "R2")
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
                               "replicate"=c("R1", "R2")))
    expect_true(length(objectMerged) == 189)
    nbNA <- apply(objectMerged@assays@data[[1]], 2,function(x) sum(is.na(x)))
    expect_identical(nbNA, c(89L, 23L))
    expect_true("StrictGInteractions" %in% class(objectMerged@interactions))
    expect_identical(
        which(objectMerged@interactions@anchor1>objectMerged@interactions@anchor2),
        vector("integer", 0))
})

test_that(".checkConditionsReplicates works as expected", {
    expect_error(.checkConditionsReplicates(c(1,2,3), c(1,2,3)),
                 NA)
    expect_error(.checkConditionsReplicates(c(1,2,3), NULL),
                 "'replicates' must be a vector of replicates.")
    expect_error(.checkReplicatesConditions(NULL, c(1,2,3)),
                 "'conditions' must be a vector of conditions")
    expect_error(.checkConditionsReplicates(c(1,NA,3), c(1,2,3)),
                 "'replicates' and 'conditions' can't containe NA values")
    expect_error(.checkConditionsReplicates(c(1,2,3), c(1,3)),
                 "'conditions' and 'replicates' must have the same length")

    expect_identical(
        .checkConditionsReplicates(factor(c(1,2,3)), factor(c(2,2,3))),
        list("conditions"=c("1","2","3"), "replicates"=c("2","2","3")))

})

test_that(".checkPaths works as expected", {
    path <- system.file("extdata", "liver_18_10M_500000.tsv", package = "HiCParser")
    tmp <- expect_error(.checkPaths(path=factor(path)),
                 NA)
    expect_identical(tmp, path)
    expect_error(.checkPaths("matrixPaths"=c(path, NA)),
                 "matrixPaths can't contain NA values.")
    expect_error(.checkPaths("matrixPaths"=c(2, 1)),
                 "matrixPaths must be a non empty vector of characters.")
    expect_error(.checkPaths("matrixPaths"="dont_exists"),
                 "'dont_exists' does not exist.")
})
