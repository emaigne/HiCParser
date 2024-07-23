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
    paths <-
        system.file("extdata", "liver_18_10M.hic", package = "HiCParser")
    binSize <- 500000
    expect_error(object <- parseHiC(rep(paths, 3),
                                    binSize=binSize,
                                    conditions=c(1,2,3),
                                    replicates=NULL),
                 "'replicates' must be a vector of replicates.")
    expect_error(object <- parseHiC(rep(paths, 3),
                                    binSize=binSize,
                                    conditions=NULL,
                                    replicates=c(1,2,3)),
                 "'conditions' must be a vector of conditions")
    expect_error(object <- parseHiC(rep(paths, 3),
                                    binSize=binSize,
                                    conditions=c(1,NA,3),
                                    replicates=c(1,2,3)),
                 "'replicates' and 'conditions' can't contain NA values")
    expect_error(object <- parseHiC(rep(paths, 3),
                                    binSize=binSize,
                                    conditions=c(1,2,3),
                                    replicates=c(1,3)),
                 "'conditions' and 'replicates' must have the same length")

    expect_error(object <- parseHiC(rep(paths, 3),
                                    binSize=binSize,
                                    conditions=factor(c(1,2,3)),
                                    replicates=factor(c(1,2,3))),
                 NA)
    expect_identical(SummarizedExperiment::colData(object),
                     S4Vectors::DataFrame("condition"=c("1", "2", "3"),
                                          "replicate"=c("1", "2", "3")))
})

test_that(".checkPaths works as expected", {
    paths <-
        system.file("extdata", "liver_18_10M.hic", package = "HiCParser")
    binSize <- 500000

    expect_error(object <- parseHiC("paths"=c(paths, NA),
                                    binSize=binSize,
                                    conditions=c(1,2),
                                    replicates=c(1,2)),
                 "paths can't contain NA values.")
    expect_error(object <- parseHiC("paths"=c(2, 1),
                                    binSize=binSize,
                                    conditions=c(1,2),
                                    replicates=c(1,2)),
                 "paths must be a non empty vector of characters.")
    expect_error(object <- parseHiC("paths"=c(2, 1),
                                    binSize=binSize,
                                    conditions=c(1,2),
                                    replicates=c(1,2)),
                 "paths must be a non empty vector of characters.")
    expect_error(object <- parseHiC("paths"=c("dont_exists"),
                                    binSize=binSize,
                                    conditions=c(1),
                                    replicates=c(1)),
                 "'dont_exists' does not exist.")
})
