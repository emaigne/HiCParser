##############################################################################
# R code to generate inst/extdata files
options(scipen=999)

interactions_vector <-
    c(
        2096, 1341, 213, 191, 142, 127, 217, 175, 157, 113, 45, 72, 73,
        86, 94, 114, 89, 91, 69, 54,
        2592, 506, 307, 231, 197, 292, 215, 160, 167, 101, 93, 115, 131,
        108, 132, 111, 103, 76, 79,
        978, 783, 310, 285, 289, 184, 131, 108, 63, 64, 74, 99, 75, 104,
        79, 69, 48, 48,
        1925, 1167, 1034, 579, 178, 143, 202, 141, 120, 119, 141, 106, 107,
        114, 89, 71, 82,
        1889, 1016, 570, 241, 170, 169, 109, 59, 102, 110, 115, 91, 118,
        85, 67, 49,
        2036, 1284, 213, 125, 227, 161, 162, 143, 161, 102, 142, 117, 70,
        74, 57,
        2662, 1138, 518, 464, 210, 187, 243, 244, 234, 234, 182, 141, 91, 83,
        2673, 1115, 408, 159, 138, 200, 227, 203, 203, 158, 126, 69, 43,
        2261, 930, 276, 184, 180, 188, 186, 191, 159, 138, 66, 56,
        2737, 2073, 721, 361, 251, 185, 218, 200, 125, 131, 112,
        3111, 1840, 287, 196, 139, 163, 168, 99, 106, 97,
        2681, 858, 344, 206, 235, 195, 143, 133, 134,
        1662, 937, 427, 336, 251, 141, 117, 94,
        2538, 1707, 565, 373, 212, 144, 129,
        2171, 939, 479, 266, 149, 113,
        2900, 1885, 428, 293, 215,
        2956, 1114, 546, 317,
        2030, 1349, 430,
        2274, 1505,
        2817
    )

################################################
# liver_18_10M_500000.tsv
tsv_data <- data.frame(
    chromosome = rep(18, 210),
    position_1 = as.character(
        unlist(
            mapply(
                FUN = function(x, y) rep(x, y),
                seq(0, 9500000, by=500000),
                20:1
            )
        )
    ),
    position_2 = as.character(
        unlist(
            lapply(seq(0, 9500000, by = 500000),
                   function(x) seq(x, 9500000, by = 500000)
            )
        )
    ),
    `1_R1` = interactions_vector
)
colnames(tsv_data) <- c("chromosome", "position 1", "position 2", "1.R1")
write.table(tsv_data, "inst/extdata/liver_18_10M_500000.tsv",
            sep="\t", quote = FALSE, row.names = FALSE)


################################################
# liver_18_10M_500000.matrix and liver_18_10M_500000.bed

mat_data <- data.frame(
    index1 = as.character(
        unlist(
            mapply(FUN = function(x, y) rep(x, y), seq(4427, 4446), 20:1)
        )
    ),
    index2 = as.character(
        unlist(lapply(seq(4427, 4446), function(x) seq(x, 4446)))
    ),
    `1_R1` = interactions_vector
)

write.table(mat_data, "inst/extdata/liver_18_10M_500000.matrix",
            sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

################################################
# liver_18_10M_500000.bed

bed_data <- data.frame(
    chromosome = 18,
    start = as.character(seq(0, 9500000, by=500000)),
    end = as.character(seq(500000, 10000000, by=500000)),
    index = 4427:4446
)

write.table(bed_data, "inst/extdata/liver_18_10M_500000_2.bed",
            sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

################################################
# liver_18_10M_500000.bed
