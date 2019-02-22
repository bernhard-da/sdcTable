context("test protectTable()")

sp <- searchpaths()
fn <- file.path(sp[grep("sdcTable", sp)], "data", "problemWithSupps.RData")

problem <- get(load(file = fn))

p1 <- protectTable(
  problem,
  method = "OPT",
  useC = FALSE,
  verbose = FALSE
)
p2 <- protectTable(
  problem,
  method = "OPT",
  useC = TRUE,
  verbose = FALSE
)
expect_equivalent(p1@finalData, p2@finalData)
expect_is(p1, "safeObj")
expect_equal(p1@nrPublishableCells, 11)
expect_equal(p1@nrSecondSupps, 3)
expect_equal(which(p1@finalData$sdcStatus != "s"), c(5, 6, 11, 12))

expect_is(p2, "safeObj")
expect_equal(p2@nrPublishableCells, 11)
expect_equal(p2@nrSecondSupps, 3)
expect_equal(which(p2@finalData$sdcStatus != "s"), c(5, 6, 11, 12))

p3 <- protectTable(problem, method = "HYPERCUBE", verbose = FALSE)
expect_equal(which(p3@finalData$sdcStatus != "s"), c(5, 6, 11, 12))

p4 <- protectTable(
  problem,
  method = "HITAS",
  useC = TRUE,
  verbose = FALSE
)
expect_equal(which(p4@finalData$sdcStatus != "s"), c(5, 6, 11, 12))

rm(problem)
fn <- file.path(sp[grep("sdcTable", sp)], "data", "problem.RData")
problem <- get(load(file = fn))

p5 <- protectTable(
  problem,
  method = "OPT",
  useC = TRUE,
  verbose = FALSE
)
expect_is(p5, "safeObj")
expect_equal(p5@nrPublishableCells, 15)
expect_equal(p5@nrSecondSupps, 0)
expect_equal(sum(p5@finalData$sdcStatus != "s"), 0)
