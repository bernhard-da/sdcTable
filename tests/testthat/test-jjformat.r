context("test jj_format")
data("microData1", package = "sdcTable")

# create hierarchies
dimList <- list(
  region = hier_create(
    root = "Total",
    nodes = LETTERS[1:4]
  ),
  gender = hier_create(
    root = "Total",
    nodes = c("male", "female")
  )
)

# third column containts a numeric variable ('val')
numVarInd <- 3

# creating an problem instance
prob <- makeProblem(
  data = microData1,
  dimList = dimList,
  numVarInd = numVarInd
)


# create inputs for jj format
inp <- createJJFormat(prob)

expect_is(inp, "jjformat")
expect_identical(length(inp), 5L)
expect_identical(inp[[1]], 0)
expect_identical(inp[[2]], 15L)

expect_is(inp[[3]], "data.table")
expect_identical(nrow(inp[[3]]), 15L)
expect_identical(ncol(inp[[3]]), 10L)

expect_identical(inp[[3]]$ind, as.numeric(0:14))
expect_identical(inp[[3]]$freqs[1], 100)
expect_identical(inp[[3]]$costs[1], 100)
expect_identical(inp[[3]]$val[1], 1284)
expect_identical(inp[[3]]$ubi[1], 150)

expect_identical(inp[[4]], 8)

expect_is(inp[[5]], "data.table")
expect_identical(nrow(inp[[5]]), 8L)
expect_identical(ncol(inp[[5]]), 4L)
expect_identical(inp[[5]]$v4[1], "0 ( -1 ) 3 ( 1 ) 6 ( 1 ) 9 ( 1 ) 12 ( 1 )")
