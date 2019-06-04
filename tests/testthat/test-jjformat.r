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

# creating an problem instance
prob <- makeProblem(
  data = microData1,
  dimList = dimList,
  numVarInd = "val"
)

# check errors
expect_error(createJJFormat(x = 5))

# create inputs for jj format
inp <- createJJFormat(prob)
expect_identical(digest::digest(inp), "8d0372b82939b2c571195d488fa70f85")

# no numvar
prob <- makeProblem(
  data = microData1,
  dimList = dimList
)
inp <- createJJFormat(prob)
expect_identical(digest::digest(inp), "4ac42d469bf1d2f8321906af1bf617a1")
