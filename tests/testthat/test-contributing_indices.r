context("test contributing_indices")
data("microData1", package = "sdcTable")

# specify hierarchies for `age` and `region`
dim_region <- hier_create(root = "Total", nodes = LETTERS[1:4])
dim_gender <- hier_create(root = "Total", nodes = c("male", "female"))
dl <- list(region = dim_region, gender = dim_gender)
prob <- makeProblem(
  data = microData1,
  dimList = dl
)

df <- sdcProb2df(prob, dimCodes = "original")

# errors
expect_error(contributing_indices(prob = prob, ids = 5))
expect_error(contributing_indices(prob = prob, ids = "01"))

# single values
ids <- contributing_indices(prob = prob, ids = "0101")

rawData <- slot(get.sdcProblem(prob, "dataObj"), "rawData")

expect_is(ids, "list")
expect_identical(names(ids), "0101")
expect_identical(ids[[1]], 1:2)

dt <- rawData[ids[[1]]]
expect_equal(unique(dt$region), "A")
expect_equal(unique(dt$gender), "female")

# all cells
ids <- contributing_indices(prob = prob)
expect_identical(digest::digest(ids), "bc9eb3d9155a36108221c0bb42904ae9")


for (i in 1:nrow(df)) {
  code_region <- df$region[i]
  if (code_region == "Total") {
    code_region <- LETTERS[1:4]
  }
  code_gender <- df$gender[i]
  if (code_gender == "Total") {
    code_gender <- c("male", "female")
  }
  ids <- contributing_indices(prob = prob, ids = df$strID[i])[[1]]
  dt <- rawData[ids]
  expect_true(all(unique(dt$region) %in% code_region))
  expect_true(all(unique(dt$gender) %in% code_gender))
}
