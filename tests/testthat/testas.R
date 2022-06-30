library(geonet)
library(spatstat)

test_that("test",{
  X <- as_gn(simplenet)
  expect_s3_class(X, "gn")
})
