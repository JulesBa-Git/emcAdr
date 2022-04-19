library(emcAdr)
test_that("HelloWorld works", {
  a = 12
  expect_snapshot_output(HelloWorld(a), a)
  expect_equal(renvoieA(),a)
})
