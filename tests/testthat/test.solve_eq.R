context("Solve equations")

test_that("find_param works when either signal strength or intercept is 0", {
  expect_equal(length(find_param(kappa = 0.1, gamma = 0, intercept = FALSE)), 3)
  expect_equal(length(find_param(kappa = 0.1, gamma = 0, intercept = TRUE)), 4)
  expect_equal(length(find_param(kappa = 0.1, gamma = sqrt(5), intercept = FALSE)), 3)
  expect_equal(length(find_param(kappa = 0.1, gamma = sqrt(5), intercept = TRUE)), 4)
})


