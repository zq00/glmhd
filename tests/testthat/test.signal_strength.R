context("Signal Strength")

test_that("signal_strength works when either signal strength or intercept is 0", {
  expect_equal(length(signal_strength(kappa_hat = 0.5, intercept = FALSE)), 2)
})
