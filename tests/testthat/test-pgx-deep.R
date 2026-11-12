skip_if_not_installed("torch")
skip_if(!torch::torch_is_installed())

make_sae <- function(...) {
  fx <- make_sae_fixture()
  net <- suppressMessages(MultiOmicsSAE$new(fx$X, fx$Y, ntop = 0, ...))
  list(net = net, fx = fx)
}

test_that("MultiOmicsSAE constructs on tiny synthetic data", {
  s <- make_sae()
  expect_s3_class(s$net, "MultiOmicsSAE")
  expect_equal(s$net$nviews, 3)
  expect_equal(names(s$net$X), c("gx", "px", "mx"))
  expect_equal(s$net$labels$class, paste0("class", 1:3))
})

test_that("fit() lowers the training cross-entropy", {
  s <- make_sae()
  capture.output(s$net$fit(niter = 10))
  h <- s$net$loss_history
  expect_length(h, 10)
  expect_true(all(is.finite(h)))
  expect_lt(h[length(h)], h[1])
})

test_that("get_dims() returns per-view encoder/decoder and predictor dims", {
  s <- make_sae()
  d <- s$net$get_dims()
  expect_named(d, c("encoder", "decoder", "integrator", "predictor"))
  expect_equal(names(d$encoder), c("gx", "px", "mx"))
  expect_equal(names(d$decoder), c("gx", "px", "mx"))
  expect_equal(names(d$predictor), "class")
  ## encoder starts at the input size; with GLU (default) the last predictor
  ## linear layer has 2 * K outputs, which the GLU halves back to K
  expect_true(all(sapply(d$encoder, head, 1) == 30))
  expect_equal(tail(d$predictor$class, 1), 2 * 3)
})

test_that("predict() returns an n x K matrix per target", {
  s <- make_sae()
  capture.output(s$net$fit(niter = 10))
  s$net$model$train(FALSE)
  p <- s$net$predict()
  expect_length(p, 1)
  expect_equal(dim(p[[1]]), c(24L, 3L))
  expect_equal(colnames(p[[1]]), paste0("class", 1:3))
  expect_equal(rownames(p[[1]]), colnames(s$fx$X[[1]]))
})

test_that("predict() is deterministic in eval mode even with input noise", {
  s <- make_sae(add_noise = 0.5)
  capture.output(s$net$fit(niter = 5))
  s$net$model$train(FALSE)
  expect_identical(s$net$predict(), s$net$predict())
})

test_that("input noise is active in train mode", {
  s <- make_sae(add_noise = 0.5)
  s$net$model$train(TRUE)
  xx <- lapply(s$net$X, function(x) torch::torch_tensor(scale(t(x))))
  a <- as.matrix(s$net$model(xx)$y_pred[[1]])
  b <- as.matrix(s$net$model(xx)$y_pred[[1]])
  expect_false(isTRUE(all.equal(a, b)))
})

test_that("fit() restores train mode after get_gradients()", {
  s <- make_sae()
  capture.output(s$net$fit(niter = 5))
  s$net$get_gradients()
  expect_false(s$net$model$training)
  capture.output(s$net$fit(niter = 5))
  expect_true(s$net$model$training)
})
