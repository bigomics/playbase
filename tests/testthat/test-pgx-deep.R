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

test_that("history records every loss component and continues across fit() calls", {
  s <- make_sae()
  capture.output(s$net$fit(niter = 5))
  capture.output(s$net$fit(niter = 3))
  h <- s$net$history
  expect_named(h, c("iteration", "split", "component", "value"))
  expect_type(h$iteration, "integer")
  expect_true(all(is.finite(h$value)))

  train <- h[h$split == "train", ]
  expect_equal(sort(unique(train$iteration)), 1:8)
  components <- c("classification", paste0("reconstruction:", c("gx", "px", "mx")), "l1", "l2", "total")
  for (i in 1:8) {
    expect_setequal(train$component[train$iteration == i], components)
  }

  test <- h[h$split == "test", ]
  expect_equal(test$iteration[test$component == "classification"], 1:8)
})

test_that("legacy loss_history/val_history are derived from history", {
  s <- make_sae()
  capture.output(s$net$fit(niter = 4))
  h <- s$net$history
  expect_equal(s$net$loss_history, h$value[h$split == "train" & h$component == "classification"])
  expect_equal(s$net$val_history, h$value[h$split == "test" & h$component == "classification"])
  expect_length(s$net$val_history, 4)
  pdf(NULL)
  on.exit(dev.off())
  expect_no_error(s$net$plot_loss())
})

test_that("validation does not touch BatchNorm running stats and keeps train mode", {
  s <- make_sae()
  capture.output(s$net$fit(niter = 2))
  s$net$model$train(TRUE)
  bufs <- function() lapply(s$net$model$buffers, as.numeric)
  before <- bufs()
  expect_true(any(grepl("running_mean", names(before))))
  v <- s$net$validate()
  expect_identical(bufs(), before)
  expect_true(s$net$model$training)
  expect_true(all(c("classification", "reconstruction:gx") %in% names(v)))
})

test_that("validate() returns NULL without a validation split", {
  s <- make_sae(validation_ratio = 0)
  expect_null(s$net$validate())
  capture.output(s$net$fit(niter = 2))
  expect_false(any(s$net$history$split == "test"))
})

## mofa.augment lives with the MOFA code; source it when playbase is not loaded
if (!exists("mofa.augment")) {
  sys.source(file.path("..", "..", "R", "pgx-mofa.R"), envir = globalenv())
}

test_that("augment = 0 keeps the previous train/test column counts", {
  s <- make_sae()
  expect_equal(nrow(s$net$x_train[[1]]), 19)
  expect_equal(nrow(s$net$x_test[[1]]), 5)
  expect_equal(nrow(s$net$xx[[1]]), 24)
  expect_named(s$net$split, c("sample", "label", "split"))
  expect_equal(s$net$split$sample, colnames(s$fx$X[[1]]))
  expect_equal(s$net$split$label, as.character(s$fx$Y$class))
  expect_equal(sum(s$net$split$split == "test"), 5)
})

test_that("test split is stratified by class", {
  s <- make_sae()
  sp <- s$net$split
  ntest <- table(factor(sp$label[sp$split == "test"], levels = unique(sp$label)))
  expected <- 5 * table(sp$label) / nrow(sp)
  expect_true(all(abs(ntest - expected[names(ntest)]) < 1))
  expect_true(all(ntest >= 1))
  expect_equal(as.numeric(s$net$y_test[[1]]), as.integer(factor(sp$label[sp$split == "test"], levels = s$net$labels$class)))
})

test_that("augmentation copies only training samples and tests on clean originals", {
  s <- make_sae(augment = 10, augment_z = 0.1, validation_ratio = 0.25)
  sp <- s$net$split
  is_test <- sp$split == "test"
  expect_equal(nrow(s$net$x_train[[1]]), 10 * sum(!is_test))
  expect_equal(nrow(s$net$x_test[[1]]), sum(is_test))
  expect_equal(as.numeric(s$net$y_train[[1]]), rep(as.numeric(factor(sp$label[!is_test], levels = s$net$labels$class)), 10))
  ## test rows are exactly the scaled originals (no noise)
  xs <- t(scale(t(s$net$X$gx)))
  expect_equal(as.matrix(s$net$x_test$gx), t(xs[, is_test]), ignore_attr = TRUE, tolerance = 1e-6)
  ## no training row is a noisy copy of a test sample: each training row is
  ## closer to some training original than to every test original
  xtr <- as.matrix(s$net$x_train$gx)
  d <- as.matrix(dist(rbind(t(xs), xtr)))[-(1:24), 1:24]
  nearest <- apply(d, 1, which.min)
  expect_false(any(is_test[nearest]))
  expect_equal(nearest, rep(which(!is_test), 10), ignore_attr = TRUE)
  capture.output(s$net$fit(niter = 2))
  expect_true(any(s$net$history$split == "test"))
})

test_that("singleton classes stay in training", {
  fx <- make_sae_fixture()
  fx$Y$class <- factor(c(as.character(fx$Y$class[-24]), "lonely"))
  net <- suppressMessages(MultiOmicsSAE$new(fx$X, fx$Y, ntop = 0))
  expect_equal(net$split$split[net$split$label == "lonely"], "train")
  tab <- table(net$split$label[net$split$split == "test"])
  expect_true(all(tab[paste0("class", 1:3)] >= 1))
})

test_that("fit() keeps one optimizer across calls and recreates it on a method change", {
  s <- make_sae()
  expect_null(s$net$optimizer)
  capture.output(s$net$fit(niter = 2))
  opt <- s$net$optimizer
  expect_s3_class(opt, "optim_adam")
  capture.output(s$net$fit(niter = 2))
  expect_identical(s$net$optimizer, opt)
  expect_message(capture.output(s$net$fit(niter = 2, optim = "sgd")), "switching optimizer")
  expect_s3_class(s$net$optimizer, "optim_sgd")
})

## noise and dropout draw from the torch RNG, which a checkpoint does not
## carry, so exact resumption is checked without them
test_that("training resumed from a saved checkpoint matches uninterrupted training", {
  weights <- function(net) lapply(net$model$state_dict(), as.array)

  a <- make_sae(add_noise = FALSE)
  capture.output(a$net$fit(niter = 10))

  b <- make_sae(add_noise = FALSE)
  capture.output(b$net$fit(niter = 5))
  ckpt <- b$net$checkpoint()
  expect_named(ckpt, c("model", "optimizer", "optim", "iteration", "history", "args"))
  expect_equal(ckpt$iteration, 5L)
  f <- tempfile(fileext = ".pt")
  on.exit(unlink(f))
  torch::torch_save(ckpt, f)
  capture.output(b$net$fit(niter = 2)) ## must not alter the saved state

  c <- make_sae(add_noise = FALSE)
  c$net$restore(torch::torch_load(f))
  expect_equal(c$net$iteration, 5L)
  capture.output(c$net$fit(niter = 5))

  expect_equal(weights(c$net), weights(a$net), tolerance = 1e-6)
  h <- c$net$history
  expect_equal(unique(h$iteration[h$split == "train"]), 1:10)
  expect_equal(c$net$loss_history, a$net$loss_history, tolerance = 1e-6)
})

test_that("checkpoint() copies tensors so later training leaves it unchanged", {
  s <- make_sae()
  capture.output(s$net$fit(niter = 2))
  ckpt <- s$net$checkpoint()
  before <- lapply(ckpt$model, as.array)
  capture.output(s$net$fit(niter = 2))
  expect_identical(lapply(ckpt$model, as.array), before)
})

test_that("restore() refuses a checkpoint built with different arguments", {
  s <- make_sae()
  capture.output(s$net$fit(niter = 2))
  ckpt <- s$net$checkpoint()
  other <- make_sae(dropout = 0.1)
  expect_error(other$net$restore(ckpt), "different arguments: dropout")
})

test_that("augmentation noise is on the raw scale, as the legacy board added it", {
  fx <- make_sae_fixture()
  ## one view with a low- and a high-variance feature block
  fx$X$gx[1:15, ] <- fx$X$gx[1:15, ] * 0.1
  fx$X$gx[16:30, ] <- fx$X$gx[16:30, ] * 10
  set.seed(3)
  net <- suppressMessages(MultiOmicsSAE$new(fx$X, fx$Y, ntop = 0, augment = 5, augment_z = 1))
  ntr <- sum(net$split$split == "train")
  xtr <- as.matrix(net$x_train$gx)
  ## every copy is noisy (mofa.augment), so compare against copy 1
  first <- xtr[seq_len(ntr), , drop = FALSE]
  noise <- xtr - first[rep(seq_len(ntr), 5), , drop = FALSE]
  sd_low <- mean(apply(noise[, 1:15], 2, sd))
  sd_high <- mean(apply(noise[, 16:30], 2, sd))
  ## one view-wide raw noise level => relatively much more noise on the
  ## low-variance features once scaled by their own SD
  expect_gt(sd_low, 10 * sd_high)
})
