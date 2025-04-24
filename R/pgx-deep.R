##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2024 BigOmics Analytics SA. All rights reserved.
##
##


## ======================================================================
## ==================== TORCH ===========================================
## ======================================================================


#'
#'
#' @export
MultiOmicsSAE <- R6::R6Class(
  "MultiOmicsSAE",
  private = list(),
  public = list(
    X = NULL,
    Y = NULL,
    nviews = NULL,
    labels = NULL,
    model = NULL,
    xx = NULL,
    x_train = NULL,
    y_train = NULL,
    x_test = NULL,
    y_test = NULL,
    sdx = NULL,
    loss_weights = NULL,
    sd_weight = NULL,
    loss_history = NULL,
    val_history = NULL,
    initialize = function(X, Y,
                          model = "MT", ntop = 1000,
                          num_layers = list(c(0.5, 40), c(0.5, 0.25), c(0.5)),
                          dropout = 0, scale = TRUE, validation_ratio = 0.2,
                          loss_weights = c(y = 1, ae = 10, l1 = 0.1, l2 = 0.01),
                          actfun = c("relu", "leaky", "gelu", "silu")[1],
                          use_glu = 1, use_bn = TRUE, add_noise = 0,
                          sd_weight = TRUE, device = "cpu") {
      if (!is.null(ntop) && ntop > 0) {
        message("reducing to maximum ", ntop, " features")
        X <- lapply(X, function(x) head(x[order(-matrixStats::rowSds(x, na.rm = TRUE)),,drop=FALSE ], ntop))
      }

      if (any(sapply(X, function(x) sum(is.na(x))) > 0)) {
        message("missing features detected! imputing.")
        X <- lapply(X, function(x) svdImpute2(x))
      }

      ## add tiny noise: for some reason this prevents NaN in some
      ## calculations.
      for (i in 1:length(X)) {
        X[[i]] <- X[[i]] + 1e-4 * matrix(
          rnorm(length(X[[i]])),
          nrow(X[[i]]), ncol(X[[i]])
        )
      }

      ## must be factor for now
      Y <- lapply(Y, function(y) factor(as.character(y)))

      ##
      self$X <- X
      self$Y <- Y
      self$nviews <- length(X)
      self$sdx <- lapply(X, matrixStats::rowSds, na.rm = TRUE)
      self$sd_weight <- sd_weight

      if (scale) {
        for (i in 1:length(X)) {
          X[[i]] <- X[[i]] - rowMeans(X[[i]], na.rm = TRUE)
          X[[i]] <- X[[i]] / self$sdx[[i]]
        }
      }

      ## normalized torch input. This is the input the network sees.
      self$xx <- lapply(X, function(x) torch::torch_tensor(t(x)))

      ## split in train and validation
      if (validation_ratio > 0) {
        r <- (1 - validation_ratio)
        xdim <- ncol(X[[1]])
        ii <- sample(1:xdim, r * xdim)
        self$x_train <- lapply(X, function(x) t(x[, ii, drop=FALSE]))
        self$x_test <- lapply(X, function(x) t(x[, -ii, drop=FALSE]))
        self$y_train <- lapply(Y, function(y) torch::torch_tensor(y[ii]))
        self$y_test <- lapply(Y, function(y) torch::torch_tensor(y[-ii]))
      } else {
        self$x_train <- lapply(X, function(x) t(x))
        self$x_test <- NULL
        self$y_train <- lapply(Y, function(y) torch::torch_tensor(y))
        self$y_test <- NULL
      }
      self$x_train <- lapply(self$x_train, torch::torch_tensor)
      self$x_test <- lapply(self$x_test, torch::torch_tensor)

      self$Y <- lapply(Y, factor)
      self$labels <- lapply(self$Y, levels)

      if (device == "gpu" && !torch::cuda_is_available()) {
        message("Cuda not available. Switching to cpu.")
        device <- "cpu"
      }

      if (model == "MT") {
        message("creating MultiBlockMultiTargetSAE")
        self$model <- MultiBlockMultiTargetSAE_module(
          self$x_train,
          self$y_train,
          num_layers = num_layers,
          dropout = dropout,
          actfun = actfun,
          use_glu = use_glu,
          use_bn = use_bn,
          add_noise = add_noise,
          device = device
        )
      } else if (model == "simple") {
        message("creating SimpleSAE")
        self$model <- SimpleSAE_module(
          self$x_train,
          self$y_train,
          latent_dim = 20,
          actfun = actfun
        )
      } else {
        stop("unknown model")
      }

      ##
      if (device == "gpu") self$to_device("gpu")

      ## set optimizer
      if (length(loss_weights) != 4) stop("loss_weights wrong length (must be 4)")
      self$loss_weights <- loss_weights
    },
    fit = function(niter = 200, optim = c("adam", "adamw", "sgd", "lbfgs")[1]) {
      model <- self$model
      optim <- optim[1]

      if (optim == "adam") {
        optimizer <- torch::optim_adam(model$parameters, lr = 0.001)
      } else if (optim == "adamw") {
        optimizer <- torch::optim_adamw(
          model$parameters,
          lr = 0.001, weight_decay = self$loss_weights["l2"]
        )
      } else if (optim == "sgd") {
        optimizer <- torch::optim_sgd(model$parameters, lr = 0.001)
      } else if (optim == "lbfgs") {
        optimizer <- torch::optim_lbfgs(model$parameters, lr = 1)
      } else {
        stop("unknown optimization method ", optim)
      }

      xx <- self$x_train
      yy <- self$y_train
      vx <- self$x_test
      vy <- self$y_test

      calc_loss <- function() {
        ## forward compute
        output <- model(xx)

        ## compute fit loss
        nnf_cross_entropy.NA <- function(a, b) {
          ii <- which(!is.na(as.numeric(b)))
          torch::nnf_cross_entropy(a[ii], b[ii])
        }
        y_loss <- Reduce("+", mapply(nnf_cross_entropy.NA, output$y_pred, yy))

        if (!is.null(vx)) {
          ## compute validation loss
          voutput <- model(vx)
          val_loss <- Reduce("+", mapply(nnf_cross_entropy.NA, voutput$y_pred, vy))
        } else {
          val_loss <- 0
        }

        ae_loss <- Reduce("+", mapply(torch::nnf_mse_loss, output$ae, xx))
        enc_weights <- lapply(model$encoder, function(m) m$parameters[[1]])
        l1_loss <- Reduce("+", sapply(enc_weights, function(w) torch::torch_norm(w, 1)))
        l2_loss <- Reduce("+", sapply(enc_weights, function(w) torch::torch_norm(w, 2)))
        ## loss <- 1*y_loss + 1e1*ae_loss + 0.1*l1_loss + 0.01*l2_loss
        w <- self$loss_weights
        loss <- w["y"] * y_loss + w["ae"] * ae_loss + w["l1"] * l1_loss + w["l2"] * l2_loss

        # Store loss for plotting
        self$loss_history <- c(self$loss_history, y_loss$item())
        self$val_history <- c(self$val_history, val_loss$item())

        optimizer$zero_grad()
        loss$backward()
        loss
      }

      for (t in 1:niter) {
        if (t %% (niter / 10) == 0) cat(".")
        optimizer$step(calc_loss)
        # calc_loss()
        # optimizer$step()
      }
    },
    predict = function(X = NULL) {
      if (is.null(X)) {
        X <- self$X
      }
      xx <- lapply(X, function(x) scale(t(x)))
      xx <- lapply(xx, torch::torch_tensor)
      output <- self$model(xx)
      y_pred <- list()
      for (i in 1:length(output$y_pred)) {
        y_pred[[i]] <- as.matrix(output$y_pred[[i]])
        rownames(y_pred[[i]]) <- colnames(X[[1]])
        colnames(y_pred[[i]]) <- self$labels[[i]]
      }
      y_pred
    },
    get_redux = function(xx = NULL) {
      if (is.null(xx)) {
        xx <- self$xx ## full matrix
      } else {
        xx <- lapply(xx, function(x) scale(t(x)))
        xx <- lapply(xx, torch::torch_tensor)
      }
      output <- self$model(xx)
      redux <- c(output$enc, "multi-omics" = output$combined)
      redux <- lapply(redux, as.matrix)
      redux
    },
    get_gradients = function(sd_weight = NULL) {
      self$get_gradients.autograd(sd_weight = sd_weight)
    },
    get_gradients.delta = function(sd_weight = NULL) {
      sdx <- self$sdx
      x_train <- self$x_train
      y_train <- self$y_train
      if (is.null(sd_weight)) sd_weight <- self$sd_weight
      ## calculate feature perturbation response
      grad <- vector("list", length(self$Y))
      names(grad) <- names(self$Y)
      k <- names(x_train)[1]
      for (k in names(x_train)) {
        x <- x_train[[k]]
        delta <- rbind(0, diag(ncol(x)))
        xx <- lapply(x_train, function(x) matrix(0, nrow(delta), ncol(x)))
        xx[[k]] <- delta
        xx <- lapply(xx, torch::torch_tensor)
        lapply(xx, dim)

        ## call forward model
        self$model$train(FALSE)
        output <- self$model(xx)
        m <- names(output$y_pred)[1]
        for (m in names(output$y_pred)) {
          z <- as.matrix(output$y_pred[[m]])
          colnames(z) <- self$labels[[m]]
          rownames(z) <- c("ref", rownames(self$X[[k]]))
          dz <- t(t(z) - z["ref", ])
          dz <- dz[-1, , drop = FALSE] ## remove ref
          if (sd_weight) {
            ## dz <- dz[names(sdx[[k]]),] * sdx[[k]]  ## SD is really important!!
            dz <- dz * sdx[[k]][rownames(dz)] ## SD is really important!!
          }
          grad[[m]][[k]] <- dz
        }
      }
      return(grad)
    },
    get_gradients.autograd = function(sd_weight = NULL) {
      sdx <- self$sdx
      x_train <- self$x_train
      y_train <- self$y_train
      if (is.null(sd_weight)) sd_weight <- self$sd_weight
      ## calculate feature perturbation response
      delta <- lapply(self$x_train, function(x) 0 * x[1, , drop = FALSE])
      for (i in 1:length(delta)) delta[[i]]$requires_grad <- TRUE
      self$model$train(FALSE)
      output <- self$model(delta)
      grad <- vector("list", length(self$Y))
      names(grad) <- names(self$Y)
      k <- 1
      for (k in 1:length(output$y_pred)) {
        y_pred <- output$y_pred[[k]]
        grad_k <- sapply(1:ncol(y_pred), function(i) {
          torch::autograd_grad(y_pred[1, i], inputs = delta, retain_graph = TRUE)
        })
        grad_k <- lapply(grad_k, as.numeric)
        ny <- length(self$labels[[k]])
        nx <- length(self$x_train)
        for (j in 1:nx) {
          m <- names(self$x_train)[j]
          ii <- j + ((1:ny) - 1) * nx
          grad_km <- do.call(cbind, grad_k[ii])
          if (sd_weight) {
            grad_km <- grad_km * self$sdx[[m]]
          }
          rownames(grad_km) <- rownames(self$X[[m]])
          colnames(grad_km) <- self$labels[[k]]
          grad[[k]][[m]] <- grad_km
        }
      }
      return(grad)
    },
    plot_loss = function() {
      x <- self$loss_history
      y <- self$val_history
      x <- pmax(x, 0.1 * min(x[x > 0], na.rm = TRUE)) ## avoid real zero
      y <- pmax(y, 0.1 * min(y[y > 0], na.rm = TRUE))
      x <- log(x)
      y <- log(y)
      ylim <- c(min(x, y, na.rm = TRUE), max(x, y, na.rm = TRUE))
      plot(x, ylim = ylim, xlab = "iteration", ylab = "loss (log)")
      points(y, col = "green3")
      legend("topright",
        legend = c("training", "validation"),
        pch = NULL, fill = c("black", "green3"), cex = 0.9
      )
    },
    to_device = function(dev) {
      if (dev == "gpu" && torch::cuda_is_available()) {
        message("moving parameters to gpu")
        self$model$cuda()
        self$xx <- lapply(self$xx, function(x) x$cuda())
        self$x_train <- lapply(self$x_train, function(x) x$cuda())
        self$x_test <- lapply(self$x_test, function(x) x$cuda())
        self$y_train <- lapply(self$y_train, function(x) x$cuda())
        self$y_test <- lapply(self$y_test, function(x) x$cuda())
      } else {
        message("moving parameters to cpu")
        self$model$cpu()
        self$xx <- lapply(self$xx, function(x) x$cpu())
        self$x_train <- lapply(self$x_train, function(x) x$cpu())
        self$x_test <- lapply(self$x_test, function(x) x$cpu())
        self$y_train <- lapply(self$y_train, function(x) x$cpu())
        self$y_test <- lapply(self$y_test, function(x) x$cpu())
      }
    },
    get_dims = function() {
      getdims <- function(e) {
        pp <- lapply(e$parameters, dim)
        pp <- pp[sapply(pp, length) == 2]
        as.vector(sapply(pp, head, 1))
      }
      encoder.dims <- lapply(self$model$encoder, function(e) getdims(e))
      input.dims <- sapply(self$xx, ncol)
      encoder.dims <- mapply(c, input.dims, encoder.dims, SIMPLIFY = FALSE)
      decoder.dims <- lapply(self$model$decoder, function(e) getdims(e))
      concat_dim <- sum(sapply(encoder.dims, tail, 1))
      integrator.dims <- unlist(c(concat_dim, getdims(self$model$integrator)))
      predictor.dims <- lapply(self$model$predictor, function(e) getdims(e))
      list(
        encoder = encoder.dims,
        decoder = decoder.dims,
        integrator = integrator.dims,
        predictor = predictor.dims
      )
    }
  )
)


#'
#' @export
nn_gaussian_noise <- torch::nn_module(
  "gaussian_noise",
  initialize = function(sd = 0.1, device = "cpu") {
    self$sd <- sd
    self$device <- device
  },
  forward = function(x) {
    rx <- self$sd * torch::torch_randn(c(dim(x)[1], dim(x)[2]))
    if (self$device == "gpu") rx <- rx$cuda()
    x + rx
  }
)

# Prepare network module
MultiBlockMultiTargetSAE_module <- torch::nn_module(
  "MultiBlockMultiTargetSAE_module",
  initialize = function(xx, yy,
                        dropout = 0,
                        actfun = c("relu", "leaky", "gelu", "silu")[1],
                        use_glu = 1,
                        use_bn = TRUE,
                        add_noise = 0,
                        num_layers = list(c(0.5, 40), c(0.5, 0.25), c(0.5, 0.25)),
                        device = "cpu") {
    self$encoder <- torch::nn_module_list()
    self$decoder <- torch::nn_module_list()
    self$predictor <- torch::nn_module_list()
    self$yy <- yy
    self$device <- device

    activationFunc <- torch::nn_relu
    if (actfun == "leaky") activationFunc <- torch::nn_leaky_relu
    if (actfun == "gelu") activationFunc <- torch::nn_gelu
    if (actfun == "silu") activationFunc <- torch::nn_silu

    ## custom input processing layers
    add_input_mods <- function(mods, dim, use_bn, dropout, add_noise) {
      if (use_bn) mods <- c(mods, torch::nn_batch_norm1d(dim))
      if (add_noise > 0) mods <- c(mods, nn_gaussian_noise(add_noise, device = device))
      if (dropout > 0) mods <- c(mods, torch::nn_dropout(dropout))
      mods
    }

    ae_dims <- list()
    k <- names(xx)[1]
    for (k in names(xx)) {
      xdim <- ncol(xx[[k]])
      mods <- list()
      dims <- num_layers[[1]]
      dims <- round(ifelse(dims < 1, dims * xdim, dims))
      dims <- pmax(dims, tail(dims, 1)) ## always larger than bottlenect
      dims <- pmin(dims, ceiling(xdim/2))  ## always smaller than input (??)
      dims <- c(xdim, dims)
      nlayers <- length(dims)
      for (i in 1:(nlayers - 1)) {
        mods <- add_input_mods(mods, dims[i], use_bn, dropout, add_noise)
        if (use_glu == 2) {
          mods <- c(mods, torch::nn_linear(dims[i], 2 * dims[i + 1])) ## input
          mods <- c(mods, torch::nn_glu())
        } else {
          mods <- c(mods, torch::nn_linear(dims[i], dims[i + 1]))
        }
      }
      ae_dims[[k]] <- dims
      if (use_glu > 1) {
        message("Creating GLU encoding layer: ", paste(ae_dims[[k]], collapse = " "))
      } else {
        message("Creating encoding layer: ", paste(ae_dims[[k]], collapse = " "))
      }
      self$encoder[[k]] <- torch::nn_sequential(!!!mods)
    }
    ae_dims <- lapply(ae_dims, rev)

    ## create decoder counterpart for auto-encoder using reverse
    ## number of layers.
    for (k in names(xx)) {
      kdim <- ncol(xx[[k]])
      mods <- list()
      num_ae <- length(ae_dims[[k]]) - 1
      for (i in 1:num_ae) {
        kdim <- ae_dims[[k]][i]
        outdim <- ae_dims[[k]][i + 1]
        mods <- add_input_mods(mods, kdim, use_bn, dropout, add_noise)
        if (use_glu == 2) {
          mods <- c(mods, torch::nn_linear(kdim, 2 * outdim)) ## input
          mods <- c(mods, torch::nn_glu())
        } else {
          mods <- c(mods, torch::nn_linear(kdim, outdim))
        }
        if (i < num_ae) mods <- c(mods, activationFunc())
      }
      if (use_glu > 1) {
        message("Creating GLU decoding layer: ", paste(ae_dims[[k]], collapse = " "))
      } else {
        message("Creating decoding layer: ", paste(ae_dims[[k]], collapse = " "))
      }
      self$decoder[[k]] <- torch::nn_sequential(!!!mods)
    }

    ## integration layers: create n-layer integration NN network
    ## concatenating latent vectors of each datatype.
    merge_dim <- sum(sapply(ae_dims, head, 1))
    integration_dim <- merge_dim
    self$integrator <- NULL
    if (length(num_layers[[2]])) {
      dims <- c(merge_dim, num_layers[[2]])
      dims <- round(ifelse(dims < 1, dims * merge_dim, dims))
      nlayers <- length(dims)
      mods <- list()
      for (i in 1:(nlayers - 1)) {
        mods <- add_input_mods(mods, dims[i], use_bn, dropout, add_noise)
        if (use_glu > 0) {
          mods <- c(mods, torch::nn_linear(dims[i], 2 * dims[i + 1])) ## input
          mods <- c(mods, torch::nn_glu())
        } else {
          mods <- c(mods, torch::nn_linear(dims[i], dims[i + 1])) ## input
        }
        mods <- c(mods, activationFunc())
      }
      if (use_glu > 0) {
        message("Creating GLU integration layer: ", paste(dims, collapse = " "))
      } else {
        message("Creating integration layer: ", paste(dims, collapse = " "))
      }
      self$integrator <- torch::nn_sequential(!!!mods)
      integration_dim <- tail(dims, 1) ## last dim from integration layers
    } else {
      message("Skipping integration layer")
    }

    ## predictor layers
    m <- names(yy)[1]
    for (m in names(yy)) {
      yval <- as.numeric(yy[[m]])
      ydim <- length(unique(setdiff(yval, NA)))
      dims <- c(integration_dim, num_layers[[3]], ydim)
      dims <- round(ifelse(dims < 1, dims * integration_dim, dims))
      ## if(length(dims)>1) dims[-1] <- pmax(dims[-1], ydim)  ## really?
      nlayers <- length(dims)
      mods <- list()
      for (i in 1:(nlayers - 1)) {
        mods <- add_input_mods(mods, dims[i], use_bn, dropout, add_noise)
        if (use_glu > 0) {
          mods <- c(mods, torch::nn_linear(dims[i], 2 * dims[i + 1])) ## input
          mods <- c(mods, torch::nn_glu())
        } else {
          mods <- c(mods, torch::nn_linear(dims[i], dims[i + 1])) ## input
        }
        if (i < nlayers) mods <- c(mods, activationFunc()) ## layer1
      }
      if (use_glu > 0) {
        message("Creating GLU prediction layer: ", paste(dims, collapse = " "))
      } else {
        message("Creating prediction layer: ", paste(dims, collapse = " "))
      }
      self$predictor[[m]] <- torch::nn_sequential(!!!mods)
    }
  },
  forward = function(xx) {
    enc_outputs <- list()
    ae_outputs <- list()
    pred_outputs <- list()
    for (k in names(xx)) {
      enc_outputs[[k]] <- self$encoder[[k]](xx[[k]])
      ae_outputs[[k]] <- self$decoder[[k]](enc_outputs[[k]])
    }
    combined_output <- torch::torch_cat(enc_outputs, dim = 2)
    if (!is.null(self$integrator)) {
      combined_output <- self$integrator(combined_output)
    }
    for (m in names(self$yy)) {
      pred_outputs[[m]] <- self$predictor[[m]](combined_output)
    }
    list(
      y_pred = pred_outputs,
      enc = enc_outputs,
      ae = ae_outputs,
      combined = combined_output
    )
  }
)

# Very simple model of supervised multi-view auto-encoder. Kept for
# illustration and baseline benchmarking.
#
#
SimpleSAE_module <- torch::nn_module(
  "SimpleSAE_module",
  initialize = function(xx, yy, latent_dim = 20,
                        actfun = c("relu", "leaky", "gelu")) {
    self$encoder <- torch::nn_module_list()
    self$decoder <- torch::nn_module_list()
    self$predictor <- torch::nn_module_list()
    self$yy <- yy

    activationFunc <- torch::nn_relu
    if (actfun == "leaky") activationFunc <- torch::nn_leaky_relu
    if (actfun == "gelu") activationFunc <- torch::nn_gelu

    for (k in names(xx)) {
      xdim <- ncol(xx[[k]])
      message("creating encoder for ", k, ": ", paste(xdim, latent_dim, collapse = " "))
      self$encoder[[k]] <- torch::nn_sequential(
        torch::nn_linear(ncol(xx[[k]]), latent_dim),
        activationFunc()
      )
    }
    for (k in names(xx)) {
      xdim <- ncol(xx[[k]])
      message("creating decoder for ", k, ": ", paste(latent_dim, xdim, collapse = " "))
      self$decoder[[k]] <- torch::nn_sequential(
        torch::nn_linear(latent_dim, xdim)
      )
    }
    nviews <- length(xx)
    nmerge <- nviews * latent_dim
    mlp_dim <- round(nmerge / 2)
    message("creating integrator layer: ", paste(nmerge, mlp_dim, collapse = " "))
    self$integrator <- torch::nn_sequential(
      torch::nn_linear(nmerge, mlp_dim),
      activationFunc()
    )
    for (m in names(yy)) {
      ydim <- length(unique(as.numeric(yy[[m]])))
      message("creating predictor layer: ", paste(mlp_dim, ydim, collapse = " "))
      self$predictor[[m]] <- torch::nn_sequential(
        torch::nn_linear(mlp_dim, ydim)
      )
    }
  },
  forward = function(xx) {
    enc_output <- list()
    ae_output <- list()
    for (k in names(xx)) {
      enc_output[[k]] <- self$encoder[[k]](xx[[k]])
      ae_output[[k]] <- self$decoder[[k]](enc_output[[k]])
    }
    combined_output <- self$integrator(torch::torch_cat(enc_output, dim = 2))
    pred_output <- list()
    for (m in names(self$yy)) {
      pred_output[[m]] <- self$predictor[[m]](combined_output)
    }
    list(
      y_pred = pred_output,
      enc = enc_output,
      ae = ae_output,
      combined = combined_output
    )
  }
)


## ======================================================================
## =========================== FUNCTIONS ================================
## ======================================================================

#'
#'
#' @export
deep.plotBiomarkerHeatmap <- function(net, datatypes = NULL, balanced = TRUE,
                                      annot = NULL, labels=NULL, ntop = 50, ...) {
  grad <- net$get_gradients()

  if (!is.null(datatypes)) {
    for (k in 1:length(grad)) {
      sel <- (names(grad[[k]]) %in% datatypes)
      grad[[k]] <- grad[[k]][sel]
    }
  }

  if (balanced) {
    ## ranking per phenotype condition
    rnk <- vector("list", length(grad[[1]]))
    for (k in 1:length(grad)) {
      gg <- mofa.prefix(grad[[k]])
      ## ranking per datatype
      rr <- lapply(gg, function(g) rank(-rowMeans(g**2)))
      rnk <- mapply(cbind, rnk, rr, SIMPLIFY = FALSE)
    }
    ## now we have ranking per datatype across multiple phenotypes
    rnk2 <- lapply(rnk, function(r) rownames(r)[order(rowMeans(r))])
    rnk2 <- lapply(rnk2, function(x) head(rep(x,ntop),ntop))
    rnk2 <- do.call(rbind,rnk2)    
    sel  <- head(unique(as.vector(rnk2)),ntop)   
  } else {
    ## not balancing between datatypes. just largest gradient.
    mgrad <- lapply(grad, function(gr) mofa.merge_data(gr))
    pgrad <- do.call(cbind, lapply(mgrad, function(m) rowMeans(m**2)))
    sel <- rownames(pgrad)[head(order(-rowMeans(pgrad)), ntop)]
  }
  X <- do.call(rbind, mofa.prefix(net$X))
  colnames(X) <- make_unique(colnames(X))
  X <- X[sel, ]

  Y <- data.frame(net$Y)
  rownames(Y) <- colnames(X)
  if (!is.null(annot)) {
    kk <- setdiff(colnames(annot), colnames(Y))
    Y <- cbind(Y, annot[, kk, drop = FALSE])
  }
  
  if(!is.null(labels)) {
    rownames(X) <- labels[match(rownames(X),names(labels))]
  }

  gx.splitmap(X,
    col.annot = Y, split = 1, mar = c(2, 6),
    show_key = TRUE, annot.cex = 1.2, ...
  )
}


#'
#'
#' @export
deep.getConfusionMatrix <- function(net, what = c("train", "test")[1]) {
  matlist <- list()
  pheno <- names(net$labels)[1]
  for (pheno in names(net$labels)) {
    labels <- net$labels[[pheno]]
    if (what == "test" && !is.null(net$x_test)) {
      inputx <- net$x_test
      y <- labels[as.numeric(net$y_test[[pheno]])]
    } else {
      inputx <- net$x_train
      y <- labels[as.numeric(net$y_train[[pheno]])]
    }
    output <- net$model(inputx)
    y_out <- as.matrix(output$y_pred[[pheno]])
    y_pred <- labels[max.col(y_out)]
    mat <- table(actual = y, y_pred, useNA = "always")
    mat <- mat[
      match(labels, rownames(mat)),
      match(labels, colnames(mat))
    ]
    rownames(mat) <- colnames(mat) <- labels
    mat[is.na(mat)] <- 0
    matlist[[pheno]] <- mat
  }
  return(matlist)
}

#'
#'
#' @export
deep.plotAutoEncoderReconstructions <- function(net, dtypes = NULL,
                                                main = NULL, par = TRUE) {
  output <- net$model(net$xx)
  str(output)
  if (is.null(dtypes)) {
    dtypes <- names(net$X)
  } else {
    intersect(dtypes, c("mixed", names(net$X)))
  }
  if (length(dtypes) == 0) {
    return(NULL)
  }
  if (par) {
    nc <- ceiling(sqrt(length(dtypes)))
    nr <- ceiling(length(dtypes) / nc)
    par(mfrow = c(nr, nc))
  }
  if ("mixed" %in% dtypes) {
    if (is.null(main)) main <- "mixed features"
    dtypes <- names(net$X)
    sX <- t(do.call(rbind, net$X[dtypes]))
    pX <- lapply(output$ae[dtypes], as.matrix)
    pX <- do.call(cbind, pX)
    sX <- scale(sX)
    jj <- sample(length(sX), 1000, replace = TRUE)
    plot(sX[jj], pX[jj],
      cex = 0.5, main = main,
      xlab = "actual", ylab = "reconstructed"
    )
  } else {
    for (k in dtypes) {
      sX <- scale(t(net$X[[k]]))
      pX <- as.matrix(output$ae[[k]])
      jj <- sample(length(sX), 1000, replace = TRUE)
      tt <- k
      if (!is.null(main)) tt <- main
      plot(sX[jj], pX[jj],
        cex = 0.5, main = tt,
        xlab = "actual", ylab = "reconstructed"
      )
    }
  }
}


#'
#'
#' @export
deep.plotRedux <- function(net, pheno = NULL, method = "tsne", views = NULL, par = TRUE, cex = 1) {
  redux <- net$get_redux(xx = NULL)
  if(any(sapply(redux,ncol)<2)) {
    sel <- which(sapply(redux,ncol)<3)
    message("[deep.plotRedux] warning: augmenting ",
            paste(names(sel),collapse=" "))
    for(j in sel) {
      Rj <- do.call(cbind,rep(list(redux[[j]]),3))
      Rj <- Rj + 1e-3*sd(Rj)*matrix(rnorm(length(Rj)),nrow(Rj),ncol(Rj))
      redux[[j]] <- cbind(redux[[j]], Rj)
    }
  }
  if (method == "tsne") {
    px <- min(30, nrow(redux[[1]]) / 4)
    redux <- lapply(redux, function(r) {
      Rtsne::Rtsne(r, perplexity = px, check_duplicates = FALSE)$Y
    })
  } else if (method == "umap") {
    px <- max(min(15, nrow(redux[[1]]) / 6), 1)
    redux <- lapply(redux, function(r) uwot::umap(r, n_neighbors = px))
  } else if (method == "pca") {
    px <- min(30, nrow(redux[[1]]) / 4)
    redux <- lapply(redux, function(x) svd(x, nv = 2)$u[, 1:2])
  } else {
    stop("Error method not supported")
  }
  y <- net$y
  if (is.null(pheno)) pheno <- 1
  if ("Y" %in% names(net)) y <- net$Y[[pheno]]
  if (!is.null(views)) redux <- redux[names(redux) %in% views]
  if (par) {
    nq <- ceiling(sqrt(length(redux)))
    par(mfrow = c(nq, nq), mar = c(4, 4, 4, 4))
  }
  for (k in names(redux)) {
    plot(redux[[k]],
      col = factor(y), main = k,
      cex = cex, pch = 19,
      xlab = paste0(method, "-x"),
      ylab = paste0(method, "-y")
    )
  }
}

#'
#'
#' @export
deep.plotMultiOmicsGradients <- function(grad, n = 20, cex.names = 1,
                                         onlypositive = FALSE,
                                         par = TRUE) {
  if (par) {
    ngroup <- ncol(grad[[1]])
    nview <- length(grad)
    par(mfrow = c(ngroup, nview), mar = c(10, 4, 1.5, 2))
  }
  for (k in 1:ncol(grad[[1]])) {
    for (i in 1:length(grad)) {
      gr <- grad[[i]][, k]
      gr[is.na(gr) | is.infinite(gr)] <- 0
      if (onlypositive) {
        gr <- head(gr[order(-gr)], n)
      } else {
        gr <- head(gr[order(-abs(gr))], n)
      }

      barplot(sort(gr), ylab = "gradient", las = 3, cex.names = cex.names)
      title(
        paste0(colnames(grad[[1]])[k],"@", names(grad)[i]),
        line = -0.5, cex.main = 1.15
      )
    }
  }
}


#'
#'
#' @export
deep.plotGradientExplained <- function(grad, main = "") {
  w <- sapply(grad, function(x) mean(x**2)**0.5)
  barplot(100 * w / sum(w),
    main = main,
    xlab = "datatype", ylab = "gradient explained (%)"
  )
}


#'
#'
#' @export
deep.plotGradientVSFoldchange <- function(grad, fc, data = FALSE, par = TRUE) {
  ## align datatypes
  dt <- intersect(names(grad), names(fc))
  grad <- grad[dt]
  fc <- fc[dt]

  ## align features
  for (k in names(grad)) {
    pp <- intersect(rownames(grad[[k]]), rownames(fc[[k]]))
    ph <- intersect(colnames(grad[[k]]), colnames(fc[[k]]))
    grad[[k]] <- grad[[k]][pp, ph, drop = FALSE]
    fc[[k]] <- fc[[k]][pp, ph, drop = FALSE]
  }

  if (data) {
    G <- do.call(rbind, mofa.prefix(grad))
    F <- do.call(rbind, mofa.prefix(fc))
    colnames(G) <- paste0("grad.", colnames(G))
    colnames(F) <- paste0("logFC.", colnames(F))
    df <- data.frame(G, F)
    avg.rnk <- rank(rowMeans(G**2)) + rank(rowMeans(F**2))
    df <- df[order(-avg.rnk), , drop = FALSE]
    return(df)
  }

  ## plot gradient vs foldchagne
  if (par) {
    ngrad <- length(grad) * ncol(grad[[1]])
    nc <- ceiling(sqrt(ngrad))
    nr <- ceiling(ngrad / nc)
    par(mfrow = c(nr, nc), mar = c(4.5, 4, 3, 2))
  }
  i <- 2
  k <- 3
  for (i in 1:ncol(grad[[1]])) {
    for (k in names(grad)) {
      f <- fc[[k]][, i]
      g <- grad[[k]][, i]
      xlab <- "logFC"
      plot(f, g, xlab = xlab, ylab = "network gradient")
      abline(h = 0, v = 0, lty = 2)
      title(paste0(colnames(grad[[k]])[i],"@", k))
      sel <- unique(c(head(order(-f), 5), head(order(-g), 5)))
      text(f[sel], g[sel], labels = names(f)[sel], col = "red", pos = 2)
      sel <- unique(c(head(order(f), 5), head(order(g), 5)))
      text(f[sel], g[sel], labels = names(f)[sel], col = "blue", pos = 2)
    }
  }
}

#'
#'
#' @export
deep.plotNeuralNet <- function(net, svgfile = NULL, rm.files=TRUE, image=NULL) {
  
  if (!dir.exists("/opt/PlotNeuralNet")) {
    message("ERROR: please install PlotNeuralNet in /opt")
    return(NULL)
  }

  if (!is.null(svgfile)) {
    if (!grepl("svg$", svgfile)) {
      message("ERROR: output file must be SVG")
      return(NULL)
    }
  }

  header_src <- glue::glue("
import sys
sys.path.append('../')
from pycore.tikzeng import *

# defined your arch
arch = [
    to_head( '..' ),
    to_cor(),
    to_begin(),
", .trim = FALSE)

  str_at <- function(at) paste0("(", at[1], ",", at[2], ",", at[3], ")")
  input_src <- function(name, at, image, caption) {
    to1 <- str_at(c(at[1] - 2.5, at[2], at[3]))
    to2 <- str_at(at)
    glue::glue("
    # input {name}
    to_input( '{image}', to='{to1}', width=9.5, height=2.2 ),
    to_Conv('{name}', '', '', offset='(-4,0,6.5)', to='{to2}', height=0, depth=0, width=0, caption='{caption}'),
", .trim = FALSE)
  }

  ## name="PX";caption="CAPTION";at=c(0,0,0);dims=c(1000,500,40,500,1000)
  autoencoder_src <- function(name, caption, at, dims) {
    to <- str_at(at)
    nlayer <- length(dims)
    ae <- paste0(name, "ae", 1:nlayer)
    w <- round((dims / dims[1])**0.66 * 45)
    midi <- (nlayer + 1) / 2
    txt <- "    # Auto-encoder {name}\n"
    for (i in 1:nlayer) {
      if (i == 1) {
        txt <- paste0(txt, "    to_Conv('{ae[1]}', '{dims[1]}', '', offset='(0,0,0)', to='{to}', height=10, depth={w[1]}, width=2, caption='{caption}'),\n")
      } else {
        if (i == midi) {
          txt <- paste0(txt, "    to_Pool('{ae[", i, "]}', offset='(2,0,0)', to='({ae[", i - 1, "]}-east)', height=10, depth={w[", i, "]}, width=2),\n")
        } else {
          txt <- paste0(txt, "    to_Conv('{ae[", i, "]}', '{dims[", i, "]}', '', offset='(2,0,0)', to='({ae[", i - 1, "]}-east)', height=10, depth={w[", i, "]}, width=2),\n")
        }
        txt <- paste0(txt, "    to_connection('{ae[", i - 1, "]}', '{ae[", i, "]}'),\n")
      }
    }
    txt <- paste0(txt, "    to_skip( of='{ae[", midi, "]}', to='sum1', pos=1.4),\n")
    glue::glue(txt, .trim = FALSE)
  }

  merge_src <- function() {
    nview <- length(views)
    to <- str_at(c(0, 0, (nview - 1) * 22 / 2))
    glue::glue("
    # Combine layer
    to_Sum('sum1', offset='(14,0,0)', to='{to}', radius=2.5, opacity=0.6),
    to_Conv('label2', '', '', offset='(15.2,0,0)', to='{to}', height=0, depth=0, width=0, caption='MERGE'),
", .trim = FALSE)
  }

  dense_src <- function(dims = c(100, 20)) {
    nlayer <- length(dims)
    w <- round((dims / dims[1])**0.5 * 35)
    txt <- "    # Dense MLP layer\n"
    for (i in 1:nlayer) {
      if (i == 1) {
        txt <- paste0(txt, "    to_Conv('dense1', '{dims[1]}', '', offset='(4,0,0)', to='(sum1-east)', height=10, depth={w[1]}, width=2, caption='integrator'),\n")
        txt <- paste0(txt, "    to_connection('sum1', 'dense", i, "'),\n")
      } else {
        txt <- paste0(txt, "    to_Conv('dense", i, "', '{dims[", i, "]}', '', offset='(2,0,0)', to='(dense", i - 1, "-east)', height=10, depth={w[", i, "]}, width=2, caption=''),\n")
        txt <- paste0(txt, "    to_connection('dense", i - 1, "', 'dense", i, "'),\n")
      }
    }
    glue::glue(txt, .trim = FALSE)
  }

  predictor_src <- function(name, dims, caption, offset, lastdense) {
    offset <- str_at(offset)
    nlayer <- length(dims)
    pred <- paste0(name, "pred", 1:nlayer)
    w <- round((dims / dims[1])**0.5 * 15)
    txt <- "    # Predictor layer {name}\n"
    for (i in 1:nlayer) {
      if (i == 1 && i != nlayer) {
        txt <- paste0(txt, "    to_Conv('{pred[1]}', '{dims[1]}', '', offset='(3,0,0)', to='(dense", lastdense, "-east)', height=10, depth={w[1]}, width=2, caption='predictor'),\n")
        txt <- paste0(txt, "    to_connection('dense", lastdense, "', '{pred[1]}'),\n")
      } else {
        if (i == nlayer && i != 1) {
          txt <- paste0(txt, "    to_ConvSoftMax( name='{pred[", i, "]}', s_filer = '{dims[", i, "]}', offset='(2,0,0)', to='({pred[", i - 1, "]}-east)', width=2, height=10, depth={w[", i, "]}, caption='{caption}'),\n")
        } else if (i == nlayer && i == 1) {
          txt <- paste0(txt, "    to_ConvSoftMax( name='{pred[", i, "]}', s_filer = '{dims[", i, "]}', offset='(2,0,0)', to='(dense", lastdense, "-east)', width=2, height=10, depth={w[", i, "]}, caption='{caption}'),\n")
        } else {
          txt <- paste0(txt, "    to_Conv('{pred[", i, "]}', '{dims[", i, "]}', '', offset='{offset}', to='({pred[", i - 1, "]}-east)', height=10, depth={w[", i, "]}, width=2, caption=''),\n")
        }
        ## add arrows
        if (i != 1) {
          txt <- paste0(txt, "    to_connection( '{pred[", i - 1, "]}', '{pred[", i, "]}'),\n")
        } else {
          txt <- paste0(txt, "    to_connection( 'dense", lastdense, "', '{pred[", i, "]}'),\n")
        }
      }
    }
    glue::glue(txt, .trim = FALSE)
  }

  end_src <- function(nview) {
    zoff <- (nview - 1) * 22
    glue::glue("
    ## margins
    to_Conv('empty1', '', '', offset='(0,0,0)', to='(-7,0,{zoff})', height=0, depth=0, width=0),
    to_Conv('empty2', '', '', offset='(14,0,0)', to='(dense1-east)', height=0, depth=0, width=0),
    to_end()
    ]

def main():
    namefile = str(sys.argv[0]).split('.')[0]
    to_generate(arch, namefile + '.tex' )

if __name__ == '__main__':
    main()
", .trim = FALSE)
  }

  ## get dimensions
  net_dims <- net$get_dims()
  views <- rev(names(net$X))
  views <- iconv(views, "latin1", "ASCII", sub = "")
  views <- gsub("[_\\$ ]","",views) ## latex special chars
  nview <- length(views)

  ## create input image if not given
  if(is.null(image) || !file.exists(image)) {
    image <- c()
    for(i in 1:nview) {
      image[[i]] = paste0("/tmp/inputimage",i,".png")
      x1 <- net$X[[i]]
      x1 <- head(x1[order(-matrixStats::rowSds(x1,na.rm=TRUE)),,drop=FALSE],100)
      x1 <- head(t(rowscale(x1)),100)
      x1 <- 1*(x1 > mean(x1,na.rm=TRUE))
      png(image[[i]],width=950,height=220)
      par(mar=c(0,0,0,0))
      gx.imagemap(x1, cex=0, col=grey.colors(64))
      dev.off()
    }
  }
  if(length(image) < nview) {
    image <- head(rep(image,nview),nview)
  }
  
  ## build code text
  ltx <- header_src
  redux <- net$get_redux()
  rdim <- ncol(redux[[1]]) ## bottleneck dimension
  targets <- names(net$Y)
  targets <- iconv(targets, "latin1", "ASCII", sub = "")
  targets <- gsub("[_\\$ ]","",targets) ## latex special chars
  
  ntargets <- length(targets)
  for (i in 1:length(views)) {
    at <- c(0, 0, (i - 1) * 22)
    caption <- toupper(views[i])
    if (max(nchar(views)) < 10) caption <- paste("datatype:~", caption)
    ltx1 <- input_src(paste0("input", views[i]), at, image[[i]], caption = caption)
    ltx <- paste(ltx, ltx1)
  }
  ltx <- paste(ltx, merge_src())
  for (i in 1:length(views)) {
    at <- c(0, 0, (i - 1) * 22)
    dims <- c(net_dims[["encoder"]][[i]], net_dims[["decoder"]][[i]])
    caption <- paste0("autoencoder~", toupper(views[i]))
    ltx1 <- autoencoder_src(views[i], at = at, caption = caption, dims = dims)
    ltx <- paste(ltx, ltx1)
  }
  ## merge dot and integration layer
  dims <- net_dims[["integrator"]]
  ltx <- paste(ltx, dense_src(dims = dims))

  ## predictors
  ntargets <- length(targets)
  lastdense <- length(net_dims$integrator)
  i <- 1
  for (i in 1:ntargets) {
    offset <- c(2, 0, 8 * ((i - 1) - (ntargets - 1) / 2))
    dims <- net_dims[["predictor"]][[i]]
    # caption <- paste0("predictor~",toupper(targets[i]))
    caption <- paste0("", toupper(targets[i]))
    ltx1 <- predictor_src(
      name = targets[i], dims = dims,
      caption = caption, offset = offset, lastdense = lastdense
    )
    ltx <- paste(ltx, ltx1)
  }

  ## end/closing text
  ltx <- paste(ltx, end_src(nview = nview))

  ## create PDF using PlotNeuralNet
  pyfile <- tempfile(fileext = ".py", tmpdir = "/opt/PlotNeuralNet/pyexamples")
  ## pyfile="/opt/PlotNeuralNet/pyexamples/model2.py"
  write(ltx, file = pyfile)

  if (!file.exists(pyfile)) {
    message("[deep.plotNeuralNet] WARNING. failed to create Python file ",pyfile)
  }

  texfile <- sub("[.]py$", ".tex", pyfile)
  auxfile <- sub("[.]py$", ".aux", pyfile)
  logfile <- sub("[.]py$", ".log", pyfile)
  pdffile <- sub("[.]py$", ".pdf", pyfile)
  if(is.null(svgfile)) {
    svgfile <- sub("[.]py$", ".svg", pyfile)
    svgfile <- paste0("/tmp/",basename(svgfile)) ## write to /tmp
  }

  if (file.exists(pyfile)) {
    cmd <- glue::glue("cd /opt/PlotNeuralNet/pyexamples/ && python {pyfile}")
    suppressMessages(system(cmd,
      ignore.stdout = TRUE, ignore.stderr = TRUE,
      intern = FALSE, show.output.on.console = FALSE
    ))
    if (!file.exists(texfile)) {
      message("[deep.plotNeuralNet] WARNING. failed to create TeX file ",texfile)
      message("[deep.plotNeuralNet] python cmd = ",cmd)
    }
  }

  if (file.exists(texfile)) {
    cmd <- glue::glue("cd /opt/PlotNeuralNet/pyexamples/ && pdflatex -interaction=batchmode {texfile}")
    suppressMessages(system(cmd,
      ignore.stdout = TRUE, ignore.stderr = TRUE,
      intern = FALSE, show.output.on.console = FALSE
    ))
    if (!file.exists(pdffile)) {      
      message("[deep.plotNeuralNet] WARNING. failed to create PDF file ",pdffile)
      message("[deep.plotNeuralNet] pdflatex cmd = ",cmd)    
    }
  }

  if (file.exists(pdffile)) {
    cmd <- glue::glue("pdf2svg {pdffile} {svgfile}")
    suppressMessages(system(cmd,
      ignore.stdout = TRUE, ignore.stderr = TRUE,
      intern = FALSE, show.output.on.console = FALSE
    ))
    if (!file.exists(svgfile)) {
      message("[deep.plotNeuralNet] WARNING. failed to create SVG file ",svgfile)
      message("[deep.plotNeuralNet] pdf2svg cmd = ",cmd)    
    }
  }

  
  if(rm.files) {
    if (file.exists(auxfile)) unlink(auxfile)
    if (file.exists(logfile)) unlink(logfile)
    if (file.exists(texfile)) unlink(texfile)
    if (file.exists(pyfile))  unlink(pyfile)
    if (file.exists(pdffile)) unlink(pdffile)
  }
  
  if (file.exists(svgfile)) {
    return(svgfile)
  } else {
    return(NULL)
  }
}

## ======================================================================
## ======================== END OF FILE =================================
## ======================================================================
