##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' Goeman's global test on gene sets
#'
#' Self-contained test as in MetaboAnalyst's QEA. Same p-values, bit for bit, as
#' \code{globaltest::gt(y, t(X), subsets, model = "logistic")}: the linear score
#' test Q = Y'XX'Y/Y'Y on Y = y - mean(y), Imhof tail, generalised-F fallback
#' for tiny p. Ported so we do not depend on the (dormant) globaltest package.
#'
#' The test is non-directional; the score takes its sign from the members' mean
#' difference. A set whose members move both ways is therefore counted as up
#' or down by its mean in the per-method \code{sig.counts}.
#'
#' @param X Numeric matrix, features (rows) by samples. NA is mean-imputed as
#'   gt does.
#' @param y Binary group per sample: 0/1, logical or two-level factor.
#' @param gmt Named list of gene sets (vectors of feature names).
#' @param standardize Scale features to unit variance first, as
#'   \code{gt(standardize = TRUE)}, so high-variance features do not dominate.
#'   Default \code{FALSE} (gt default).
#'
#' @return Data frame in \code{gmt} order, rownames \code{names(gmt)}:
#'   \code{score} (sign of mean difference * -log10 p, p floored at 1e-300),
#'   \code{p.value}, \code{q.value} (BH), \code{SIZE} (members in X). NA if no
#'   member is in X.
#'
#' @export
gset.globaltest <- function(X, y, gmt, standardize = FALSE) {
  y <- if (is.factor(y)) as.integer(y) - 1 else as.numeric(y)
  X <- as.matrix(X)
  stopifnot(ncol(X) == length(y), !anyNA(y))
  ## members as row positions (gmt order, as intersect); one rownames hash
  idx <- match(unlist(gmt, use.names = FALSE), rownames(X))
  idx <- split(idx, factor(rep(seq_along(gmt), lengths(gmt)), seq_along(gmt)))
  idx <- lapply(idx, function(i) unique(i[!is.na(i)]))
  ok <- lengths(idx) > 0
  dx <- rowMeans(X[, y == 1, drop = FALSE], na.rm = TRUE) -
    rowMeans(X[, y == 0, drop = FALSE], na.rm = TRUE)
  pv <- score <- rep(NA_real_, length(gmt))
  if (any(ok)) {
    pv[ok] <- .gt.pvalues(X, y, idx[ok], standardize)
    dir <- sapply(idx[ok], function(i) sign(mean(dx[i], na.rm = TRUE)))
    score[ok] <- dir * -log10(pmax(pv[ok], 1e-300))
  }
  nm <- if (is.null(names(gmt))) NULL else make.unique(names(gmt))
  data.frame(
    score = score,
    p.value = pv,
    q.value = stats::p.adjust(pv, method = "BH"),
    SIZE = replace(lengths(idx), !ok, NA),
    row.names = nm
  )
}

## gt's p-value per set of row indices of X. Null fit, centring and X'Y are done
## once; the BLAS calls mirror globaltest, so eigenvalues and p-values are
## bit-identical.
.gt.pvalues <- function(X, y, idx, standardize = FALSE) {
  n <- length(y)
  ## gt imputes the lm(x ~ 1) fit; all-NA -> 0
  for (j in which(rowSums(is.na(X)) > 0)) {
    ok <- !is.na(X[j, ])
    X[j, !ok] <- if (any(ok)) {
      stats::lm.fit(matrix(1, sum(ok), 1), X[j, ok])$coefficients
    } else {
      0
    }
  }
  if (standardize) { ## population sd, zero variance -> 1
    vars <- apply(X, 1, stats::var) * (n - 1) / n
    X <- X * sqrt(rep(1, nrow(X)) / replace(vars, vars == 0, 1))
  }
  Z <- matrix(1, n, 1)
  fit <- suppressWarnings(stats::glm.fit(Z, y, family = stats::binomial()))
  Y <- fit$y - fit$fitted.values
  sumYY <- sum(Y * Y)
  X <- t(X)
  ## centre; samples x features
  X <- X - Z %*% solve(crossprod(Z), crossprod(Z, X))
  csm <- colSums(X * X)
  ## constant features -> exactly 0. Upstream quirk, kept for bit-identity: the
  ## cut-off is global, so a feature ~1e7x smaller in SD than the largest
  ## anywhere is zeroed too (p = 1 for its sets). Moot on log data or with
  ## standardize = TRUE.
  X[, csm < max(csm) * 1e-14] <- 0
  zero <- colSums(X != 0) == 0
  Xt <- t(X)
  xy <- crossprod(X, Y)
  one <- function(i) {
    if (all(zero[i])) {
      return(1)
    }
    big <- length(i) > n ## use the smaller Gram matrix: same non-zero spectrum
    G <- crossprod(if (big) Xt[i, , drop = FALSE] else X[, i, drop = FALSE])
    S <- if (sumYY == 0) {
      0
    } else if (big) {
      sum(Y * (G %*% Y)) / sumYY
    } else {
      sum(xy[i] * xy[i]) / sumYY
    }
    lams <- c(
      eigen(G, symmetric = TRUE, only.values = TRUE)$values,
      numeric(max(0, n - length(i)))
    )
    lams[-n] <- lams[-n] - S ## P(Q >= S) = P(sum lams_i chi2_1 >= 0)
    .gt.getP(lams)
  }
  d <- duplicated(idx) ## identical sets are tested once
  pv <- numeric(length(idx))
  pv[!d] <- vapply(idx[!d], one, 1)
  ## match() on lists is slow: twins only
  dd <- which(d | duplicated(idx, fromLast = TRUE))
  pv[dd] <- pv[dd][match(idx[dd], idx[dd])]
  pv
}

## P(sum lams_i chi2_1 >= 0): F for 2 distinct values, else Imhof, else
## generalised F.
.gt.getP <- function(lams) {
  if (all(lams >= 0)) {
    return(1)
  }
  if (all(lams <= 0)) {
    return(0)
  }
  lams <- lams[lams != 0]
  if (length(unique(lams)) == 2) {
    return(.gt.genF(lams))
  }
  m <- length(lams)
  lams2 <- lams^2
  K <- sum(abs(lams)) / 20
  g <- function(v) { ## Imhof (1961) integrand at 0, (0, Inf) mapped to (0, 1)
    u <- -log(1 - v) / K
    U <- rep(u, each = m)
    th <- 0.5 * .colSums(atan(lams * U), m, length(u))
    rho <- exp(.colSums(0.25 * log(1 + lams2 * U^2), m, length(u)))
    f <- sin(th) / (u * rho)
    0.5 + replace(f, u == 0, sum(lams) / 2) / (pi * K * (1 - v))
  }
  im <- tryCatch(
    stats::integrate(g, 0, 1, rel.tol = sqrt(.Machine$double.eps)),
    error = function(e) NULL
  )
  if (is.null(im) || im$value <= 10 * im$abs.error) .gt.genF(lams) else im$value
}

## Generalised F: each side a Robbins-Pitman (Ruben) mixture of scaled chi2.
## Small eigenvalues are weeded into blocks >= max/acc, coarser if no beta gives
## a mixture.
.gt.genF <- function(lams, eps = 1e-10) {
  mix <- function(lam, acc = 50) {
    repeat {
      lam <- sort(lam, decreasing = TRUE)
      for (k in (length(lam) - 1):0) {
        rest <- lam[seq_len(k)]
        join <- lam[seq_along(lam) > k]
        nl <- if (k == 0) 1 else min(sum(join) %/% (lam[1] / acc), length(join))
        join <- rep(sum(join) / nl, nl)
        if (k == 0 || (length(join) > 0 && rest[k] >= join[1])) break
      }
      lam <- c(rest, join)
      b0 <- 2 * max(lam) * min(lam) / (max(lam) + min(lam))
      bb <- seq(b0, min(lam), length = 10)
      for (b in bb) {
        cc <- .gt.getC(lam, b, eps)
        if (!is.null(cc)) {
          return(list(c = cc, beta = b, m = length(lam)))
        }
      }
      acc <- mean(c(acc, 1))
    }
  }
  P <- mix(lams[lams > 0])
  N <- mix(-lams[lams < 0])
  df1 <- P$m + 2 * (seq_along(P$c) - 1)
  df2 <- rep(N$m + 2 * (seq_along(N$c) - 1), each = length(df1))
  pfs <- stats::pf(N$beta / P$beta * (df2 / df1), df1, df2, lower.tail = FALSE)
  drop(crossprod(P$c, matrix(pfs, length(df1))) %*% N$c)
}

## Ruben series coefficients; NULL if not a valid (non-negative) mixture.
.gt.getC <- function(lams, beta, eps = 1e-10) {
  lams <- sort(lams)
  cc <- prod(sqrt(beta / lams))
  if (cc == 0 || cc > 1) {
    return(NULL)
  }
  rest <- 1 - cc
  mx <- cc
  db <- 1 - beta / lams
  d <- NULL
  ix <- 1
  repeat {
    d <- c(d, 0.5 * sum(db^ix))
    cc <- c(cc, mean(cc[1:ix] * d[ix:1]))
    rest <- if (rest > 100 * .Machine$double.neg.eps) {
      rest - cc[ix + 1]
    } else {
      rest * lams[1] / lams[length(lams)]
    }
    mx <- max(mx, cc[ix + 1])
    if (rest / mx < eps || cc[ix + 1] < 0) break
    ix <- ix + 1
  }
  if (rest > 0) cc <- c(cc, rest)
  if (any(cc < 0)) NULL else cc
}
