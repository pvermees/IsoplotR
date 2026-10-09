# Pure-R non-negative least squares (Lawson & Hanson, 1974, Algorithm NNLS)
#
# Solves  min ||A x - b||_2   subject to  x >= 0
#
# Drop-in replacement for nnls::nnls(A, b). Returns an object with the same
# main components: x, deviance, residuals, residualNorm, reduced (the dual /
# gradient vector w = A'(b - Ax)), and mode (1 = success, 2 = bad dimensions,
# 3 = iteration limit reached).

nnls_r <- function(A, b, max_iter = 3L * ncol(A)) {
  A <- as.matrix(A)
  b <- as.numeric(b)
  m <- nrow(A)
  n <- ncol(A)

  if (m != length(b) || m == 0L || n == 0L) {
    stop("nnls_r: dimensions of A and b do not match (or are empty)")
  }
  if (anyNA(A) || anyNA(b)) stop("nnls_r: missing values are not allowed")

  eps <- .Machine$double.eps
  col_norm <- sqrt(colSums(A^2))

  # Unconstrained least squares on the passive (free) set, zeros elsewhere
  solve_passive <- function(P) {
    z <- numeric(n)
    if (!any(P)) return(z)        # empty passive set: all zeros
    AP <- A[, P, drop = FALSE]
    cf <- qr.coef(qr(AP), b)
    cf[is.na(cf)] <- 0            # rank-deficient: aliased coefficients -> 0
    z[P] <- cf
    z
  }

  x <- numeric(n)
  P <- logical(n)                 # TRUE = in passive set, FALSE = held at zero
  w <- as.numeric(crossprod(A, b))  # gradient (x = 0 initially)
  excluded <- logical(n)          # columns rejected as linearly dependent
  iter <- 0L
  mode <- 1L

  repeat {
    # ---- Outer loop: pick the zero-set variable with the largest gradient ----
    cand <- which(!P & !excluded)
    if (length(cand) == 0L || all(w[cand] <= 0)) break

    j <- cand[which.max(w[cand])]

    # Reject the column if it is (numerically) dependent on the passive set
    if (any(P)) {
      AP <- A[, P, drop = FALSE]
      cf <- replace_na0(qr.coef(qr(AP), A[, j]))
      r <- A[, j] - AP %*% cf
      indep <- sqrt(sum(r^2)) > 100 * eps * col_norm[j]
    } else {
      indep <- col_norm[j] > 0
    }
    if (!indep) {
      excluded[j] <- TRUE
      next
    }

    P[j] <- TRUE
    iter <- iter + 1L
    if (iter > max_iter) { mode <- 3L; P[j] <- FALSE; break }

    z <- solve_passive(P)

    # ---- Inner loop: restore feasibility if any passive coefficient <= 0 ----
    while (any(P & z <= 0)) {
      iter <- iter + 1L
      if (iter > max_iter) { mode <- 3L; break }

      idx <- which(P & z <= 0)
      ratio <- x[idx] / (x[idx] - z[idx])
      alpha <- min(ratio)
      x <- x + alpha * (z - x)

      # Move to the zero set every passive variable that has hit zero
      drop <- P & (x <= 0 | seq_len(n) == idx[which.min(ratio)])
      x[drop] <- 0
      P[drop] <- FALSE

      z <- solve_passive(P)
    }
    if (mode == 3L) break

    x <- z
    excluded[] <- FALSE           # a new passive set may make columns usable
    w <- as.numeric(crossprod(A, b - A %*% x))
  }

  fitted <- as.numeric(A %*% x)
  resid <- b - fitted
  w <- as.numeric(crossprod(A, resid))

  structure(
    list(
      x = x,
      deviance = sum(resid^2),
      residuals = resid,
      residualNorm = sqrt(sum(resid^2)),
      reduced = w,
      mode = mode,
      passive = which(P)
    ),
    class = "nnls"
  )
}

replace_na0 <- function(v) { v[is.na(v)] <- 0; v }
