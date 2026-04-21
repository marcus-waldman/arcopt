# Low-Rank Symmetric Matrix Factorization
# ========================================
#
# Problem: minimize  f(U) = 0.5 * || U U' - M ||_F^2
#          over      U in R^{d x r}
#
# with M in R^{d x d} symmetric. When rank(M) = r the global minimum is
# f* = 0 achieved on the orbit { U : U U' = M } (a continuous manifold
# with O(r) rotational ambiguity). The origin U = 0 is a symmetric
# saddle: gradient vanishes, Hessian has 2*rank(M) negative eigenvalues
# per ambient direction.
#
# Parameter count n = d * r. For d = 50, r = 10 this is n = 500.
#
# Difficulty:
#   - Nonconvex (symmetric saddles at U = 0 and along orbit of permuted
#     eigenvector signs)
#   - Gradient of pure BFGS updates destroys negative-curvature info
#   - At n = 500 the Hessian is 500 x 500 (still tractable via
#     eigendecomposition solver, which is the point)
#
# References:
#   Ge, Jin, Netrapalli, Sidford (2015). "Escaping From Saddle Points"
#   Srebro, Jaakkola (2003). "Weighted low-rank approximations"

make_matrix_factorization_problem <- function(d = 50L, r = 10L,
                                              noise_sd = 0,
                                              u_true_seed = 1L) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else {
    NULL
  }
  set.seed(u_true_seed)
  u_true <- matrix(stats::rnorm(d * r), d, r)
  m_target <- tcrossprod(u_true)
  if (noise_sd > 0) {
    noise <- matrix(stats::rnorm(d * d, sd = noise_sd), d, d)
    m_target <- m_target + 0.5 * (noise + t(noise))
  }
  if (!is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = .GlobalEnv)
  }

  fn <- function(x) {
    u_mat <- matrix(x, d, r)
    resid_mat <- tcrossprod(u_mat) - m_target
    0.5 * sum(resid_mat * resid_mat)
  }

  gr <- function(x) {
    u_mat <- matrix(x, d, r)
    resid_mat <- tcrossprod(u_mat) - m_target
    as.vector(2 * resid_mat %*% u_mat)
  }

  hess <- function(x) {
    u_mat <- matrix(x, d, r)
    gram <- crossprod(u_mat)
    resid_mat <- tcrossprod(u_mat) - m_target
    n <- d * r
    h <- matrix(0, n, n)
    eye_d <- diag(d)
    for (i1 in seq_len(r)) {
      row_idx <- ((i1 - 1L) * d + 1L):(i1 * d)
      for (i2 in seq_len(r)) {
        col_idx <- ((i2 - 1L) * d + 1L):(i2 * d)
        block <- 2 * gram[i1, i2] * eye_d +
          2 * tcrossprod(u_mat[, i2], u_mat[, i1])
        if (i1 == i2) {
          block <- block + 2 * resid_mat
        }
        h[row_idx, col_idx] <- block
      }
    }
    h
  }

  list(
    d = d,
    r = r,
    n = d * r,
    u_true = u_true,
    m_target = m_target,
    fn = fn,
    gr = gr,
    hess = hess
  )
}

make_start_point <- function(d, r, seed, sd = 0.1) {
  set.seed(seed)
  as.vector(matrix(stats::rnorm(d * r, sd = sd), d, r))
}
