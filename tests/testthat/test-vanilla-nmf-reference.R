if (!exists("spectral_unmix", mode = "function")) {
  candidate_paths <- c(
    file.path(getwd(), "R", "spectral_unmix.R"),
    file.path(getwd(), "..", "..", "R", "spectral_unmix.R")
  )
  source_path <- candidate_paths[file.exists(candidate_paths)][1]
  if (is.na(source_path)) {
    stop("Could not locate R/spectral_unmix.R for source-based tests.", call. = FALSE)
  }
  sys.source(source_path, envir = environment())
}

set_all_seeds <- function(seed) {
  set.seed(seed)
  if (requireNamespace("torch", quietly = TRUE)) {
    torch::torch_manual_seed(seed)
  }
}

display_positive <- function(x) {
  s <- max(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) {
    return(rep(0, length(x)))
  }

  x / s
}

spectral_distance <- function(truth, fit) {
  outer(
    seq_len(nrow(truth)),
    seq_len(nrow(fit)),
    Vectorize(function(i, j) {
      mean((display_positive(truth[i, ]) - display_positive(fit[j, ]))^2)
    })
  )
}

best_assignment <- function(dist_mat) {
  k <- nrow(dist_mat)
  perms <- as.matrix(expand.grid(rep(list(seq_len(k)), k)))
  perms <- perms[apply(perms, 1, function(z) length(unique(z)) == k), , drop = FALSE]

  scores <- apply(perms, 1, function(p) {
    sum(vapply(seq_len(k), function(i) dist_mat[i, p[i]], numeric(1)))
  })

  perms[which.min(scores), ]
}

reorder_fit <- function(fit, perm) {
  fit$spatial <- fit$spatial[, perm, drop = FALSE]
  fit$abundance <- fit$spatial
  fit$spectra <- fit$spectra[perm, , drop = FALSE]
  fit$basis <- t(fit$spectra)
  fit$coef <- t(fit$spatial)
  fit$reconstruction <- fit$spatial %*% fit$spectra
  fit$fitted <- fit$reconstruction
  fit
}

fit_best_spectral_unmix <- function(x, seeds, ...) {
  best_fit <- NULL
  best_mse <- Inf

  for (seed in seeds) {
    set_all_seeds(seed)
    fit <- spectral_unmix(x, ...)
    mse <- mean((x - fit$reconstruction)^2)
    if (mse < best_mse) {
      best_fit <- fit
      best_mse <- mse
    }
  }

  best_fit
}

fit_best_reference <- function(x, rank, seeds) {
  best_fit <- NULL
  best_mse <- Inf

  for (seed in seeds) {
    set.seed(seed)
    raw_fit <- suppressMessages(suppressWarnings(
      NMF::nmf(x, rank = rank, method = "lee", nrun = 1, seed = seed)
    ))
    fit <- list(
      spatial = as.matrix(NMF::basis(raw_fit)),
      abundance = as.matrix(NMF::basis(raw_fit)),
      spectra = as.matrix(NMF::coef(raw_fit))
    )
    fit$reconstruction <- fit$spatial %*% fit$spectra
    fit$fitted <- fit$reconstruction

    mse <- mean((x - fit$reconstruction)^2)
    if (mse < best_mse) {
      best_fit <- fit
      best_mse <- mse
    }
  }

  best_fit
}

spectra_mse <- function(truth, fit) {
  mean((t(apply(truth, 1, display_positive)) - t(apply(fit, 1, display_positive)))^2)
}

test_that("spectral_unmix recovers an exact low-rank toy factorization", {
  skip_if_not_installed("torch")

  toy <- simulate_identifiable_nmf_toy(nx = 8, ny = 8, n_wave = 60, noise = 0)
  fit_seeds <- c(11, 23, 37, 41, 53, 67)

  fit <- fit_best_spectral_unmix(
    toy$matrix,
    seeds = fit_seeds,
    k = 3,
    lambda_smooth = 0,
    niter = 2000,
    lr = 0.03,
    tol = 0
  )

  fit <- reorder_fit(fit, best_assignment(spectral_distance(toy$spectra, fit$spectra)))

  recon_mse <- mean((toy$matrix - fit$reconstruction)^2)
  spec_mse <- spectra_mse(toy$spectra, fit$spectra)

  expect_lt(recon_mse, 1e-5)
  expect_lt(spec_mse, 5e-2)
})

test_that("spectral_unmix is comparable to the NMF package on the same toy problem", {
  skip_if_not_installed("torch")
  skip_if_not_installed("NMF")

  toy <- simulate_identifiable_nmf_toy(nx = 8, ny = 8, n_wave = 60, noise = 0)
  fit_seeds <- c(11, 23, 37, 41, 53, 67)

  fit_su <- fit_best_spectral_unmix(
    toy$matrix,
    seeds = fit_seeds,
    k = 3,
    lambda_smooth = 0,
    niter = 2000,
    lr = 0.03,
    tol = 0
  )
  fit_ref <- fit_best_reference(toy$matrix, rank = 3, seeds = fit_seeds)

  fit_su <- reorder_fit(fit_su, best_assignment(spectral_distance(toy$spectra, fit_su$spectra)))
  fit_ref <- reorder_fit(fit_ref, best_assignment(spectral_distance(toy$spectra, fit_ref$spectra)))

  recon_su <- mean((toy$matrix - fit_su$reconstruction)^2)
  recon_ref <- mean((toy$matrix - fit_ref$reconstruction)^2)
  spec_su <- spectra_mse(toy$spectra, fit_su$spectra)
  spec_ref <- spectra_mse(toy$spectra, fit_ref$spectra)

  expect_lte(recon_su, recon_ref * 1.25 + 1e-6)
  expect_lte(spec_su, spec_ref * 2 + 1e-6)
})

test_that("spectral_unmix reports convergence diagnostics and supports early stopping", {
  skip_if_not_installed("torch")

  toy <- simulate_identifiable_nmf_toy(nx = 8, ny = 8, n_wave = 60, noise = 0)

  set_all_seeds(101)
  fit_full <- spectral_unmix(
    toy$matrix,
    k = 3,
    lambda_smooth = 0,
    niter = 400,
    lr = 0.03,
    tol = 0,
    patience = 400
  )

  set_all_seeds(101)
  fit_early <- spectral_unmix(
    toy$matrix,
    k = 3,
    lambda_smooth = 0,
    niter = 400,
    lr = 0.03,
    tol = 1e-2,
    patience = 5
  )

  expect_type(fit_early$loss, "double")
  expect_identical(fit_early$niter_run, length(fit_early$loss))
  expect_true(is.logical(fit_early$converged))
  expect_true(fit_early$converged)
  expect_lt(fit_early$niter_run, fit_full$niter_run)
  expect_lt(utils::tail(fit_early$loss, 1L), fit_early$loss[1L])
  expect_lt(mean((toy$matrix - fit_early$reconstruction)^2), 5e-4)
})

test_that("summary and print expose convergence information", {
  fake <- list(
    spatial = matrix(runif(12), 4, 3),
    abundance = matrix(runif(12), 4, 3),
    spectra = matrix(runif(15), 3, 5),
    reconstruction = matrix(runif(20), 4, 5),
    fitted = matrix(runif(20), 4, 5),
    loss = c(3, 2, 1),
    niter_run = 3L,
    converged = TRUE
  )
  class(fake) <- "spectral_unmix"

  smry <- summary(fake)

  expect_identical(smry$n_iter, 3L)
  expect_true(smry$converged)
  expect_match(paste(capture.output(print(fake)), collapse = "\n"), "converged: yes")
  expect_match(paste(capture.output(print(smry)), collapse = "\n"), "converged: yes")
})
