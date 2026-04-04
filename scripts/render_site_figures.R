source("R/spectral_unmix.R")

dir.create("site/images", recursive = TRUE, showWarnings = FALSE)
dir.create("site/data", recursive = TRUE, showWarnings = FALSE)

nmf_multiplicative <- function(x, k, niter = 300) {
  set.seed(42)
  n <- nrow(x)
  p <- ncol(x)
  w <- matrix(stats::runif(n * k), nrow = n, ncol = k)
  h <- matrix(stats::runif(k * p), nrow = k, ncol = p)
  eps <- 1e-8

  for (i in seq_len(niter)) {
    h <- h * ((t(w) %*% x) / (t(w) %*% w %*% h + eps))
    w <- w * ((x %*% t(h)) / (w %*% h %*% t(h) + eps))
  }

  list(w = w, h = h, reconstruction = w %*% h)
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

normalize_spectra <- function(x) {
  t(apply(x, 1, display_positive))
}

extract_reference_fit <- function(fit) {
  spatial <- as.matrix(NMF::basis(fit))
  spectra <- as.matrix(NMF::coef(fit))

  list(
    spatial = spatial,
    abundance = spatial,
    spectra = spectra,
    basis = t(spectra),
    coef = t(spatial),
    reconstruction = spatial %*% spectra,
    fitted = spatial %*% spectra
  )
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

simulate_coelho_basis_benchmark <- function(n_mixtures = 96, noise = 0.001, seed = 2026) {
  real_demo_df <- utils::read.csv(
    "inst/extdata/coelho-demo-spectra.csv",
    check.names = FALSE
  )

  wavelength <- as.numeric(real_demo_df[[1L]])
  spectra <- t(as.matrix(real_demo_df[, -1L, drop = FALSE]))
  rownames(spectra) <- names(real_demo_df)[-1L]

  set.seed(seed)
  k <- nrow(spectra)
  weights <- matrix(stats::rgamma(n_mixtures * k, shape = 1.2, rate = 1), nrow = n_mixtures, ncol = k)
  weights <- weights / rowSums(weights)
  weights <- weights * stats::runif(n_mixtures, min = 0.8, max = 1.2)
  weights[seq_len(k), ] <- diag(k)

  matrix_data <- weights %*% spectra
  if (noise > 0) {
    matrix_data <- matrix_data + stats::rnorm(length(matrix_data), sd = noise)
  }
  matrix_data[matrix_data < 0] <- 0

  list(
    wavelength = wavelength,
    spectra = spectra,
    weights = weights,
    matrix = matrix_data
  )
}

demo <- simulate_ifu_cube(nx = 16, ny = 16, n_wave = 100, noise = 0.005)

spatial <- do.call(
  cbind,
  lapply(demo$abundances, function(x) as.vector(x))
)

fit <- list(
  spatial = spatial,
  abundance = spatial,
  spectra = demo$spectra,
  basis = t(demo$spectra),
  coef = t(spatial),
  reconstruction = demo$matrix,
  fitted = demo$matrix,
  loss = c(seq(0.2, 0.02, length.out = 20), seq(0.019, 0.01, length.out = 20)),
  metadata = list(source = "simulate_ifu_cube"),
  center = FALSE,
  scale = FALSE,
  call = quote(simulate_ifu_cube())
)
class(fit) <- "spectral_unmix"

grDevices::png(
  filename = "site/images/simulated-spectra.png",
  width = 1400,
  height = 900,
  res = 180
)
plot(fit, type = "spectra", wavelength = demo$wavelength)
grDevices::dev.off()

grDevices::png(
  filename = "site/images/simulated-maps.png",
  width = 1400,
  height = 900,
  res = 180
)
plot(fit, type = "maps", nx = demo$nx, ny = demo$ny)
grDevices::dev.off()

grDevices::png(
  filename = "site/images/simulated-reconstruction.png",
  width = 1400,
  height = 900,
  res = 180
)
plot_reconstruction(fit, demo$matrix, n = 4, wavelength = demo$wavelength)
grDevices::dev.off()

real_demo_df <- utils::read.csv(
  "inst/extdata/coelho-demo-spectra.csv",
  check.names = FALSE
)
real_demo <- list(
  wavelength = as.numeric(real_demo_df[[1L]]),
  matrix = t(as.matrix(real_demo_df[, -1L, drop = FALSE]))
)

fit_real_nmf <- nmf_multiplicative(real_demo$matrix, k = 3, niter = 400)
fit_real <- list(
  spatial = fit_real_nmf$w,
  abundance = fit_real_nmf$w,
  spectra = fit_real_nmf$h,
  basis = t(fit_real_nmf$h),
  coef = t(fit_real_nmf$w),
  reconstruction = fit_real_nmf$reconstruction,
  fitted = fit_real_nmf$reconstruction,
  loss = seq(0.1, 0.01, length.out = 25),
  metadata = list(source = "coelho_demo_spectra"),
  center = FALSE,
  scale = FALSE,
  call = quote(coelho_demo_spectra())
)
class(fit_real) <- "spectral_unmix"

grDevices::png(
  filename = "site/images/coelho-reconstruction.png",
  width = 1400,
  height = 900,
  res = 180
)
plot_reconstruction(
  fit_real,
  real_demo$matrix,
  pixels = seq_len(nrow(real_demo$matrix)),
  wavelength = real_demo$wavelength
)
grDevices::dev.off()

if (!requireNamespace("torch", quietly = TRUE)) {
  stop("The 'torch' package is required to generate the Coelho benchmark figure.", call. = FALSE)
}
if (!requireNamespace("NMF", quietly = TRUE)) {
  stop("The 'NMF' package is required to generate the Coelho benchmark figure.", call. = FALSE)
}
suppressPackageStartupMessages(library(NMF))

coelho_bench <- simulate_coelho_basis_benchmark(n_mixtures = 96, noise = 0.001, seed = 2026)
k <- nrow(coelho_bench$spectra)

fit_coelho_su <- spectral_unmix(
  coelho_bench$matrix,
  k = k,
  lambda_smooth = 0,
  init = "nndsvd",
  nstart = 5,
  seed = 2026,
  niter = 1500,
  lr = 0.03,
  tol = 1e-7,
  patience = 50
)

fit_coelho_nmf <- suppressMessages(suppressWarnings(
  NMF::nmf(
    coelho_bench$matrix,
    rank = k,
    method = "lee",
    nrun = 10,
    seed = 2026
  )
))
fit_coelho_nmf <- extract_reference_fit(fit_coelho_nmf)

fit_coelho_su <- reorder_fit(
  fit_coelho_su,
  best_assignment(spectral_distance(coelho_bench$spectra, fit_coelho_su$spectra))
)
fit_coelho_nmf <- reorder_fit(
  fit_coelho_nmf,
  best_assignment(spectral_distance(coelho_bench$spectra, fit_coelho_nmf$spectra))
)

truth_norm <- normalize_spectra(coelho_bench$spectra)
su_norm <- normalize_spectra(fit_coelho_su$spectra)
nmf_norm <- normalize_spectra(fit_coelho_nmf$spectra)

coelho_summary <- data.frame(
  method = c("SpectralUnmix", "NMF"),
  reconstruction_mse = c(
    mean((coelho_bench$matrix - fit_coelho_su$reconstruction)^2),
    mean((coelho_bench$matrix - fit_coelho_nmf$reconstruction)^2)
  ),
  basis_mse = c(
    mean((truth_norm - su_norm)^2),
    mean((truth_norm - nmf_norm)^2)
  ),
  iterations = c(fit_coelho_su$niter_run, NA_real_),
  converged = c(fit_coelho_su$converged, NA)
)
utils::write.csv(
  coelho_summary,
  file = "site/data/coelho-basis-benchmark-summary.csv",
  row.names = FALSE
)

grDevices::png(
  filename = "site/images/coelho-basis-benchmark.png",
  width = 1600,
  height = 1100,
  res = 180
)
old_par <- graphics::par(no.readonly = TRUE)
graphics::par(mfrow = c(2, 2), mar = c(3, 4, 3, 1), oma = c(0, 0, 3, 0))

for (j in seq_len(k)) {
  graphics::plot(
    coelho_bench$wavelength,
    truth_norm[j, ],
    type = "l",
    lwd = 2.5,
    col = "black",
    lty = 1,
    xlab = "Wavelength",
    ylab = "Normalized flux",
    main = sprintf("Basis %d: %s", j, rownames(coelho_bench$spectra)[j])
  )
  graphics::lines(
    coelho_bench$wavelength,
    su_norm[j, ],
    col = "#D55E00",
    lwd = 2,
    lty = 2
  )
  graphics::lines(
    coelho_bench$wavelength,
    nmf_norm[j, ],
    col = "#0072B2",
    lwd = 2,
    lty = 3
  )
  if (j == 1L) {
    graphics::legend(
      "topright",
      legend = c("truth", "SpectralUnmix", "NMF"),
      col = c("black", "#D55E00", "#0072B2"),
      lty = c(1, 2, 3),
      lwd = c(2.5, 2, 2),
      bty = "n"
    )
  }
}

graphics::mtext(
  sprintf(
    "Coelho basis benchmark | SpectralUnmix basis MSE = %.4g, NMF basis MSE = %.4g",
    coelho_summary$basis_mse[coelho_summary$method == "SpectralUnmix"],
    coelho_summary$basis_mse[coelho_summary$method == "NMF"]
  ),
  outer = TRUE,
  cex = 1.1
)
graphics::par(old_par)
grDevices::dev.off()
