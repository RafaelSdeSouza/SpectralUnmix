source("R/spectral_unmix.R")

dir.create("site/images", recursive = TRUE, showWarnings = FALSE)

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
