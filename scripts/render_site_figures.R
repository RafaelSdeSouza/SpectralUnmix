source("R/spectral_unmix.R")

dir.create("site/images", recursive = TRUE, showWarnings = FALSE)

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
