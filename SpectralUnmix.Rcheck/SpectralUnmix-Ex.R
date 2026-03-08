pkgname <- "SpectralUnmix"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SpectralUnmix')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("cube_helpers")
### * cube_helpers

flush(stderr()); flush(stdout())

### Name: cube_to_matrix
### Title: Helpers for Cube Reshaping and Map Extraction
### Aliases: cube_to_matrix matrix_to_cube component_map component_spectrum
###   component_reconstruction

### ** Examples

demo <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
x <- cube_to_matrix(demo$cube)
rebuilt <- matrix_to_cube(x, demo$nx, demo$ny)
dim(rebuilt)
## Not run: 
##D fit <- spectral_unmix(x, k = 3, niter = 50)
##D img <- component_map(fit, demo$nx, demo$ny, 2)
##D spec2 <- component_spectrum(fit, 2)
##D comp2 <- component_reconstruction(fit, 2, nx = demo$nx, ny = demo$ny)
## End(Not run)



cleanEx()
nameEx("demos")
### * demos

flush(stderr()); flush(stdout())

### Name: simulate_ifu_cube
### Title: Synthetic IFU Demos
### Aliases: simulate_ifu_cube demo_ifu_unmix

### ** Examples

demo <- simulate_ifu_cube()
str(demo, max.level = 1)
## Not run: 
##D result <- demo_ifu_unmix(plot = TRUE, niter = 150)
## End(Not run)



cleanEx()
nameEx("model_methods")
### * model_methods

flush(stderr()); flush(stdout())

### Name: basis
### Title: Standard Model Methods for SpectralUnmix Fits
### Aliases: basis basis.spectral_unmix coef.spectral_unmix
###   fitted.spectral_unmix predict.spectral_unmix residuals.spectral_unmix
###   summary.spectral_unmix print.summary.spectral_unmix

### ** Examples

fake <- list(
  spatial = matrix(runif(12), 4, 3),
  spectra = matrix(runif(15), 3, 5),
  reconstruction = matrix(runif(20), 4, 5),
  loss = c(3, 2, 1)
)
class(fake) <- "spectral_unmix"

basis(fake)
coef(fake)
fitted(fake)
summary(fake)

x <- matrix(runif(20), 4, 5)
residuals(fake, x = x)



cleanEx()
nameEx("spectral_unmix")
### * spectral_unmix

flush(stderr()); flush(stdout())

### Name: spectral_unmix
### Title: Spectral Unmixing for Hyperspectral and IFU Data
### Aliases: spectral_unmix print.spectral_unmix

### ** Examples

cube <- simulate_ifu_cube(nx = 12, ny = 12, n_wave = 80)
## Not run: 
##D fit <- spectral_unmix(cube$matrix, k = 3, niter = 150, lr = 0.03)
##D print(fit)
##D component_map(fit, nx = cube$nx, ny = cube$ny, component = 1)
## End(Not run)



cleanEx()
nameEx("visualization")
### * visualization

flush(stderr()); flush(stdout())

### Name: plot.spectral_unmix
### Title: Visualization Helpers for Spectral Unmixing Results
### Aliases: plot.spectral_unmix plot_reconstruction

### ** Examples

## Not run: 
##D demo <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
##D fit <- spectral_unmix(demo$matrix, k = 3, niter = 50)
##D 
##D plot(fit, type = "spectra", wavelength = demo$wavelength)
##D plot(fit, type = "maps", nx = demo$nx, ny = demo$ny)
##D plot(fit, type = "loss")
##D plot_reconstruction(fit, demo$matrix, n = 4, wavelength = demo$wavelength)
## End(Not run)



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
