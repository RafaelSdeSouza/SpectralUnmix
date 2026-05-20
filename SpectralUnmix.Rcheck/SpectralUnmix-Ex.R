pkgname <- "SpectralUnmix"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('SpectralUnmix')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("basis")
### * basis

flush(stderr()); flush(stdout())

### Name: basis
### Title: NMF-style Basis Matrix
### Aliases: basis

### ** Examples

fake <- list(spectra = matrix(runif(15), 3, 5))
class(fake) <- "spectral_unmix"
basis(fake)



cleanEx()
nameEx("coef.spectral_unmix")
### * coef.spectral_unmix

flush(stderr()); flush(stdout())

### Name: coef.spectral_unmix
### Title: Extract Abundance Coefficients
### Aliases: coef.spectral_unmix

### ** Examples

fake <- list(spatial = matrix(runif(12), 4, 3))
class(fake) <- "spectral_unmix"
coef(fake)



cleanEx()
nameEx("coelho_demo_spectra")
### * coelho_demo_spectra

flush(stderr()); flush(stdout())

### Name: coelho_demo_spectra
### Title: Demo Stellar Spectra Extracted From the Coelho Library
### Aliases: coelho_demo_spectra

### ** Examples

demo <- coelho_demo_spectra()
dim(demo$matrix)
demo$metadata



cleanEx()
nameEx("coelho_stellar_subset")
### * coelho_stellar_subset

flush(stderr()); flush(stdout())

### Name: coelho_stellar_subset
### Title: Curated Stellar Library Subset for NMF Demos
### Aliases: coelho_stellar_subset

### ** Examples

lib <- coelho_stellar_subset()
dim(lib$matrix)
table(lib$metadata$type)



cleanEx()
nameEx("component_map")
### * component_map

flush(stderr()); flush(stdout())

### Name: component_map
### Title: Extract a Component Map
### Aliases: component_map

### ** Examples

cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
fit <- spectral_unmix(cube$matrix, k = 3, niter = 50)
img <- component_map(fit, cube$nx, cube$ny, 2)
dim(img)



cleanEx()
nameEx("component_reconstruction")
### * component_reconstruction

flush(stderr()); flush(stdout())

### Name: component_reconstruction
### Title: Reconstruct a Single Component
### Aliases: component_reconstruction

### ** Examples

## Not run: 
##D cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
##D fit <- spectral_unmix(cube$matrix, k = 3, niter = 50)
##D comp1 <- component_reconstruction(fit, 1)
##D comp1_cube <- component_reconstruction(fit, 1, nx = cube$nx, ny = cube$ny)
## End(Not run)



cleanEx()
nameEx("component_spectrum")
### * component_spectrum

flush(stderr()); flush(stdout())

### Name: component_spectrum
### Title: Extract a Component Spectrum
### Aliases: component_spectrum

### ** Examples

## Not run: 
##D cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
##D fit <- spectral_unmix(cube$matrix, k = 3, niter = 50)
##D spec1 <- component_spectrum(fit, 1)
## End(Not run)



cleanEx()
nameEx("cube_helpers")
### * cube_helpers

flush(stderr()); flush(stdout())

### Name: cube_to_matrix
### Title: Helpers for Cube Reshaping and Map Extraction
### Aliases: cube_to_matrix matrix_to_cube component_map component_spectrum
###   component_reconstruction cube_metadata

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
##D cube_metadata(comp2)
## End(Not run)



cleanEx()
nameEx("cube_metadata")
### * cube_metadata

flush(stderr()); flush(stdout())

### Name: cube_metadata
### Title: Return Stored Metadata
### Aliases: cube_metadata

### ** Examples

demo <- simulate_ifu_cube()
cube_metadata(demo$cube)



cleanEx()
nameEx("cube_to_matrix")
### * cube_to_matrix

flush(stderr()); flush(stdout())

### Name: cube_to_matrix
### Title: Reshape a Spectral Cube Into a Matrix
### Aliases: cube_to_matrix

### ** Examples

cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
x <- cube_to_matrix(cube$cube)
dim(x)



cleanEx()
nameEx("demo_ifu_unmix")
### * demo_ifu_unmix

flush(stderr()); flush(stdout())

### Name: demo_ifu_unmix
### Title: Run the Built-in IFU Demonstration
### Aliases: demo_ifu_unmix

### ** Examples

## Not run: 
##D demo <- demo_ifu_unmix(plot = TRUE, niter = 150)
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
nameEx("fitted.spectral_unmix")
### * fitted.spectral_unmix

flush(stderr()); flush(stdout())

### Name: fitted.spectral_unmix
### Title: Return the Fitted Reconstruction
### Aliases: fitted.spectral_unmix

### ** Examples

fake <- list(reconstruction = matrix(runif(20), 4, 5))
class(fake) <- "spectral_unmix"
fitted(fake)



cleanEx()
nameEx("matrix_to_cube")
### * matrix_to_cube

flush(stderr()); flush(stdout())

### Name: matrix_to_cube
### Title: Reshape a Matrix Back Into a Cube
### Aliases: matrix_to_cube

### ** Examples

cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
rebuilt <- matrix_to_cube(cube$matrix, cube$nx, cube$ny)
dim(rebuilt)



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
nameEx("plot.spectral_unmix")
### * plot.spectral_unmix

flush(stderr()); flush(stdout())

### Name: plot.spectral_unmix
### Title: Plot a Spectral Unmixing Result
### Aliases: plot.spectral_unmix

### ** Examples

## Not run: 
##D cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
##D fit <- spectral_unmix(cube$matrix, k = 3, niter = 50)
##D plot(fit, type = "spectra", wavelength = cube$wavelength)
##D plot(fit, type = "maps", nx = cube$nx, ny = cube$ny)
##D plot(fit, type = "loss")
## End(Not run)



cleanEx()
nameEx("plot_reconstruction")
### * plot_reconstruction

flush(stderr()); flush(stdout())

### Name: plot_reconstruction
### Title: Plot Data Versus Reconstruction
### Aliases: plot_reconstruction

### ** Examples

## Not run: 
##D cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
##D fit <- spectral_unmix(cube$matrix, k = 3, niter = 50)
##D plot_reconstruction(fit, cube$matrix, n = 4, wavelength = cube$wavelength)
## End(Not run)



cleanEx()
nameEx("predict.spectral_unmix")
### * predict.spectral_unmix

flush(stderr()); flush(stdout())

### Name: predict.spectral_unmix
### Title: Predict From a Fitted Spectral Unmixing Model
### Aliases: predict.spectral_unmix

### ** Examples

## Not run: 
##D demo <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
##D fit <- spectral_unmix(demo$matrix, k = 3, niter = 50)
##D fitted_x <- predict(fit)
##D new_weights <- predict(fit, newdata = demo$cube, type = "spatial")
## End(Not run)



cleanEx()
nameEx("residuals.spectral_unmix")
### * residuals.spectral_unmix

flush(stderr()); flush(stdout())

### Name: residuals.spectral_unmix
### Title: Compute Residuals
### Aliases: residuals.spectral_unmix

### ** Examples

fake <- list(reconstruction = matrix(runif(20), 4, 5))
class(fake) <- "spectral_unmix"
x <- matrix(runif(20), 4, 5)
residuals(fake, x = x)



cleanEx()
nameEx("simulate_ifu_cube")
### * simulate_ifu_cube

flush(stderr()); flush(stdout())

### Name: simulate_ifu_cube
### Title: Simulate a Small IFU-like Cube
### Aliases: simulate_ifu_cube

### ** Examples

demo <- simulate_ifu_cube()
str(demo, max.level = 1)



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
nameEx("summary.spectral_unmix")
### * summary.spectral_unmix

flush(stderr()); flush(stdout())

### Name: summary.spectral_unmix
### Title: Summarize a Spectral Unmixing Fit
### Aliases: summary.spectral_unmix

### ** Examples

fake <- list(
  spatial = matrix(runif(12), 4, 3),
  spectra = matrix(runif(15), 3, 5),
  loss = c(3, 2, 1)
)
class(fake) <- "spectral_unmix"
summary(fake)



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
