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
### Aliases: cube_to_matrix matrix_to_cube component_map

### ** Examples

demo <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
x <- cube_to_matrix(demo$cube)
rebuilt <- matrix_to_cube(x, demo$nx, demo$ny)
dim(rebuilt)
## Not run: 
##D fit <- spectral_unmix(x, k = 3, niter = 50)
##D img <- component_map(fit, demo$nx, demo$ny, 2)
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
