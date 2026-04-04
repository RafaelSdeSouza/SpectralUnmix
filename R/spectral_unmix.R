#' Spectral Unmixing for Hyperspectral and IFU Data
#'
#' Fits a non-negative matrix factorization model with an optional smoothness
#' penalty on the recovered spectra.
#'
#' The factorization follows the linear mixing model
#' \deqn{X \approx A S}
#' where rows of \eqn{X} correspond to spectra from individual spatial elements,
#' \eqn{A} contains non-negative abundance weights, and \eqn{S} contains
#' non-negative spectral components.
#'
#' @param x Numeric matrix-like object with spectra in rows and wavelength
#'   channels in columns.
#' @param k Number of components to estimate.
#' @param lambda_smooth Non-negative weight for the first-difference smoothness
#'   penalty applied to the recovered spectra.
#' @param init Initialization scheme for the non-negative factors.
#'   \code{"nndsvd"} usually gives a more stable starting point than
#'   \code{"random"}.
#' @param nstart Number of optimization restarts. The best solution by final
#'   objective value is returned.
#' @param seed Optional integer seed used to make restart generation
#'   reproducible.
#' @param lr Learning rate passed to \code{torch::optim_adam()}.
#' @param niter Number of optimization iterations.
#' @param tol Relative objective-improvement tolerance used for early stopping.
#' @param patience Number of consecutive iterations with improvement below
#'   \code{tol} before optimization stops early.
#' @param center Logical; if \code{TRUE}, center wavelength channels before
#'   fitting.
#' @param scale Logical; if \code{TRUE}, scale wavelength channels before
#'   fitting.
#' @param cuda Logical; if \code{TRUE}, request CUDA execution.
#' @param verbose Logical; if \code{TRUE}, emit progress every 100 iterations.
#' @param metadata Optional metadata list to carry through the fitted object,
#'   reconstructed cubes, and prediction outputs.
#'
#' @return A list with class \code{"spectral_unmix"} containing:
#' \itemize{
#'   \item \code{spatial}: abundance matrix with one column per component.
#'   \item \code{abundance}: alias of \code{spatial}.
#'   \item \code{spectra}: component spectra with one row per component.
#'   \item \code{basis}: component spectra arranged as wavelength-by-component.
#'   \item \code{coef}: component weights arranged as component-by-spaxel.
#'   \item \code{reconstruction}: reconstructed matrix \eqn{AS}.
#'   \item \code{fitted}: alias of \code{reconstruction}.
#'   \item \code{loss}: objective history over the optimization.
#'   \item \code{niter_run}: number of iterations actually executed.
#'   \item \code{converged}: logical flag indicating whether early-stopping
#'     convergence was detected before \code{niter}.
#'   \item \code{metadata}: optional metadata carried from the input matrix or
#'     provided explicitly.
#'   \item \code{center}: centering values used during preprocessing, or
#'     \code{FALSE}.
#'   \item \code{scale}: scaling values used during preprocessing, or
#'     \code{FALSE}.
#'   \item \code{call}: matched function call.
#' }
#'
#' @details
#' This implementation is aimed at hyperspectral cubes and integral-field
#' spectroscopy data after reshaping the spatial dimensions into a two-dimensional
#' matrix of spaxels by wavelength. The objective function is
#' \deqn{\|X - AS\|^2 + \lambda \|D S\|^2}
#' where the second term penalizes rapid channel-to-channel variations in each
#' component spectrum.
#'
#' For astronomy workflows, a typical sequence is:
#' \enumerate{
#'   \item reshape a cube with [cube_to_matrix()];
#'   \item fit the model with \code{spectral_unmix()};
#'   \item recover component maps with [component_map()].
#' }
#'
#' @references
#' Lee, D. D., and Seung, H. S. (1999). Learning the parts of objects by
#' non-negative matrix factorization. \emph{Nature}, 401, 788-791.
#'
#' @examples
#' cube <- simulate_ifu_cube(nx = 12, ny = 12, n_wave = 80)
#' fit <- spectral_unmix(cube$matrix, k = 3, niter = 150, lr = 0.03)
#' print(fit)
#'
#' map1 <- component_map(fit, nx = cube$nx, ny = cube$ny, component = 1)
#' dim(map1)
#'
#' @export
spectral_unmix <- function(x,
                           k = 3,
                           lambda_smooth = 0.01,
                           init = c("nndsvd", "random"),
                           nstart = 1,
                           seed = NULL,
                           lr = 0.02,
                           niter = 2000,
                           tol = 1e-6,
                           patience = 20,
                           center = FALSE,
                           scale = FALSE,
                           cuda = FALSE,
                           verbose = FALSE,
                           metadata = NULL) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("The 'torch' package must be installed to use spectral_unmix().", call. = FALSE)
  }

  input_metadata <- merge_metadata(extract_matrix_metadata(x), metadata)
  x <- validate_input_matrix(x)
  init <- match.arg(init)

  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 1L) {
    stop("'k' must be a positive integer.", call. = FALSE)
  }
  k <- as.integer(k)

  if (k > min(dim(x))) {
    stop("'k' must not exceed the smaller input dimension.", call. = FALSE)
  }
  if (!is.numeric(lambda_smooth) || length(lambda_smooth) != 1L || is.na(lambda_smooth) ||
      lambda_smooth < 0) {
    stop("'lambda_smooth' must be a non-negative number.", call. = FALSE)
  }
  if (!is.numeric(nstart) || length(nstart) != 1L || is.na(nstart) || nstart < 1) {
    stop("'nstart' must be a positive integer.", call. = FALSE)
  }
  nstart <- as.integer(nstart)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed))) {
    stop("'seed' must be NULL or a single integer.", call. = FALSE)
  }
  if (!is.numeric(lr) || length(lr) != 1L || is.na(lr) || lr <= 0) {
    stop("'lr' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(niter) || length(niter) != 1L || is.na(niter) || niter < 1) {
    stop("'niter' must be a positive integer.", call. = FALSE)
  }
  niter <- as.integer(niter)
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0) {
    stop("'tol' must be a non-negative number.", call. = FALSE)
  }
  if (!is.numeric(patience) || length(patience) != 1L || is.na(patience) || patience < 1) {
    stop("'patience' must be a positive integer.", call. = FALSE)
  }
  patience <- as.integer(patience)

  x <- base::scale(x, center = center, scale = scale)
  cen <- attr(x, "scaled:center")
  sc <- attr(x, "scaled:scale")

  if (any(!is.finite(x))) {
    stop("Preprocessing produced non-finite values.", call. = FALSE)
  }
  if (any(x < 0)) {
    warning(
      "Input contains negative values after preprocessing. Interpretability of ",
      "non-negative components may be reduced.",
      call. = FALSE
    )
  }

  device <- select_torch_device(cuda)
  fit_core <- fit_spectral_unmix_with_restarts(
    x = x,
    k = k,
    init = init,
    nstart = nstart,
    seed = seed,
    lr = lr,
    niter = niter,
    tol = tol,
    patience = patience,
    lambda_smooth = lambda_smooth,
    device = device,
    verbose = verbose
  )

  spatial <- fit_core$spatial
  spectra <- fit_core$spectra
  reconstruction <- fit_core$reconstruction

  fit <- list(
    spatial = spatial,
    abundance = spatial,
    spectra = spectra,
    basis = t(spectra),
    coef = t(spatial),
    reconstruction = reconstruction,
    fitted = reconstruction,
    loss = fit_core$loss,
    init = init,
    nstart = nstart,
    seed = if (is.null(seed)) NULL else as.integer(seed),
    restart_seed = fit_core$restart_seed,
    restart_losses = fit_core$restart_losses,
    niter_run = fit_core$niter_run,
    converged = fit_core$converged,
    metadata = input_metadata,
    center = if (is.null(cen)) FALSE else cen,
    scale = if (is.null(sc)) FALSE else sc,
    call = match.call()
  )

  class(fit) <- "spectral_unmix"
  fit
}

#' Reshape a Spectral Cube Into a Matrix
#'
#' Converts a three-dimensional array with dimensions \code{nx x ny x n_wave}
#' into a matrix suitable for [spectral_unmix()], with one row per spaxel.
#'
#' @param cube Numeric 3D array, or a list-like object containing a cube in
#'   \code{$imDat} or \code{$cube}.
#'
#' @return A numeric matrix with \code{nx * ny} rows and \code{n_wave} columns.
#' @examples
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' x <- cube_to_matrix(cube$cube)
#' dim(x)
#' @export
cube_to_matrix <- function(cube) {
  metadata <- extract_cube_metadata(cube)
  cube <- resolve_cube(cube)

  if (!is.array(cube) || length(dim(cube)) != 3L) {
    stop("'cube' must be a numeric 3D array.", call. = FALSE)
  }
  if (!is.numeric(cube) || any(!is.finite(cube))) {
    stop("'cube' must contain finite numeric values.", call. = FALSE)
  }

  dims <- dim(cube)
  x <- matrix(cube, nrow = dims[1L] * dims[2L], ncol = dims[3L], byrow = FALSE)
  attr(x, "spectral_unmix_metadata") <- metadata
  attr(x, "spectral_unmix_dim") <- dims[1:2]
  x
}

#' Reshape a Matrix Back Into a Cube
#'
#' Reconstructs a \code{nx x ny x n_wave} cube from a matrix with one row per
#' spaxel.
#'
#' @param x Numeric matrix with one row per spatial pixel.
#' @param nx Number of x-axis pixels.
#' @param ny Number of y-axis pixels.
#' @param metadata Optional metadata list to attach to the returned cube.
#'
#' @return A numeric 3D array.
#' @examples
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' rebuilt <- matrix_to_cube(cube$matrix, cube$nx, cube$ny)
#' dim(rebuilt)
#' @export
matrix_to_cube <- function(x, nx, ny, metadata = NULL) {
  x <- validate_input_matrix(x)

  if (!is.numeric(nx) || length(nx) != 1L || nx < 1L ||
      !is.numeric(ny) || length(ny) != 1L || ny < 1L) {
    stop("'nx' and 'ny' must be positive integers.", call. = FALSE)
  }

  nx <- as.integer(nx)
  ny <- as.integer(ny)

  if (nrow(x) != nx * ny) {
    stop("nrow(x) must equal nx * ny.", call. = FALSE)
  }

  cube <- array(x, dim = c(nx, ny, ncol(x)))
  cube_metadata <- merge_metadata(extract_matrix_metadata(x), metadata)
  if (length(cube_metadata) > 0L) {
    attr(cube, "spectral_unmix_metadata") <- cube_metadata
  }
  cube
}

#' Extract a Component Map
#'
#' Reshapes one abundance column from a fitted model into an image plane.
#'
#' @param fit Result from [spectral_unmix()].
#' @param nx Number of x-axis pixels.
#' @param ny Number of y-axis pixels.
#' @param component Component index to extract.
#'
#' @return A numeric matrix of size \code{nx x ny}.
#' @examples
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spectral_unmix(cube$matrix, k = 3, niter = 50)
#' img <- component_map(fit, cube$nx, cube$ny, 2)
#' dim(img)
#' @export
component_map <- function(fit, nx, ny, component = 1) {
  if (!inherits(fit, "spectral_unmix")) {
    stop("'fit' must be a result from spectral_unmix().", call. = FALSE)
  }
  if (!is.numeric(nx) || length(nx) != 1L || nx < 1L ||
      !is.numeric(ny) || length(ny) != 1L || ny < 1L) {
    stop("'nx' and 'ny' must be positive integers.", call. = FALSE)
  }
  if (!is.numeric(component) || length(component) != 1L || component < 1L ||
      component > ncol(fit$spatial)) {
    stop("'component' is out of bounds.", call. = FALSE)
  }

  component <- as.integer(component)
  matrix(fit$spatial[, component], nrow = as.integer(nx), ncol = as.integer(ny))
}

#' Extract a Component Spectrum
#'
#' Returns one recovered component spectrum from a fitted model.
#'
#' @param fit Result from [spectral_unmix()].
#' @param component Component index to extract.
#'
#' @return A numeric vector.
#' @examples
#' \dontrun{
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spectral_unmix(cube$matrix, k = 3, niter = 50)
#' spec1 <- component_spectrum(fit, 1)
#' }
#' @export
component_spectrum <- function(fit, component = 1) {
  if (!inherits(fit, "spectral_unmix")) {
    stop("'fit' must be a result from spectral_unmix().", call. = FALSE)
  }
  if (!is.numeric(component) || length(component) != 1L || component < 1L ||
      component > nrow(fit$spectra)) {
    stop("'component' is out of bounds.", call. = FALSE)
  }

  fit$spectra[as.integer(component), ]
}

#' Reconstruct a Single Component
#'
#' Builds the contribution from one component, either as a spaxel-by-wavelength
#' matrix or, when dimensions are supplied, as a spectral cube.
#'
#' @param fit Result from [spectral_unmix()].
#' @param component Component index to reconstruct.
#' @param nx Optional x-axis size for cube output.
#' @param ny Optional y-axis size for cube output.
#'
#' @return A numeric matrix, or a numeric 3D array when \code{nx} and
#'   \code{ny} are supplied.
#' @examples
#' \dontrun{
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spectral_unmix(cube$matrix, k = 3, niter = 50)
#' comp1 <- component_reconstruction(fit, 1)
#' comp1_cube <- component_reconstruction(fit, 1, nx = cube$nx, ny = cube$ny)
#' }
#' @export
component_reconstruction <- function(fit, component = 1, nx = NULL, ny = NULL) {
  if (!inherits(fit, "spectral_unmix")) {
    stop("'fit' must be a result from spectral_unmix().", call. = FALSE)
  }
  if (!is.numeric(component) || length(component) != 1L || component < 1L ||
      component > nrow(fit$spectra)) {
    stop("'component' is out of bounds.", call. = FALSE)
  }

  component <- as.integer(component)
  reconstruction <- fit$spatial[, component, drop = FALSE] %*%
    fit$spectra[component, , drop = FALSE]

  if (is.null(nx) && is.null(ny)) {
    return(reconstruction)
  }
  if (is.null(nx) || is.null(ny)) {
    stop("Supply both 'nx' and 'ny' to return a cube.", call. = FALSE)
  }

  matrix_to_cube(reconstruction, nx = nx, ny = ny, metadata = fit$metadata)
}

#' Simulate a Small IFU-like Cube
#'
#' Generates a synthetic cube with smooth continuum structure and simple
#' emission features. This is intended for examples, demos, and documentation.
#'
#' @param nx Number of x-axis pixels.
#' @param ny Number of y-axis pixels.
#' @param n_wave Number of wavelength channels.
#' @param noise Standard deviation of additive Gaussian noise.
#'
#' @return A list with the simulated \code{cube}, its matrix form, the true
#' component \code{spectra}, abundance maps \code{abundances}, wavelength grid,
#' and spatial dimensions.
#' @examples
#' demo <- simulate_ifu_cube()
#' str(demo, max.level = 1)
#' @export
simulate_ifu_cube <- function(nx = 20, ny = 20, n_wave = 120, noise = 0.01) {
  if (!is.numeric(nx) || !is.numeric(ny) || !is.numeric(n_wave) || !is.numeric(noise)) {
    stop("All arguments must be numeric.", call. = FALSE)
  }

  nx <- as.integer(nx)
  ny <- as.integer(ny)
  n_wave <- as.integer(n_wave)

  wavelength <- seq(3600, 7200, length.out = n_wave)
  xgrid <- seq(-1, 1, length.out = nx)
  ygrid <- seq(-1, 1, length.out = ny)
  xy <- expand.grid(x = xgrid, y = ygrid)

  g1 <- exp(-((xy$x + 0.35)^2 + (xy$y + 0.1)^2) / 0.10)
  g2 <- exp(-((xy$x - 0.25)^2 + (xy$y - 0.3)^2) / 0.18)
  disk <- exp(-sqrt(xy$x^2 + xy$y^2) / 0.9)

  spectra <- rbind(
    0.8 + 0.6 * (wavelength / max(wavelength))^-1.2,
    0.15 + 0.25 * exp(-0.5 * ((wavelength - 5007) / 40)^2) +
      0.18 * exp(-0.5 * ((wavelength - 6563) / 55)^2),
    0.25 + 0.35 * exp(-0.5 * ((wavelength - 4300) / 180)^2)
  )

  abundances <- cbind(
    0.6 * disk + 0.25 * g1,
    0.8 * g1,
    0.7 * g2
  )

  matrix_data <- abundances %*% spectra
  matrix_data <- matrix_data + stats::rnorm(length(matrix_data), sd = noise)
  matrix_data[matrix_data < 0] <- 0

  list(
    cube = matrix_to_cube(matrix_data, nx = nx, ny = ny),
    matrix = matrix_data,
    spectra = spectra,
    abundances = lapply(seq_len(ncol(abundances)), function(i) {
      matrix(abundances[, i], nrow = nx, ncol = ny)
    }),
    wavelength = wavelength,
    nx = nx,
    ny = ny
  )
}

#' Run the Built-in IFU Demonstration
#'
#' Generates a synthetic cube, fits the model, and optionally plots the result.
#'
#' @param k Number of components.
#' @param nx Number of x-axis pixels in the synthetic cube.
#' @param ny Number of y-axis pixels in the synthetic cube.
#' @param n_wave Number of wavelength channels.
#' @param niter Number of optimization iterations.
#' @param plot Logical; if \code{TRUE}, display a compact base-R summary plot.
#'
#' @return A list with the simulated data and fitted model.
#' @examples
#' \dontrun{
#' demo <- demo_ifu_unmix(plot = TRUE, niter = 150)
#' }
#' @export
demo_ifu_unmix <- function(k = 3, nx = 18, ny = 18, n_wave = 120, niter = 300, plot = interactive()) {
  demo <- simulate_ifu_cube(nx = nx, ny = ny, n_wave = n_wave)
  fit <- spectral_unmix(demo$matrix, k = k, niter = niter)

  if (isTRUE(plot)) {
    plot(fit, type = "spectra", wavelength = demo$wavelength)
    plot(fit, type = "maps", nx = demo$nx, ny = demo$ny)
  }

  list(data = demo, fit = fit)
}

#' @export
print.spectral_unmix <- function(x, ...) {
  niter_run <- if (!is.null(x$niter_run)) x$niter_run else length(x$loss)
  converged <- isTRUE(x$converged)

  cat("<spectral_unmix>\n", sep = "")
  cat(sprintf("  spaxels: %d\n", nrow(x$spatial)))
  cat(sprintf("  wavelength channels: %d\n", ncol(x$spectra)))
  cat(sprintf("  components: %d\n", nrow(x$spectra)))
  if (!is.null(x$init)) {
    cat(sprintf("  init: %s\n", x$init))
  }
  if (!is.null(x$nstart)) {
    cat(sprintf("  restarts: %d\n", x$nstart))
  }
  cat(sprintf("  iterations: %d\n", niter_run))
  cat(sprintf("  converged: %s\n", if (converged) "yes" else "no"))
  cat(sprintf("  final loss: %.6f\n", utils::tail(x$loss, 1L)))
  invisible(x)
}

#' Return Stored Metadata
#'
#' Retrieves metadata carried by matrices, cubes, or fitted
#' \code{"spectral_unmix"} objects.
#'
#' @param object Matrix, cube, or fitted object.
#'
#' @return A metadata list, or \code{NULL} when no metadata is present.
#' @examples
#' demo <- simulate_ifu_cube()
#' cube_metadata(demo$cube)
#' @export
cube_metadata <- function(object) {
  if (inherits(object, "spectral_unmix")) {
    return(object$metadata)
  }

  attr(object, "spectral_unmix_metadata", exact = TRUE)
}

#' Summarize a Spectral Unmixing Fit
#'
#' @param object Result from [spectral_unmix()].
#' @param ... Unused.
#'
#' @return A list with class \code{"summary.spectral_unmix"}.
#' @examples
#' fake <- list(
#'   spatial = matrix(runif(12), 4, 3),
#'   spectra = matrix(runif(15), 3, 5),
#'   loss = c(3, 2, 1)
#' )
#' class(fake) <- "spectral_unmix"
#' summary(fake)
#' @export
summary.spectral_unmix <- function(object, ...) {
  if (!inherits(object, "spectral_unmix")) {
    stop("'object' must be a result from spectral_unmix().", call. = FALSE)
  }

  summary_object <- list(
    n_spaxels = nrow(object$spatial),
    n_wave = ncol(object$spectra),
    rank = nrow(object$spectra),
    final_loss = utils::tail(object$loss, 1L),
    n_iter = if (!is.null(object$niter_run)) object$niter_run else length(object$loss),
    init = object$init,
    nstart = object$nstart,
    converged = isTRUE(object$converged)
  )

  class(summary_object) <- "summary.spectral_unmix"
  summary_object
}

#' @export
print.summary.spectral_unmix <- function(x, ...) {
  cat("SpectralUnmix summary\n")
  cat(sprintf("  spaxels: %d\n", x$n_spaxels))
  cat(sprintf("  wavelength channels: %d\n", x$n_wave))
  cat(sprintf("  components: %d\n", x$rank))
  if (!is.null(x$init)) {
    cat(sprintf("  init: %s\n", x$init))
  }
  if (!is.null(x$nstart)) {
    cat(sprintf("  restarts: %d\n", x$nstart))
  }
  cat(sprintf("  iterations: %d\n", x$n_iter))
  cat(sprintf("  converged: %s\n", if (isTRUE(x$converged)) "yes" else "no"))
  cat(sprintf("  final loss: %.6f\n", x$final_loss))
  invisible(x)
}

#' NMF-style Basis Matrix
#'
#' Returns the component spectra arranged as wavelength-by-component, similar to
#' a basis matrix in NMF software interfaces.
#'
#' @param object Result from [spectral_unmix()].
#' @param ... Unused.
#'
#' @return A numeric matrix with wavelengths in rows and components in columns.
#' @examples
#' fake <- list(spectra = matrix(runif(15), 3, 5))
#' class(fake) <- "spectral_unmix"
#' basis(fake)
#' @export
basis <- function(object, ...) {
  UseMethod("basis")
}

#' @export
basis.spectral_unmix <- function(object, ...) {
  if (!inherits(object, "spectral_unmix")) {
    stop("'object' must be a result from spectral_unmix().", call. = FALSE)
  }

  t(object$spectra)
}

#' Extract Abundance Coefficients
#'
#' Returns the component weights arranged as component-by-spaxel, similar to the
#' coefficient matrix used by NMF interfaces.
#'
#' @param object Result from [spectral_unmix()].
#' @param ... Unused.
#'
#' @return A numeric matrix with components in rows and spaxels in columns.
#' @examples
#' fake <- list(spatial = matrix(runif(12), 4, 3))
#' class(fake) <- "spectral_unmix"
#' coef(fake)
#' @export
coef.spectral_unmix <- function(object, ...) {
  if (!inherits(object, "spectral_unmix")) {
    stop("'object' must be a result from spectral_unmix().", call. = FALSE)
  }

  t(object$spatial)
}

#' Return the Fitted Reconstruction
#'
#' @param object Result from [spectral_unmix()].
#' @param nx Optional x-axis size for cube output.
#' @param ny Optional y-axis size for cube output.
#' @param ... Unused.
#'
#' @return A numeric matrix, or a numeric cube when \code{nx} and \code{ny} are
#' supplied.
#' @examples
#' fake <- list(reconstruction = matrix(runif(20), 4, 5))
#' class(fake) <- "spectral_unmix"
#' fitted(fake)
#' @export
fitted.spectral_unmix <- function(object, nx = NULL, ny = NULL, ...) {
  if (!inherits(object, "spectral_unmix")) {
    stop("'object' must be a result from spectral_unmix().", call. = FALSE)
  }

  if (is.null(nx) && is.null(ny)) {
    return(object$reconstruction)
  }
  if (is.null(nx) || is.null(ny)) {
    stop("Supply both 'nx' and 'ny' to return a cube.", call. = FALSE)
  }

  matrix_to_cube(object$reconstruction, nx = nx, ny = ny, metadata = object$metadata)
}

#' Compute Residuals
#'
#' @param object Result from [spectral_unmix()].
#' @param x Numeric matrix originally used for the fit.
#' @param nx Optional x-axis size for cube output.
#' @param ny Optional y-axis size for cube output.
#' @param ... Unused.
#'
#' @return A numeric matrix, or a numeric cube when \code{nx} and \code{ny} are
#' supplied.
#' @examples
#' fake <- list(reconstruction = matrix(runif(20), 4, 5))
#' class(fake) <- "spectral_unmix"
#' x <- matrix(runif(20), 4, 5)
#' residuals(fake, x = x)
#' @export
residuals.spectral_unmix <- function(object, x, nx = NULL, ny = NULL, ...) {
  if (!inherits(object, "spectral_unmix")) {
    stop("'object' must be a result from spectral_unmix().", call. = FALSE)
  }

  x <- validate_input_matrix(x)
  if (nrow(x) != nrow(object$reconstruction) || ncol(x) != ncol(object$reconstruction)) {
    stop("'x' must have the same dimensions as object$reconstruction.", call. = FALSE)
  }
  residual <- x - object$reconstruction

  if (is.null(nx) && is.null(ny)) {
    return(residual)
  }
  if (is.null(nx) || is.null(ny)) {
    stop("Supply both 'nx' and 'ny' to return a cube.", call. = FALSE)
  }

  matrix_to_cube(residual, nx = nx, ny = ny, metadata = object$metadata)
}

#' Predict From a Fitted Spectral Unmixing Model
#'
#' Uses the fitted component spectra to estimate abundance weights for new data,
#' or returns the in-sample reconstruction when \code{newdata} is omitted.
#'
#' @param object Result from [spectral_unmix()].
#' @param newdata Optional new matrix or cube.
#' @param type Prediction target: \code{"reconstruction"}, \code{"spatial"},
#'   or \code{"cube"}.
#' @param nx Optional x-axis size for cube output.
#' @param ny Optional y-axis size for cube output.
#' @param lr Learning rate used when estimating abundance weights for new data.
#' @param niter Number of optimization iterations for new data.
#' @param tol Relative objective-improvement tolerance used for early stopping
#'   during prediction.
#' @param patience Number of consecutive iterations with improvement below
#'   \code{tol} before prediction stops early.
#' @param init Initialization scheme for abundance optimization when predicting
#'   new data. \code{"least_squares"} uses a clipped least-squares warm start.
#' @param nstart Number of prediction restarts. The best solution by final
#'   objective value is returned.
#' @param seed Optional integer seed used to make prediction restarts
#'   reproducible.
#' @param cuda Logical; if \code{TRUE}, request CUDA execution.
#' @param verbose Logical; if \code{TRUE}, emit progress every 100 iterations.
#' @param ... Unused.
#'
#' @return A numeric matrix or cube, depending on \code{type}.
#' @examples
#' \dontrun{
#' demo <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spectral_unmix(demo$matrix, k = 3, niter = 50)
#' fitted_x <- predict(fit)
#' new_weights <- predict(fit, newdata = demo$cube, type = "spatial")
#' }
#' @export
predict.spectral_unmix <- function(object,
                                   newdata = NULL,
                                   type = c("reconstruction", "spatial", "cube"),
                                   nx = NULL,
                                   ny = NULL,
                                   lr = 0.05,
                                   niter = 500,
                                   tol = 1e-6,
                                   patience = 20,
                                   init = c("least_squares", "random"),
                                   nstart = 1,
                                   seed = NULL,
                                   cuda = FALSE,
                                   verbose = FALSE,
                                   ...) {
  if (!inherits(object, "spectral_unmix")) {
    stop("'object' must be a result from spectral_unmix().", call. = FALSE)
  }

  type <- match.arg(type)
  init <- match.arg(init)

  if (is.null(newdata)) {
    if (type == "spatial") {
      return(object$spatial)
    }
    if (type == "reconstruction") {
      return(object$reconstruction)
    }
    if (is.null(nx) || is.null(ny)) {
      stop("'nx' and 'ny' are required when type = 'cube'.", call. = FALSE)
    }
    return(matrix_to_cube(object$reconstruction, nx = nx, ny = ny, metadata = object$metadata))
  }

  prediction_metadata <- extract_matrix_metadata(newdata)
  cube_dims <- infer_cube_dims(newdata)
  if (!is.null(cube_dims)) {
    nx <- cube_dims[1L]
    ny <- cube_dims[2L]
    newdata <- cube_to_matrix(newdata)
    prediction_metadata <- merge_metadata(prediction_metadata, extract_matrix_metadata(newdata))
  }

  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("The 'torch' package must be installed to use predict().", call. = FALSE)
  }

  newdata <- validate_input_matrix(newdata)
  if (ncol(newdata) != ncol(object$spectra)) {
    stop("newdata must have the same number of columns as the fitted spectra.", call. = FALSE)
  }
  if (!is.numeric(lr) || length(lr) != 1L || is.na(lr) || lr <= 0) {
    stop("'lr' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(niter) || length(niter) != 1L || is.na(niter) || niter < 1L) {
    stop("'niter' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(tol) || length(tol) != 1L || is.na(tol) || tol < 0) {
    stop("'tol' must be a non-negative number.", call. = FALSE)
  }
  if (!is.numeric(patience) || length(patience) != 1L || is.na(patience) || patience < 1L) {
    stop("'patience' must be a positive integer.", call. = FALSE)
  }
  patience <- as.integer(patience)
  if (!is.numeric(nstart) || length(nstart) != 1L || is.na(nstart) || nstart < 1L) {
    stop("'nstart' must be a positive integer.", call. = FALSE)
  }
  nstart <- as.integer(nstart)
  if (!is.null(seed) && (!is.numeric(seed) || length(seed) != 1L || is.na(seed))) {
    stop("'seed' must be NULL or a single integer.", call. = FALSE)
  }

  device <- select_torch_device(cuda)
  fit_pred <- fit_spectral_unmix_prediction_with_restarts(
    x = newdata,
    spectra = object$spectra,
    init = init,
    nstart = nstart,
    seed = seed,
    lr = lr,
    niter = as.integer(niter),
    tol = tol,
    patience = patience,
    device = device,
    verbose = verbose
  )

  spatial <- fit_pred$spatial
  reconstruction <- spatial %*% object$spectra

  if (type == "spatial") {
    return(spatial)
  }
  if (type == "reconstruction") {
    return(reconstruction)
  }

  if (is.null(nx) || is.null(ny)) {
    stop("'nx' and 'ny' are required when type = 'cube'.", call. = FALSE)
  }

  matrix_to_cube(
    reconstruction,
    nx = nx,
    ny = ny,
    metadata = merge_metadata(prediction_metadata, object$metadata)
  )
}

#' Plot a Spectral Unmixing Result
#'
#' Convenience plots for component spectra, spatial maps, and optimization loss.
#'
#' @param x Result from [spectral_unmix()].
#' @param type Plot type: \code{"spectra"}, \code{"maps"}, or \code{"loss"}.
#' @param nx Number of x-axis pixels when plotting maps.
#' @param ny Number of y-axis pixels when plotting maps.
#' @param components Optional vector of component indices to display.
#' @param wavelength Optional wavelength vector for spectra.
#' @param ... Additional graphical arguments passed to base plotting functions.
#'
#' @return Invisibly returns \code{x}.
#' @examples
#' \dontrun{
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spectral_unmix(cube$matrix, k = 3, niter = 50)
#' plot(fit, type = "spectra", wavelength = cube$wavelength)
#' plot(fit, type = "maps", nx = cube$nx, ny = cube$ny)
#' plot(fit, type = "loss")
#' }
#' @export
plot.spectral_unmix <- function(x,
                                type = c("spectra", "maps", "loss"),
                                nx = NULL,
                                ny = NULL,
                                components = NULL,
                                wavelength = NULL,
                                ...) {
  if (!inherits(x, "spectral_unmix")) {
    stop("'x' must be a result from spectral_unmix().", call. = FALSE)
  }

  type <- match.arg(type)

  if (is.null(components)) {
    components <- seq_len(nrow(x$spectra))
  }
  components <- as.integer(components)
  if (any(is.na(components)) || any(components < 1L) || any(components > nrow(x$spectra))) {
    stop("'components' contains invalid indices.", call. = FALSE)
  }

  if (type == "loss") {
    graphics::plot(
      seq_along(x$loss),
      x$loss,
      type = "l",
      xlab = "Iteration",
      ylab = "Loss",
      main = "Optimization history",
      ...
    )
    return(invisible(x))
  }

  if (type == "spectra") {
    if (is.null(wavelength)) {
      wavelength <- seq_len(ncol(x$spectra))
      xlab <- "Wavelength channel"
    } else {
      xlab <- "Wavelength"
    }

    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(length(components), 1L), mar = c(3, 4, 2, 1))

    for (component in components) {
      graphics::plot(
        wavelength,
        component_spectrum(x, component),
        type = "l",
        lwd = 2,
        xlab = xlab,
        ylab = "Flux",
        main = sprintf("Component spectrum %d", component),
        ...
      )
    }

    return(invisible(x))
  }

  if (is.null(nx) || is.null(ny)) {
    stop("'nx' and 'ny' are required when type = 'maps'.", call. = FALSE)
  }

  n_panels <- length(components)
  n_col <- min(3L, n_panels)
  n_row <- ceiling(n_panels / n_col)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(n_row, n_col), mar = c(3, 3, 2, 1))

  for (component in components) {
    graphics::image(
      component_map(x, nx = nx, ny = ny, component = component),
      xlab = "x",
      ylab = "y",
      main = sprintf("Component map %d", component),
      col = grDevices::hcl.colors(32, "YlOrRd", rev = TRUE),
      ...
    )
  }

  invisible(x)
}

#' Plot Data Versus Reconstruction
#'
#' Plots observed and reconstructed spectra for selected spaxels.
#'
#' @param fit Result from [spectral_unmix()].
#' @param x Numeric matrix used as input to the fit.
#' @param pixels Optional vector of spaxel indices to plot.
#' @param n Number of random spaxels to plot when \code{pixels} is omitted.
#' @param wavelength Optional wavelength vector.
#' @param ... Additional graphical arguments passed to \code{plot()}.
#'
#' @return Invisibly returns \code{fit}.
#' @examples
#' \dontrun{
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' fit <- spectral_unmix(cube$matrix, k = 3, niter = 50)
#' plot_reconstruction(fit, cube$matrix, n = 4, wavelength = cube$wavelength)
#' }
#' @export
plot_reconstruction <- function(fit, x, pixels = NULL, n = 6, wavelength = NULL, ...) {
  if (!inherits(fit, "spectral_unmix")) {
    stop("'fit' must be a result from spectral_unmix().", call. = FALSE)
  }

  x <- validate_input_matrix(x)

  if (nrow(x) != nrow(fit$reconstruction) || ncol(x) != ncol(fit$reconstruction)) {
    stop("'x' must have the same dimensions as fit$reconstruction.", call. = FALSE)
  }

  if (is.null(pixels)) {
    n <- min(as.integer(n), nrow(x))
    pixels <- sort(sample.int(nrow(x), n))
  } else {
    pixels <- as.integer(pixels)
  }
  if (any(is.na(pixels)) || any(pixels < 1L) || any(pixels > nrow(x))) {
    stop("'pixels' contains invalid indices.", call. = FALSE)
  }

  if (is.null(wavelength)) {
    wavelength <- seq_len(ncol(x))
    xlab <- "Wavelength channel"
  } else {
    xlab <- "Wavelength"
  }

  n_panels <- length(pixels)
  n_col <- min(2L, n_panels)
  n_row <- ceiling(n_panels / n_col)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(n_row, n_col), mar = c(3, 4, 2, 1))

  for (pixel in pixels) {
    graphics::plot(
      wavelength,
      x[pixel, ],
      type = "l",
      lwd = 2,
      col = "black",
      lty = 1,
      xlab = xlab,
      ylab = "Flux",
      main = sprintf("Spaxel %d", pixel),
      ...
    )
    graphics::lines(
      wavelength,
      fit$reconstruction[pixel, ],
      col = "red",
      lwd = 2,
      lty = 2
    )
    graphics::legend(
      "topright",
      legend = c("data", "reconstruction"),
      col = c("black", "red"),
      lty = c(1, 2),
      bty = "n"
    )
  }

  invisible(fit)
}

select_torch_device <- function(cuda) {
  if (!isTRUE(cuda)) {
    return(torch::torch_device("cpu"))
  }

  if (!torch::cuda_is_available()) {
    stop("CUDA was requested but torch reports that CUDA is not available.", call. = FALSE)
  }

  torch::torch_device("cuda:0")
}

validate_input_matrix <- function(x) {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("'x' must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (nrow(x) < 2L || ncol(x) < 2L) {
    stop("'x' must have at least two rows and two columns.", call. = FALSE)
  }
  if (any(!is.finite(x))) {
    stop("'x' must contain only finite numeric values.", call. = FALSE)
  }

  x
}

resolve_cube <- function(cube) {
  if (is.list(cube) && !is.null(cube$imDat)) {
    return(cube$imDat)
  }
  if (is.list(cube) && !is.null(cube$cube)) {
    return(cube$cube)
  }

  cube
}

infer_cube_dims <- function(x) {
  cube <- resolve_cube(x)
  if (!is.array(cube) || length(dim(cube)) != 3L) {
    return(NULL)
  }

  dim(cube)[1:2]
}

extract_cube_metadata <- function(x) {
  if (inherits(x, "spectral_unmix")) {
    return(x$metadata)
  }

  metadata <- attr(x, "spectral_unmix_metadata", exact = TRUE)
  if (!is.null(metadata)) {
    return(metadata)
  }

  if (is.list(x)) {
    keep <- setdiff(names(x), c("imDat", "cube"))
    if (length(keep) > 0L) {
      return(x[keep])
    }
  }

  NULL
}

set_optimization_seed <- function(seed) {
  if (is.null(seed)) {
    return(invisible(NULL))
  }

  set.seed(as.integer(seed))
  if (requireNamespace("torch", quietly = TRUE)) {
    torch::torch_manual_seed(as.integer(seed))
  }

  invisible(NULL)
}

generate_restart_seeds <- function(nstart, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  as.integer(sample.int(.Machine$integer.max, nstart))
}

random_nmf_factors <- function(x, k, seed = NULL, eps = 1e-6) {
  if (!is.null(seed)) {
    set.seed(as.integer(seed))
  }

  x_pos <- pmax(as.matrix(x), 0)
  scale_hint <- sqrt(max(mean(x_pos), eps) / max(k, 1L))
  hi <- max(scale_hint, eps * 10)

  list(
    spatial = matrix(stats::runif(nrow(x_pos) * k, min = eps, max = hi), nrow = nrow(x_pos), ncol = k),
    spectra = matrix(stats::runif(k * ncol(x_pos), min = eps, max = hi), nrow = k, ncol = ncol(x_pos))
  )
}

initialize_nmf_factors <- function(x, k, init = c("nndsvd", "random"), seed = NULL, eps = 1e-6) {
  init <- match.arg(init)

  if (init == "random") {
    return(random_nmf_factors(x, k, seed = seed, eps = eps))
  }

  x_pos <- pmax(as.matrix(x), 0)
  if (!any(x_pos > 0)) {
    return(random_nmf_factors(x, k, seed = seed, eps = eps))
  }

  rank_max <- min(dim(x_pos))
  svd_fit <- base::svd(x_pos, nu = min(rank_max, k), nv = min(rank_max, k))
  if (!length(svd_fit$d)) {
    return(random_nmf_factors(x, k, seed = seed, eps = eps))
  }

  W <- matrix(0, nrow = nrow(x_pos), ncol = k)
  H <- matrix(0, nrow = k, ncol = ncol(x_pos))

  u1 <- abs(svd_fit$u[, 1L])
  v1 <- abs(svd_fit$v[, 1L])
  W[, 1L] <- sqrt(svd_fit$d[1L]) * u1
  H[1L, ] <- sqrt(svd_fit$d[1L]) * v1

  rank_used <- min(k, length(svd_fit$d))
  if (rank_used >= 2L) {
    for (j in 2:rank_used) {
      uj <- svd_fit$u[, j]
      vj <- svd_fit$v[, j]
      up <- pmax(uj, 0)
      un <- pmax(-uj, 0)
      vp <- pmax(vj, 0)
      vn <- pmax(-vj, 0)

      up_norm <- sqrt(sum(up^2))
      un_norm <- sqrt(sum(un^2))
      vp_norm <- sqrt(sum(vp^2))
      vn_norm <- sqrt(sum(vn^2))

      pos_score <- up_norm * vp_norm
      neg_score <- un_norm * vn_norm

      if (pos_score >= neg_score && up_norm > 0 && vp_norm > 0) {
        uvec <- up / up_norm
        vvec <- vp / vp_norm
        sigma <- pos_score
      } else if (un_norm > 0 && vn_norm > 0) {
        uvec <- un / un_norm
        vvec <- vn / vn_norm
        sigma <- neg_score
      } else {
        next
      }

      coef <- sqrt(max(svd_fit$d[j] * sigma, 0))
      W[, j] <- coef * uvec
      H[j, ] <- coef * vvec
    }
  }

  fill <- max(mean(x_pos[x_pos > 0], na.rm = TRUE) * 1e-4, eps)
  if (!is.finite(fill) || fill <= 0) {
    fill <- eps
  }
  W[W <= 0] <- fill
  H[H <= 0] <- fill

  if (k > rank_used) {
    rand_fill <- random_nmf_factors(x, k - rank_used, seed = seed, eps = fill)
    W[, (rank_used + 1L):k] <- rand_fill$spatial
    H[(rank_used + 1L):k, ] <- rand_fill$spectra
  }

  list(spatial = W, spectra = H)
}

fit_spectral_unmix_with_restarts <- function(x,
                                             k,
                                             init = c("nndsvd", "random"),
                                             nstart = 1L,
                                             seed = NULL,
                                             ...) {
  init <- match.arg(init)
  restart_seeds <- generate_restart_seeds(nstart, seed = seed)
  restart_losses <- numeric(nstart)
  best_fit <- NULL
  best_loss <- Inf

  for (i in seq_len(nstart)) {
    init_factors <- initialize_nmf_factors(x, k, init = init, seed = restart_seeds[i])
    set_optimization_seed(restart_seeds[i])
    fit <- fit_spectral_unmix_adam(
      x = x,
      k = k,
      init_spatial = init_factors$spatial,
      init_spectra = init_factors$spectra,
      ...
    )

    restart_losses[i] <- utils::tail(fit$loss, 1L)
    if (restart_losses[i] < best_loss) {
      best_loss <- restart_losses[i]
      best_fit <- fit
      best_fit$restart_seed <- restart_seeds[i]
    }
  }

  best_fit$restart_losses <- restart_losses
  best_fit
}

fit_spectral_unmix_adam <- function(x,
                                    k,
                                    init_spatial,
                                    init_spectra,
                                    lr,
                                    niter,
                                    tol,
                                    patience,
                                    lambda_smooth,
                                    device,
                                    verbose = FALSE) {
  X <- torch::torch_tensor(as.matrix(x), device = device)
  n_wave <- ncol(x)
  A <- torch::torch_tensor(as.matrix(init_spatial), device = device, requires_grad = TRUE)
  S <- torch::torch_tensor(as.matrix(init_spectra), device = device, requires_grad = TRUE)
  opt <- torch::optim_adam(list(A, S), lr = lr)
  loss_history <- numeric(niter)
  stale_steps <- 0L
  converged <- FALSE
  best_loss <- Inf
  best_spatial <- NULL
  best_spectra <- NULL

  for (i in seq_len(niter)) {
    opt$zero_grad()

    Xhat <- A$matmul(S)
    recon <- torch::torch_mean((X - Xhat)^2)

    if (n_wave > 1L) {
      diff <- S[, 2:n_wave] - S[, 1:(n_wave - 1L)]
      smooth <- torch::torch_mean(diff^2)
    } else {
      smooth <- torch::torch_tensor(0, device = device)
    }

    loss <- recon + lambda_smooth * smooth
    loss$backward()
    opt$step()

    A$data()$clamp_(min = 0)
    S$data()$clamp_(min = 0)

    loss_history[i] <- as.numeric(loss$item())
    if (is.finite(loss_history[i]) && loss_history[i] < best_loss) {
      best_loss <- loss_history[i]
      best_spatial <- as.matrix(A$detach()$cpu())
      best_spectra <- as.matrix(S$detach()$cpu())
    }
    if (verbose && (i %% 100L == 0L || i == 1L || i == niter)) {
      message(sprintf("Iteration %d/%d, loss = %.6f", i, niter, loss_history[i]))
    }

    if (i > 1L) {
      rel_improvement <- (loss_history[i - 1L] - loss_history[i]) / max(loss_history[i - 1L], 1e-12)
      if (!is.finite(rel_improvement) || rel_improvement <= tol) {
        stale_steps <- stale_steps + 1L
      } else {
        stale_steps <- 0L
      }
      if (stale_steps >= patience) {
        converged <- TRUE
        loss_history <- loss_history[seq_len(i)]
        break
      }
    }
  }

  if (is.null(best_spatial) || is.null(best_spectra)) {
    spatial <- as.matrix(A$detach()$cpu())
    spectra <- as.matrix(S$detach()$cpu())
  } else {
    spatial <- best_spatial
    spectra <- best_spectra
  }
  reconstruction <- spatial %*% spectra

  list(
    spatial = spatial,
    spectra = spectra,
    reconstruction = reconstruction,
    loss = loss_history,
    niter_run = length(loss_history),
    converged = converged
  )
}

initialize_fixed_spectra_weights <- function(x,
                                             spectra,
                                             init = c("least_squares", "random"),
                                             seed = NULL,
                                             eps = 1e-8) {
  init <- match.arg(init)

  if (init == "random") {
    if (!is.null(seed)) {
      set.seed(as.integer(seed))
    }
    scale_hint <- sqrt(max(mean(pmax(x, 0)), eps) / max(nrow(spectra), 1L))
    hi <- max(scale_hint, eps * 10)
    return(matrix(
      stats::runif(nrow(x) * nrow(spectra), min = eps, max = hi),
      nrow = nrow(x),
      ncol = nrow(spectra)
    ))
  }

  gram <- spectra %*% t(spectra)
  gram <- gram + diag(eps, nrow(gram))
  init_weights <- x %*% t(spectra) %*% solve(gram)
  init_weights <- pmax(init_weights, eps)

  if (any(!is.finite(init_weights))) {
    return(initialize_fixed_spectra_weights(x, spectra, init = "random", seed = seed, eps = eps))
  }

  init_weights
}

fit_spectral_unmix_prediction_adam <- function(x,
                                               spectra,
                                               init_spatial,
                                               lr,
                                               niter,
                                               tol,
                                               patience,
                                               device,
                                               verbose = FALSE) {
  X <- torch::torch_tensor(as.matrix(x), device = device)
  S_fixed <- torch::torch_tensor(as.matrix(spectra), device = device)
  A <- torch::torch_tensor(as.matrix(init_spatial), device = device, requires_grad = TRUE)
  opt <- torch::optim_adam(list(A), lr = lr)
  loss_history <- numeric(niter)
  stale_steps <- 0L
  converged <- FALSE
  best_loss <- Inf
  best_spatial <- NULL

  for (i in seq_len(niter)) {
    opt$zero_grad()
    Xhat <- A$matmul(S_fixed)
    loss <- torch::torch_mean((X - Xhat)^2)
    loss$backward()
    opt$step()
    A$data()$clamp_(min = 0)

    loss_history[i] <- as.numeric(loss$item())
    if (is.finite(loss_history[i]) && loss_history[i] < best_loss) {
      best_loss <- loss_history[i]
      best_spatial <- as.matrix(A$detach()$cpu())
    }

    if (verbose && (i %% 100L == 0L || i == 1L || i == niter)) {
      message(sprintf("Prediction iteration %d/%d, loss = %.6f", i, niter, loss_history[i]))
    }

    if (i > 1L) {
      rel_improvement <- (loss_history[i - 1L] - loss_history[i]) / max(loss_history[i - 1L], 1e-12)
      if (!is.finite(rel_improvement) || rel_improvement <= tol) {
        stale_steps <- stale_steps + 1L
      } else {
        stale_steps <- 0L
      }
      if (stale_steps >= patience) {
        converged <- TRUE
        loss_history <- loss_history[seq_len(i)]
        break
      }
    }
  }

  if (is.null(best_spatial)) {
    best_spatial <- as.matrix(A$detach()$cpu())
  }

  list(
    spatial = best_spatial,
    reconstruction = best_spatial %*% spectra,
    loss = loss_history,
    niter_run = length(loss_history),
    converged = converged
  )
}

fit_spectral_unmix_prediction_with_restarts <- function(x,
                                                        spectra,
                                                        init = c("least_squares", "random"),
                                                        nstart = 1L,
                                                        seed = NULL,
                                                        lr,
                                                        niter,
                                                        tol,
                                                        patience,
                                                        device,
                                                        verbose = FALSE) {
  init <- match.arg(init)
  restart_seeds <- generate_restart_seeds(nstart, seed = seed)
  restart_losses <- numeric(nstart)
  best_fit <- NULL
  best_loss <- Inf

  for (i in seq_len(nstart)) {
    init_spatial <- initialize_fixed_spectra_weights(
      x = x,
      spectra = spectra,
      init = init,
      seed = restart_seeds[i]
    )
    set_optimization_seed(restart_seeds[i])
    fit <- fit_spectral_unmix_prediction_adam(
      x = x,
      spectra = spectra,
      init_spatial = init_spatial,
      lr = lr,
      niter = niter,
      tol = tol,
      patience = patience,
      device = device,
      verbose = verbose
    )

    restart_losses[i] <- utils::tail(fit$loss, 1L)
    if (restart_losses[i] < best_loss) {
      best_loss <- restart_losses[i]
      best_fit <- fit
      best_fit$restart_seed <- restart_seeds[i]
    }
  }

  best_fit$restart_losses <- restart_losses
  best_fit
}

extract_matrix_metadata <- function(x) {
  metadata <- attr(x, "spectral_unmix_metadata", exact = TRUE)
  dims <- attr(x, "spectral_unmix_dim", exact = TRUE)

  if (!is.null(dims)) {
    metadata <- merge_metadata(metadata, list(dim_xy = dims))
  }

  metadata
}

merge_metadata <- function(x, y) {
  if (is.null(x) || length(x) == 0L) {
    return(y)
  }
  if (is.null(y) || length(y) == 0L) {
    return(x)
  }

  utils::modifyList(x, y)
}
