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
#' @param lr Learning rate passed to \code{torch::optim_adam()}.
#' @param niter Number of optimization iterations.
#' @param center Logical; if \code{TRUE}, center wavelength channels before
#'   fitting.
#' @param scale Logical; if \code{TRUE}, scale wavelength channels before
#'   fitting.
#' @param cuda Logical; if \code{TRUE}, request CUDA execution.
#' @param verbose Logical; if \code{TRUE}, emit progress every 100 iterations.
#'
#' @return A list with class \code{"spectral_unmix"} containing:
#' \itemize{
#'   \item \code{spatial}: abundance matrix with one column per component.
#'   \item \code{spectra}: component spectra with one row per component.
#'   \item \code{reconstruction}: reconstructed matrix \eqn{AS}.
#'   \item \code{loss}: objective history over the optimization.
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
                           lr = 0.02,
                           niter = 2000,
                           center = FALSE,
                           scale = FALSE,
                           cuda = FALSE,
                           verbose = FALSE) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop("The 'torch' package must be installed to use spectral_unmix().", call. = FALSE)
  }

  x <- validate_input_matrix(x)

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
  if (!is.numeric(lr) || length(lr) != 1L || is.na(lr) || lr <= 0) {
    stop("'lr' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(niter) || length(niter) != 1L || is.na(niter) || niter < 1) {
    stop("'niter' must be a positive integer.", call. = FALSE)
  }
  niter <- as.integer(niter)

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
  X <- torch::torch_tensor(as.matrix(x), device = device)

  n_spectra <- nrow(x)
  n_wave <- ncol(x)

  A <- torch::torch_rand(n_spectra, k, device = device, requires_grad = TRUE)
  S <- torch::torch_rand(k, n_wave, device = device, requires_grad = TRUE)
  opt <- torch::optim_adam(list(A, S), lr = lr)

  loss_history <- numeric(niter)

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
    if (verbose && (i %% 100L == 0L || i == 1L || i == niter)) {
      message(sprintf("Iteration %d/%d, loss = %.6f", i, niter, loss_history[i]))
    }
  }

  fit <- list(
    spatial = as.matrix(A$detach()$cpu()),
    spectra = as.matrix(S$detach()$cpu()),
    reconstruction = as.matrix(A$detach()$cpu()) %*% as.matrix(S$detach()$cpu()),
    loss = loss_history,
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
#' @param cube Numeric 3D array.
#'
#' @return A numeric matrix with \code{nx * ny} rows and \code{n_wave} columns.
#' @examples
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' x <- cube_to_matrix(cube$cube)
#' dim(x)
#' @export
cube_to_matrix <- function(cube) {
  if (!is.array(cube) || length(dim(cube)) != 3L) {
    stop("'cube' must be a numeric 3D array.", call. = FALSE)
  }
  if (!is.numeric(cube) || any(!is.finite(cube))) {
    stop("'cube' must contain finite numeric values.", call. = FALSE)
  }

  dims <- dim(cube)
  matrix(cube, nrow = dims[1L] * dims[2L], ncol = dims[3L], byrow = FALSE)
}

#' Reshape a Matrix Back Into a Cube
#'
#' Reconstructs a \code{nx x ny x n_wave} cube from a matrix with one row per
#' spaxel.
#'
#' @param x Numeric matrix with one row per spatial pixel.
#' @param nx Number of x-axis pixels.
#' @param ny Number of y-axis pixels.
#'
#' @return A numeric 3D array.
#' @examples
#' cube <- simulate_ifu_cube(nx = 6, ny = 5, n_wave = 20)
#' rebuilt <- matrix_to_cube(cube$matrix, cube$nx, cube$ny)
#' dim(rebuilt)
#' @export
matrix_to_cube <- function(x, nx, ny) {
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

  array(x, dim = c(nx, ny, ncol(x)))
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
  if (!is.numeric(component) || length(component) != 1L || component < 1L ||
      component > ncol(fit$spatial)) {
    stop("'component' is out of bounds.", call. = FALSE)
  }

  component <- as.integer(component)
  matrix(fit$spatial[, component], nrow = as.integer(nx), ncol = as.integer(ny))
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
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)

    graphics::par(mfrow = c(2, k), mar = c(3, 3, 2, 1))
    for (i in seq_len(k)) {
      graphics::matplot(
        demo$wavelength,
        fit$spectra[i, ],
        type = "l",
        lty = 1,
        xlab = "Wavelength",
        ylab = "Flux",
        main = sprintf("Spectrum %d", i)
      )
    }
    for (i in seq_len(k)) {
      graphics::image(
        component_map(fit, nx = demo$nx, ny = demo$ny, component = i),
        main = sprintf("Map %d", i),
        xlab = "x",
        ylab = "y",
        col = grDevices::hcl.colors(32, "YlOrRd", rev = TRUE)
      )
    }
  }

  list(data = demo, fit = fit)
}

#' @export
print.spectral_unmix <- function(x, ...) {
  cat("<spectral_unmix>\n", sep = "")
  cat(sprintf("  spaxels: %d\n", nrow(x$spatial)))
  cat(sprintf("  wavelength channels: %d\n", ncol(x$spectra)))
  cat(sprintf("  components: %d\n", nrow(x$spectra)))
  cat(sprintf("  final loss: %.6f\n", tail(x$loss, 1L)))
  invisible(x)
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
