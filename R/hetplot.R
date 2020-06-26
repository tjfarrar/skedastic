#' Graphical Methods for Detecting Heteroskedasticity in a Linear Regression Model
#'
#' This function creates various two-dimensional scatter plots that can aid in
#' detecting heteroskedasticity in a linear regression model.
#'
#' The variable plotted on the horizontal axis could be the original
#' data indices, one of the explanatory variables, the OLS predicted (fitted)
#' values, or any other numeric vector specified by the user. The variable
#' plotted on the vertical axis is some function of the OLS residuals or
#' transformed version thereof such as the BLUS residuals
#' \insertCite{Theil68;textual}{skedastic} or standardised or studentised
#' residuals as discussed in \insertCite{Cook83;textual}{skedastic}. A separate
#' plot is created for each (\code{horzvar}, \code{vertvar}, \code{vertfun})
#' combination.
#'
#' @param horzvar A character vector describing the variable(s) to plot on
#'    horizontal axes (\code{"index"} for the data index \eqn{i},
#'    \code{"fitted.values"} for the OLS predicted values \eqn{\hat{y}_i},
#'    \code{"fitted.values2"} for transformed OLS predicted values
#'    \eqn{m_{ii}\hat{y}_i}, and/or \code{names} of explanatory variable columns).
#'    \code{"explanatory"} passes all explanatory variable columns. \code{"log"}
#'    concatenated with \code{names} of explanatory variable columns passes logs
#'    of those explanatory variables. \code{"log_explanatory"} passes logs of
#'    all explanatory variables. If more than one variable is specified,
#'    a separate plot is created for each.
#' @param vertvar A character vector describing the variable to plot on the
#'    vertical axis ("res" for OLS residuals [the default], "res_blus" for
#'    \code{\link[=blus]{BLUS}} residuals, "res_stand" for standardised OLS
#'    residuals: \eqn{e_i/\hat{\sigma}}, "res_constvar" for OLS
#'    residuals transformed to have constant variance: \eqn{e_i/\sqrt{m_{ii}}},
#'    "res_stud" for studentised OLS residuals: \eqn{e_i/(s\sqrt{m_{ii}})}.
#'    If more than one value is specified, a separate plot is created for each.
#' @param vertfun A character vector giving the names of functions to apply to
#'    the \code{vertvar} variable. Numerals such as "2" are taken to be powers
#'    to which \code{vertvar} should be set. If multiple values are specified,
#'    they are all applied to each element of \code{vertvar}.
#' @param filetype A character giving the type of image file to which the
#'    plot(s) should be written. Values can be \code{"png"},
#'    \code{"bmp"}, \code{"jpeg"}, or \code{"tiff"}. Image files are written
#'    to a subdirectory called "hetplot" within the R session's temporary
#'    directory, which can be located using \code{tempdir()}. The files should
#'    be moved or copied to another location if they are needed after the R
#'    session is ended. Default filenames contain timestamps for uniqueness.
#'    If \code{NULL} (the default), no image files are written, and in this
#'    case, if there are multiple plots, they are plotted on a single device
#'    using the \code{"mfrow"} graphical parameter. If many plots are requested
#'    at once, it is advisable to write them to image files.
#' @param values A logical. Should the sequences corresponding to the
#'    horizontal and vertical variable(s) be returned in a \code{list} object?
#' @param ... Arguments to be passed to methods, such as graphical parameters
#'    (see \code{\link{par}}), parameters for \code{\link[graphics]{plot}},
#'    for \code{\link[grDevices:png]{graphics devices}},
#'    and/or the \code{omit} argument for function \code{\link{blus}},
#'    if BLUS residuals are being plotted. If it is desired to pass the
#'    \code{type} argument to a graphics device, use \code{gtype = }, since
#'    a \code{type} argument will be passed to \code{\link[graphics]{plot}}.
#' @inheritParams breusch_pagan
#'
#' @return A list containing two \code{\link[base:data.frame]{data frames}}, one
#'    for vectors plotted on horizontal axes and one for vectors plotted
#'    on vertical axes.
#' @references{\insertAllCited{}}
#' @importFrom Rdpack reprompt
#' @export
#' @seealso \code{\link[stats]{plot.lm}}
#'
#' @examples
#' mtcars_lm <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' # Generates 2 x 2 matrix of plots in console
#' hetplot(mtcars_lm, horzvar = c("index", "fitted.values"),
#' vertvar = c("res_blus"), vertfun = c("2", "abs"), filetype = NULL)
#'
#' # Generates 84 png files in tempdir() folder
#' ## Not run
#' # hetplot(mainlm = mtcars_lm, horzvar = c("explanatory", "log_explanatory",
#' #  "fitted.values2"), vertvar = c("res", "res_stand", "res_stud",
#' #  "res_constvar"), vertfun = c("identity", "abs", "2"), filetype = "png")
#'

hetplot <- function (mainlm, horzvar = 1:n, vertvar = "res", vertfun = "identity",
                     filetype = NULL, values = FALSE, ...) {

  args <- list(...)
  grdevargnames <- c("filename", "width", "height", "units", "pointsize", "bg",
                   "res", "family", "restoreConsole", "compression", "antialias",
                   "quality", "gtype")
  plotargnames <- c("type", "main", "sub", "xlab", "ylab", "asp")
  if ("omit" %in% names(args)) {
    omitarg <- args[["omit"]]
    args[["omit"]] <- NULL
  }
  parargs <- args[setdiff(names(args), c(grdevargnames, plotargnames))]
  plotargs <- args[intersect(names(args), plotargnames)]
  grdevargs <- args[intersect(names(args), grdevargnames)]
  names(grdevargs)[which(names(grdevargs) == "gtype")] <- "type"

  filename_passed <- ("filename" %in% names(grdevargs))
  mar_passed <- ("mar" %in% names(parargs))
  las_passed <- ("las" %in% names(parargs))
  cex_passed <- ("cex" %in% names(parargs))
  main_passed <- ("main" %in% names(plotargs))
  xlab_passed <- ("xlab" %in% names(plotargs))
  ylab_passed <- ("ylab" %in% names(plotargs))

  processmainlm(m = mainlm, needyhat = TRUE)

  n <- nrow(X)
  M <- fastM(X, n)
  sigma_hat_sq <- sum(e ^ 2) / n
  s_sq <- sum(e ^ 2) / (n - p)

  hasintercept <- columnof1s(X)
  if (hasintercept[[1]]) {
    if (hasintercept[[2]] != 1) stop("Column of 1's must be first column of design matrix")
    colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
  } else {
    colnames(X) <- paste0("X", 1:p)
  }

  if ("explanatory" %in% horzvar) {
    horzvar <- c(horzvar[-which(horzvar == "explanatory")], colnames(X))
    if ("(Intercept)" %in% horzvar) horzvar <- horzvar[-which(horzvar == "(Intercept)")]
  }
  if ("log_explanatory" %in% horzvar) {
    horzvar <- c(horzvar[-which(horzvar == "log_explanatory")], paste0("log", colnames(X)))
    if ("log(Intercept)" %in% horzvar) horzvar <- horzvar[-which(horzvar == "log(Intercept)")]
  }

  x_hor <- as.data.frame(lapply(horzvar, function(x) if (x %in% colnames(X)) {X[, x]}
                                    else if (regexpr("log", substr(x, 1, 3)) != -1) {
                                      log(X[, substr(x, 4, nchar(x))])
                                    }
                                    else if (x == "fitted.values") {yhat}
                                    else if (x == "fitted.values2") {diag(M) * yhat}
                                    else if (x == "index") {1:n}))
  names(x_hor) <- horzvar

  vertfun_function <- lapply(vertfun, function(x) {
      if (suppressWarnings(is.na(as.integer(x)))) {get(x)}
      else {function(y) `^`(y, as.integer(x)) }
    })

  y_ver_nofunc <- lapply(vertvar, function(x) if (x == "res") {e}
                                else if (x == "res_blus") {
                                  if (exists("omitarg")) {
                                    blus(mainlm, omit = omitarg, keepNA = TRUE)
                                  } else { blus(mainlm, keepNA = TRUE)}
                                }
              else if (x == "res_stand") {e / sqrt(sigma_hat_sq)}
              else if (x == "res_constvar") {e / sqrt(diag(M))}
              else if (x == "res_stud") {e / sqrt(diag(M) * s_sq)})

  y_ver <- as.data.frame(unlist(lapply(vertfun_function,
                                function(x) lapply(y_ver_nofunc, x)),
                                recursive = FALSE))
  names(y_ver) <- unlist(lapply(vertfun, function(x) paste(vertvar, x,
                                                           sep = "_")))

  yline_mtext <- unlist(lapply(names(y_ver), function(l) {
    ylinebase <- ifelse(is.null(filetype), 3, 4.25)
    if (regexpr("stud", l) != -1) {
      ylinebase <- ylinebase - 0.5
    } else if (regexpr("const", l) != -1) {
      ylinebase <- ylinebase - 0.25
    }
    if (regexpr("identity", l) == -1) {
      ylinebase <- ylinebase - 0.5
    }
    ylinebase
  }))

  opar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(opar))
  numplots <- length(horzvar) * length(vertvar) * length(vertfun)
  if (is.null(filetype)) {
    mfrow_dim <- plotdim(numplots)
    if (max(mfrow_dim) >= 5) warning("Large number of plots in one dimension")
    graphics::par(mfrow = mfrow_dim)
    if (!las_passed) parargs$las <- 1
    if (!mar_passed) parargs$mar <- c(4, 5.5, 0.5, 1.5)
    do.call(graphics::par, parargs)
    mapply(function(y, ynames, yline) mapply(function(x, xnames, y, ynames, yline) {
      if (!main_passed) plotargs$main <- ""
      if (!xlab_passed) plotargs$xlab <- ""
      if (!ylab_passed) plotargs$ylab <- ""
      plotargs$x <- x
      plotargs$y <- y
      do.call(graphics::plot, plotargs)
      if (!xlab_passed) graphics::mtext(parselabels(xnames, theaxis = "x"),
                              side = 1, line = 2.5, las = 1,
                              cex = ifelse(cex_passed, parargs$cex, 0.8))
      if (!ylab_passed) graphics::mtext(parselabels(ynames, theaxis = "y"),
                              side = 2, line = yline, las = 1,
                              cex = ifelse(cex_passed, parargs$cex, 0.8))
    }, x_hor, names(x_hor), MoreArgs = list(y, ynames, yline), SIMPLIFY = FALSE),
    y_ver, names(y_ver), yline_mtext, SIMPLIFY = FALSE)
  } else {
    if (!dir.exists(paste0(tempdir(),"/hetplot"))) {
      dir.create(paste0(tempdir(),"/hetplot"))
    }
    mapply(function(y, ynames, yline) mapply(function(x, xnames, y, ynames, yline) {
      if (!filename_passed) {
        grdevargs$filename <- paste0(tempdir(),"/hetplot/",xnames, "_",
                          ynames, "_", gsub("[[:space:]]|[[:punct:]]",
                                  "_", Sys.time()), ".", filetype)
      }
      do.call(get(filetype), grdevargs)
      if (!las_passed) parargs$las <- 1
      if (!mar_passed) parargs$mar <- c(4, 7.25, 0.1, 0.1)
      do.call(graphics::par, parargs)
      if (!main_passed) plotargs$main <- ""
      if (!xlab_passed) plotargs$xlab <- ""
      if (!ylab_passed) plotargs$ylab <- ""
      plotargs$x <- x
      plotargs$y <- y
      do.call(graphics::plot, plotargs)
      if (!xlab_passed) graphics::mtext(parselabels(xnames, theaxis = "x"),
          side = 1, line = 2.5, las = 1, cex = ifelse(cex_passed, parargs$cex, 1.2))
      if (!ylab_passed) graphics::mtext(parselabels(ynames, theaxis = "y"),
          side = 2, line = yline, las = 1, cex = ifelse(cex_passed, parargs$cex, 1.2))
      if (length(grDevices::dev.list()) > 59) {
        grDevices::graphics.off()
      }
    }, x_hor, names(x_hor), MoreArgs = list(y, ynames, yline), SIMPLIFY = FALSE),
           y_ver, names(y_ver), yline_mtext, SIMPLIFY = FALSE)
    grDevices::graphics.off()
  }
  if (values) list("horizontal" = x_hor, "vertical" = y_ver)
}
