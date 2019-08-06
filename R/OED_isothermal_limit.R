
#'
#' OED of isothermal microbial inactivation with detection limit
#'
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the nominal model parameters.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @param n_points numerical stating the number of data points.
#' @param min_time numerical stating the lower limit for the time points.
#' @param max_time numerical stating the upper limit for the time points.
#' @param min_temp numerical stating the lower limit for the temperature.
#' @param max_temp numerical stating the upper limit for the temperature.
#' @param criterion character string defining the criterion to use.
#' @param opts options for the MEIGO algorithm. By default, a maximum of 2000
#' function evaluations with local finish with the DHC algorithm
#' (see help from MEIGO).
#' @param x_0 initial point for the MEIGO algorithm. By default, it is NULL.
#'
#' @return A MEIGO object
#'
#' @import tidyverse
#'
#' @export
#'
#' @examples
#' pars <- list(z = 4.2, D_R = 3.9, temp_ref = 55)
#' opts <- list(maxeval=500,local_finish="DHC")
#' OED <- isothermal_OED_limit("Bigelow", pars, n_points = 5, criterion = "E-mod", limit = 6,
#'                       min_time = 0, max_time = 100, min_temp = 52.5, max_temp = 60,
#'                       opts = opts)
#' plot(OED)
#'
isothermal_OED_limit <- function(model, pars, limit,
                                 n_points, min_time, max_time, min_temp, max_temp, criterion = "D",
                                 opts = NULL, x_0=NULL) {

  if (min_time <= 0) {
    min_time <- 1e-6
    print("NOTE: min_time has been set to 1e-6 to avoid singularities in Weibullian models")
  }

  tgt_function = switch(criterion,
                        D = criterium_D_iso,
                        `E-mod` = criterium_Emod_iso,
                        # E = criterium_E_iso,
                        `A-mod` = criterium_Amod_iso,
                        # A = criterium_A_iso,
                        stop(paste("Unknown criterion:", criterion))
  )

  problem <- list(f = tgt_function,
                  x_0 = x_0,
                  x_L = c(rep(min_time, n_points),rep(min_temp, n_points)),
                  x_U = c(rep(max_time, n_points),rep(max_temp, n_points))
  )

  if (is.null(opts)) {

    opts <- list(maxeval=2000,local_finish="DHC")
  }

  result <- MEIGO(problem, opts, algorithm="ESS",
                  model = model, pars = pars, limit = limit)

  ## Map the results back

  half <- length(result$xbest)/2

  time_points <- result$xbest[1:half]
  temp_points <- result$xbest[(half+1):length(result$xbest)]

  my_design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select("times", "temperature")

  ## Return

  out <- list(
    optim = result,
    model = model,
    pars = pars,
    criterion = "D",
    optim_algorithm = "MEIGO",
    optim_design = my_design,
    limit = limit
  )

  class(out) <- c("OEDisothermal", class(out))

  out

}
