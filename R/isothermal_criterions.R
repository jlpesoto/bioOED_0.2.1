#'
#' Objective function for D-optimal OED with detection limit
#' 
#' Points outside of the allowable area are moved back in time to the 
#' detection limit
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @import tidyselect
#' @importFrom dplyr mutate 
#' @return Numeric value of the objective function for criterium D, which is a determinant of the FIM.
#'
criterium_D_iso <- function(x, model, pars, limit){
  
  half <- length(x)/2
  
  time_points <- x[1:half]
  temp_points <- x[(half+1):length(x)]
  
  ## read and map
  
  design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select("times", "temperature")
  
  -det(calculate_isothermal_FIM(model, design, pars))
}

#'
#' Objective function for E modified-optimal OED with detection limit
#' 
#' Points outside of the allowable area are moved back in time to the 
#' detection limit
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @importFrom dplyr mutate
#' @return Numeric value of the objective function for criterium E modified, which is a determinant of the FIM.
#' 
criterium_Emod_iso <- function(x, model, pars, limit) {
  tol_eigen <- 1e-100
  half <- length(x)/2
  
  time_points <- x[1:half]
  temp_points <- x[(half+1):length(x)]
  
  ## read and map
  
  design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select("times", "temperature")
  eigenvalues <- eigen(calculate_isothermal_FIM(model, design, pars))
  if(abs(min(eigenvalues$values))-tol_eigen<0){
    return(1e6)
  }
  else if(abs(min(eigenvalues$values))-tol_eigen>0) {
    return(abs(max(eigenvalues$values)/min(eigenvalues$values)))
  }
}

#'
#' Objective function for E-optimal OED with detection limit
#' 
#' Points outside of the allowable area are moved back in time to the 
#' detection limit
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @importFrom dplyr mutate
#' @return Numeric value of the objective function for criterium E, which is a determinant of the FIM.
#'
criterium_E_iso <- function(x, model, pars, limit) {
  tol_det <- 1e-5
  half <- length(x)/2
  
  time_points <- x[1:half]
  temp_points <- x[(half+1):length(x)]
  
  ## read and map
  
  design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select("times", "temperature")
  if(abs(det(calculate_isothermal_FIM(model, design, pars)))-tol_det<0) {
    return(1e100)
  }
  else {
    eigenvalues <- eigen(solve(calculate_isothermal_FIM(model, design, pars)))
    return(max(eigenvalues$values))
    
  }
}

#'
#' Objective function for A modified-optimal OED with detection limit
#' 
#' Points outside of the allowable area are moved back in time to the 
#' detection limit
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @importFrom dplyr mutate
#' @return Numeric value of the objective function for criterium A modified, which is a determinant of the FIM.
#' 
criterium_Amod_iso <- function(x, model, pars, limit) {
  half <- length(x)/2
  time_points <- x[1:half]
  temp_points <- x[(half+1):length(x)]
  
  ## read and map
  
  design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select("times", "temperature")
  return(-sum(diag(calculate_isothermal_FIM(model, design, pars))))
  
}

#'
#' Objective function for A-optimal OED with detection limit
#' 
#' Points outside of the allowable area are moved back in time to the 
#' detection limit
#' 
#' @param x a numeric vector of length \code{n} defining the design matrix.
#' The first n/2 elements are the time points and the last n/2 are the 
#' temperatures of these points.
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @importFrom dplyr mutate
#' @return Numeric value of the objective function for criterium A, which is a determinant of the FIM.
#'
criterium_A_iso <- function(x, model, pars, limit) {
  tol_det <- 1e-5
  half <- length(x)/2
  
  time_points <- x[1:half]
  temp_points <- x[(half+1):length(x)]
  
  ## read and map
  
  design <- data.frame(times = time_points, temperature = temp_points) %>%
    mutate(det_limit = get_detection(model, pars, .data$temperature, limit)) %>%
    mutate(times = ifelse(.data$times > .data$det_limit, .data$det_limit, .data$times)) %>%
    select("times", "temperature")
  if(abs(det(calculate_isothermal_FIM(model, design, pars)))-tol_det<0) {
    return(1e100)
  }
  else {
    return(sum(diag(solve(calculate_isothermal_FIM(model, design, pars)))))
    
  }
}
