#'
#' Detection limit of the Bigelow model
#' 
#' Calculation of the detection limit for the Bigelow model
#' 
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param temperature numerical value that describes the temperature at which the detection limit will be calculated 
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @return Numerical value that indicates the limit of detection for that temperature for the Bigelow model
#'
#' @examples  
#' pars <- list(temp_ref = 55,
#'         z = 5.18 ,
#'         D_R = 12.10 )
#' detection_bigelow( pars, temperature = 57, limit=7)
#' @export
#' 


detection_bigelow <- function(pars, temperature, limit) {
    
    limit * pars$D_R * 10^(-(temperature - pars$temp_ref)/pars$z) 
}

#' Detection limit of the Peleg model
#' 
#' Calculation of the detection limit for the Peleg model
#' 
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param temperature numerical value that describes the temperature at which the detection limit will be calculated
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @return Numerical value that indicates the limit of detection for that temperature for the Peleg model
#'
#' @examples  
#' pars <- list(temp_crit = 56.95,
#'         k_b = 0.58 ,
#'         n = 1 )
#' detection_peleg( pars, temperature = 57, limit=7)
#' @export
#' 
detection_peleg <- function(pars, temperature, limit) {
    
    ( limit/log(1 + exp(pars$k_b*(temperature - pars$temp_crit)) ) )^(1/pars$n)
    
}

#' Detection limit of the Mafart model
#' 
#' Calculation of the detection limit for the Mafart model
#' 
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param temperature numerical value that describes the temperature at which the detection limit will be calculated 
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#' @return Numerical value that indicates the limit of detection for that temperature for the Mafart model
#'
#' @examples  
#' pars <- list(temp_ref = 55,
#'         z = 5.18 ,
#'         p = 0.99 ,
#'         delta_ref = 11.96)
#' detection_mafart( pars, temperature = 57, limit=7)
#' @export
#' 
detection_mafart <- function(pars, temperature, limit) {
    
    limit^(1/pars$p) * pars$delta_ref * 10^(-(temperature - pars$temp_ref)/pars$z) 
}

#' Calculate detection limit
#' 
#' Calculation of the detection limit depending on the model.
#' 
#' @param model character string defining the inactivation model to use.
#' @param pars list defining the model parameters according to the rules defined in the bioinactivation package.
#' @param temperature numerical value that describes the temperature at which the detection limit will be calculated
#' @param limit numerical value describing the maximum number of log-reductions
#' that can be identified in the experiment limit = logDL - logN0, where DL
#' is the detection limit.
#'
#' @return Numerical value that indicates the limit of detection
#'
#' @export
#' 
get_detection <- function(model, pars, temperature, limit) {
    
    switch(model,
           Bigelow = detection_bigelow(pars, temperature, limit),
           Mafart = detection_mafart(pars, temperature, limit),
           Peleg = detection_peleg(pars, temperature, limit)
           )
}







