
## Local sensitivities

#'
#' Local sensitivities of the Bigelow model
#' 
#' @importFrom rlang .data
#' @importFrom tidyr %>%
#' @importFrom dplyr mutate
#' 
#' @param exp_design data.frame with two columns named \code{times} and 
#' \code{temperature} describing the experiment design.
#' @param pars list defining the model parameters according to the rules defined
#' in the \code{bioinactivation} package.
#' 
#' @return A data frame with the same number of rows as \code{exp_design} with 
#' additional column for local sensitivities. These are named D_R and z for local
#' sensitivities and D_R_scaled and z_scaled for scaled locan sensitivities.
#' 
#'
sensitivities_Bigelow <- function(exp_design, pars) {
    
    exp_design %>%
        mutate(D_R = .data$times * 10^((.data$temperature - pars$temp_ref)/pars$z)/pars$D_R^2,
               z = .data$times*log(10)*(pars$temp_ref - .data$temperature)*10^(-(pars$temp_ref - .data$temperature)/pars$z)/pars$D_R/pars$z^2
        ) %>%
        mutate(D_R_scaled = .data$D_R/pars$D_R,
               z_scaled = .data$z/pars$z
        )
}


#'
#' Local sensitivities of the Peleg model
#' 
#' @inheritParams sensitivities_Bigelow
#' 
#' @importFrom rlang .data
#' @importFrom dplyr mutate 
#' @importFrom tidyr %>%
#'
sensitivities_Peleg <- function(exp_design, pars) {
    
    denom <- exp(exp_design$temperature * pars$k_b) + exp(pars$k_b * pars$temp_crit)
    
    exp_design %>%
        mutate(
            k_b = (pars$temp_crit - .data$temperature)*.data$times^pars$n*exp(.data$temperature*pars$k_b)/denom,
            temp_crit = pars$k_b*exp(pars$k_b*.data$temperature)*.data$times^pars$n/denom,
            n = -.data$times^pars$n * log(.data$times)*log(1 + exp(pars$k_b*(.data$temperature - pars$temp_crit)))
        ) %>%
        mutate(
            k_b_scaled = .data$k_b/pars$k_b,
            temp_crit_scaled = .data$temp_crit/pars$temp_crit,
            n_scaled = .data$n/pars$n
        )
}

#'
#' Local sensitivities of the Mafart model
#' 
#' @inheritParams sensitivities_Bigelow
#'
sensitivities_Mafart <- function(exp_design, pars) {
    
    z_def <- 10^(-(pars$temp_ref - exp_design$temperature)/pars$z)
    my_exp <-  (exp_design$times * z_def/pars$delta_ref)^(pars$p - 1)
    
    exp_design %>%
        mutate(
            delta_ref = pars$p*.data$times*z_def*my_exp/pars$delta_ref^2,
            z = -pars$p*.data$times*log(10)*(pars$temp_ref - .data$temperature)*z_def*my_exp/pars$delta_ref/pars$z^2,
            p = -(.data$times * z_def/pars$delta_ref)^pars$p*log(.data$times * z_def/pars$delta_ref)
            
        ) %>%
        mutate(
            delta_ref_scaled = .data$delta_ref/pars$delta_ref,
            z_scaled = .data$z/pars$z,
            p_scaled = .data$p/pars$p
        )
}

#' Local sensitivites of isothermal microbial inactivation
#' 
#' @inheritParams sensitivities_Bigelow
#' @param model character defining the inactivation model according to the rules
#' in the \code{bioinactivation} package.
#' 
#' @return A list of class \code{"IsoSensitivities"} with 3 entries:
#' \describe{
#'   \item{model}{Inactivation model.}
#'   \item{pars}{Model parameters used for the calculations.}
#'   \item{sensitivities}{data.frame adding columns to exp_design with the 
#'   calculated sensitivities. Local sensitivities are named as the parameters,
#'   scaled sensitivities as parameter_name+_scaled.}
#' }
#' 
#' @export
#' 
#' @examples 
#' 
#' library("tidyverse")
#' time_profile <- seq(0, 50, length = 20)
#' Temp_profile <- seq(52.5,60, length = 3)
#' 
#' exp_design <- expand.grid(time_profile,Temp_profile) %>%
#'   rename(times = Var1, temperature = Var2)
#' 
#' pars <- list(z = 4.2, D_R = 3.9, temp_ref = 55)
#' 
#' my_sensitivities <- isothermal_sensitivities("Bigelow", exp_design, pars)
#' plot(my_sensitivities)
#' plot(my_sensitivities, limit = 6)
#'
isothermal_sensitivities <- function(model, exp_design, pars) {
    
    sens <- switch(model,
                   Bigelow = sensitivities_Bigelow(exp_design, pars),
                   Mafart = sensitivities_Mafart(exp_design, pars),
                   Peleg = sensitivities_Peleg(exp_design, pars)
    )
    
    out <- list(model = model,
                pars = pars,
                sensitivities = sens)
    
    class(out) <- c("IsoSensitivities", class(out))
    
    out
    
}

## Correlation functions

#' Parameter correlation for isothermal inactivation experiments
#' 
#' @inheritParams isothermal_sensitivities
#' 
#' @importFrom dplyr select
#' 
#' @export
#'
#'@examples
#'
#' library(tidyverse)
#' time_profile <- seq(0, 50, length = 20)
#' Temp_profile <- seq(52.5,60, length = 3)
#' 
#' exp_design <- expand.grid(time_profile,Temp_profile) %>%
#'   rename(times = Var1, temperature = Var2)
#' 
#' pars <- list(temp_crit = 55,
#'              n = 1.5,
#'              k_b = 0.1)
#' 
#' get_isothermal_correlation("Peleg", exp_design, pars )
#' 

get_isothermal_correlation <- function(model, exp_design, pars) {
    
    out <- isothermal_sensitivities(model, exp_design, pars)$sensitivities %>%
        select(ends_with("_scaled")) %>%
        na.omit() %>%
        cor()
    
    class(out) <- c("IsoSensitivitiesCor", class(out))
    
    out
}

## Fisher Information Matrix

#'
#' Fisher Information Matrix for isothermal experiments
#' 
#' @inheritParams isothermal_sensitivities
#' 
#' @importFrom dplyr select
#' 
#' @export
#' 
#' @examples 
#' 
#' library("dplyr")
#' time_profile <- seq(0, 50, length = 20)
#' Temp_profile <- seq(52.5,60, length = 3)
#' 
#' exp_design <- expand.grid(time_profile,Temp_profile) %>%
#'   rename(times = Var1, temperature = Var2)
#' 
#' pars <- list(temp_crit = 55,
#'              n = 1.5,
#'              k_b = 0.1)
#' 
#' calculate_isothermal_FIM("Peleg", exp_design, pars )

calculate_isothermal_FIM <- function(model, exp_design, pars) {
    
    sens <- isothermal_sensitivities(model, exp_design, pars)$sensitivities %>%
        select(-ends_with("_scaled"), -"times", -"temperature") %>%
        as.matrix()
    
    t(sens) %*% sens
    
    
}
