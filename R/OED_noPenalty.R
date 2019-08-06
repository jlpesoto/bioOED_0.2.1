
#'
#' Optimum Experimental Design of Microbial Inactivation
#' 
#' Performs an optimum experimental design for the settings selected. The
#' OED is based on the FIM, estimated using the local sensitivity functions
#' provided by \code{\link{sensitivity_inactivation}}.
#' 
#' @importFrom MEIGOR MEIGO
#' @importFrom stats optim
#' @importFrom stats runif
#' @importFrom bioinactivation predict_inactivation
#' 
#' @param inactivation_model Character string defining the inacivation model.
#' @param parms Named numeric vector defining the model parameters. They must
#'        be named according to the needs of \code{\link{predict_inactivation}}.
#' @param temp_profile Data frame defining the temperature profile. It must 
#'        contain a column named \code{time} and a column named
#'        \code{temperature}.
#' @param parms_fix Named numeric vector defining the model parameters to be 
#'        omitted during the calculation of the local sensitivities.
#' @param n_points Number of measurements which will be taken during the
#'        experiment.
#' @param criteria Character defining the criteria for the OED. Either 
#'        \code{D} (default) or \code{E-mod}.
#' @param n_times Integer defining th enumber of discrete time points used for
#'        the interpolation of the local sensitivities.
#' @param sensvar Character defining the variable to use for the OED. Either
#'        \code{logN} (default) or \code{N}.
#' @param optim_algorithm Character defining the type of algorithm to use for
#'        the optimization. Either \code{global} (default) or \code{local}.
#' @param opts_global List defining the options for the global optimization
#'        algorithm (see \code{\link{MEIGO}}). By default, global solver with
#'        a maximum of 50000 function evaluations and printout on every step.
#' 
#' @export
#' 
#' @return A list of class \code{OEDinactivation} with the following items:
#'      \itemize{
#'          \item optim: Objetc returned by the optimization function.
#'          \item model: Inactivation model used for the calculations.
#'          \item parms: Nominal model parameters.
#'          \item parms_fix: Model parameters not considered for the
#'                sensitivity calculation.
#'          \item criteria: Criteria used for the OED.
#'          \item sensvar: Variable used for the OED.
#'          \item optim_algorithm: Type of optimization algorithm.
#'          \item optim_times: Optimum measurement times calculated.
#'          \item penalty: Logical indicating whether penalty function was
#'                used.
#'          \item temp_profile: Temperature profile of the experiment.
#'          }
#' 
#' @examples
#' ## Definition of input variables
#' 
#' parms_fix <- c(temp_ref = 57.5)
#' parms <- c(delta_ref = 3.9,
#'            z = 4.2,
#'            p = 1,
#'            N0 = 1e6
#'            )
#' 
#' temp_profile <- data.frame(time = c(0, 60), temperature = c(30, 60))
#' 
#' n_points <- 5
#' 
#' ## OED with local optimization
#' 
#' set.seed(191210)
#' 
#' local_OED <- inactivation_OED("Mafart", parms, temp_profile, parms_fix,
#'                       n_points, criteria = "E-mod", sensvar = "logN",
#'                       optim_algorithm = "local")
#' 
#' print(local_OED$optim_times)
#' plot(local_OED)
#'
inactivation_OED <- function(inactivation_model, parms, temp_profile, parms_fix,
                     n_points, criteria = "D",
                     n_times = 100, sensvar = "logN",
                     optim_algorithm = "global",
                     opts_global = NULL) {
    
    ## Calculate sensitivities
    
    sensitivities <- sensitivity_inactivation(inactivation_model, parms,
                                              temp_profile, parms_fix,
                                              n_times = n_times, sensvar = sensvar)
    
    ## Prepare for optimization
    
    if (grepl(criteria, "D")) {
        
        tgt_function <- objective_D
        
    } else if (grepl(criteria, "E-mod")) {
        
        tgt_function <- objective_Emod
        
    } else {
        
        stop(paste("Unknown criteria:", criteria))
        
    }
    
    max_time <- max(temp_profile$time)
    problem <- list(f=tgt_function,
                    x_L=rep(1e-3, n_points),  # 1e-3 to avoid singularity
                    x_U=rep(max_time, n_points)
                    )
    
    if (is.null(opts_global)) {
        opts_global <- list(maxeval=50000,  local_solver=0,
                            local_finish="DHC", local_iterprint=1)
    }
    
    ## Make optimization
    
    if (grepl(optim_algorithm, "local")) {
        
        times <- runif(n_points, 0, max_time)
        
        results <- optim(times, tgt_function,
                         lower = problem$x_L, upper = problem$x_U,
                         sensitivities = sensitivities,
                         method = "L-BFGS-B")
        
        best_points <- results$par
        
    } else if (grepl(optim_algorithm, "global")) {
        results <- MEIGO(problem, opts_global, algorithm="ESS",
                         sensitivities = sensitivities)
        
        best_points <- results$xbest
        
    } else {
        stop(paste("Unknown optimization algorithm:", optim_algorithm))
    }
    
    ## Return results
    
    out <- list(optim = results)
    
    class(out) <- c("OEDinactivation", class(out))
    
    out$model <- inactivation_model
    out$parms <- parms
    out$parms_fix <- parms_fix
    out$criteria <- criteria
    out$sensvar <- sensvar
    out$optim_algorithm <- optim_algorithm
    out$optim_times <- sort(best_points)
    out$penalty <- FALSE
    out$time_min <- NA
    out$temp_profile <- temp_profile
    
    out
}

#'
#' Objective Function for the D Criterium
#' 
#' @param times A numeric vector of points where the FIM will be calculated.
#' @param sensitivities An object returned by sensitivity_inactivation.
#' 
objective_D <- function(times, sensitivities) {
    
    FIM <- calculate_FIM(sensitivities, times)
    
    out <- criterium_D(FIM)
    
    out
    
}

#'
#' Objective Function for the modified-E Criterium
#' 
#' @param times A numeric vector of points where the FIM will be calculated.
#' @param sensitivities An object returned by sensitivity_inactivation.
#' 
objective_Emod <- function(times, sensitivities) {
    
    FIM <- calculate_FIM(sensitivities, times)
    
    out <- criterium_modE(FIM)
    
    out
    
}








