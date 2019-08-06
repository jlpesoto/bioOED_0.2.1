
#'
#' Calculation of Fisher Information Matrix
#' 
#' The sensitivities at the different \code{times} are calculated by 
#' linear interpolation of the results provided in \code{sensitivities}.
#' 
#' @param sensitivities data.frame of class \code{sensFun} as returned by
#' \code{\link{sensitivity_inactivation}}.
#' @param times Numeric vector of time points where observations will be taken.
#' 
#' @return Matrix with the estimation of the Fisher Information Matrix.
#' 
#' @importFrom dplyr select_
#' @importFrom dplyr %>%
#' @importFrom stats approx
#' 
#' @export
#' 
calculate_FIM <- function(sensitivities, times) {
    
    old_times <- sensitivities$x
    
    sensitivities <- select_(sensitivities, quote("-x"), quote("-var")) # %>%
        # as.matrix()
    
    interp_sensitivities <- matrix(NA, length(times), ncol(sensitivities))
    
    for (i in 1:ncol(sensitivities)) {
        interp_sensitivities[, i] <- approx(old_times, sensitivities[,i],
                                            times)$y
    }
    
    FIM <- t(interp_sensitivities) %*% interp_sensitivities
    FIM
}

#'
#' D Optimality Criterium
#' 
#' @param FIM Matrix with the values of the Fisher Information Matrix
#' 
criterium_D <- function(FIM) {
    
    -det(FIM)
    
}

#'
#' Modified-E Optimality Criterium
#' 
#' @param FIM Matrix with the values of the Fisher Information Matrix
#' @param eig_tol Tolerance for the eigen values. If any eigen value is lower than this
#' value, the FIM is singular and a high value (1e20) is returned. 1e-10 by default.
#' 
criterium_modE <- function(FIM, eig_tol = 1e-10) {

    eig_vals <- eigen(FIM, only.values=TRUE)
    
    min_eig <- min(eig_vals$values)
    
    if (abs(min_eig) < eig_tol) {
        
        return(1e20)
        
    } else {
        
        abs(max(eig_vals$values)/min_eig)
        
    }
}
















