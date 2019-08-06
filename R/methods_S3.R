
#'
#' Correlation Plot of Parameter Sensitivities
#' 
#' Makes a correlation plot of the sensitivities between model
#' parameters.
#' 
#' @param x Instance of parCorrelation
#' @param y Ignored
#' @param ... Ignored
#' 
#' @importFrom corrplot corrplot
#' 
#' @method plot parCorrelation
#' 
#' @export
#' 
plot.parCorrelation <- function(x, y = NULL, ...) {
    
    corrplot(x, method="color",  
             type="upper", order="hclust", 
             addCoef.col = "black", # Add coefficient of correlation
             tl.col="black", tl.srt=45, #Text label color and rotation
             # hide correlation coefficient on the principal diagonal
             diag=FALSE )
}

#'
#' Plot of OEDinactivation 
#' 
#' @importFrom ggplot2 geom_vline
#' @importFrom bioinactivation predict_inactivation
#' 
#' @param x An instance of OEDinactivation
#' @param y Ignored
#' @param ... Ignored
#' 
#' @importFrom graphics plot
#' 
#' @method plot OEDinactivation
#' 
#' @export
#' 
plot.OEDinactivation <- function(x, y = NULL, ...) {
    
    temp_profile <- x$temp_profile
    times <- seq(0, max(temp_profile$time), length = 100)
    
    prediction <- predict_inactivation(x$model, times,
                                       c(x$parms, x$parms_fix),
                                       temp_profile)
    
    plot(prediction) + geom_vline(xintercept = x$optim_times,
                                  linetype = 2, colour = "red", size = 1)
}






