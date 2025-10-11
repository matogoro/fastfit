#' Create Forest Plot with Automatic Model Detection
#'
#' Automatically detects the model type and creates an appropriate forest plot
#' using the corresponding specialized function (glmforest, coxforest, or lmforest).
#'
#' @param model A model object (glm, lm, coxph, or clogit)
#' @param data Data frame or data.table containing the original data
#' @param title Character string for plot title. If NULL, generates appropriate default
#' @param ... Additional arguments passed to the specific forest plot function
#'
#' @return A ggplot object containing the forest plot
#'
#' @examples
#' \dontrun{
#' # Automatically detects and plots any supported model type
#' model1 <- glm(outcome ~ x1 + x2, family = binomial, data = mydata)
#' forest(model1)
#' 
#' model2 <- coxph(Surv(time, status) ~ x1 + x2, data = mydata)
#' forest(model2)
#' 
#' model3 <- lm(y ~ x1 + x2, data = mydata)
#' forest(model3)
#' }
#'
#' @export
autoforest <- function(model, data = NULL, title = NULL, ...) {
    
    ## Detect model class
    model_class <- class(model)[1]
    
    ## Generate appropriate title if not provided
    if (is.null(title)) {
        title <- switch(model_class,
                        "coxph" = "Cox Proportional Hazards Model",
                        "clogit" = "Conditional Logistic Regression",
                        "glm" = {
                            if (model$family$family == "binomial") {
                                "Logistic Regression Model"
                            } else if (model$family$family == "poisson") {
                                "Poisson Regression Model" 
                            } else {
                                "Generalized Linear Model"
                            }
                        },
                        "lm" = "Linear Regression Model",
                        "Model Results"  # Generic fallback
                        )
    }
    
    ## Route to appropriate function
    if (model_class %in% c("coxph", "clogit")) {
        coxforest(model = model, data = data, title = title, ...)
        
    } else if (model_class == "glm") {
        glmforest(model = model, data = data, title = title, ...)
        
    } else if (model_class == "lm") {
        ## Check it's not actually a GLM
        if (inherits(model, "glm")) {
            glmforest(model = model, data = data, title = title, ...)
        } else {
            lmforest(model = model, data = data, title = title, ...)
        }
        
    } else {
        ## Check if it might be a supported class with different name
        if (inherits(model, "coxph")) {
            coxforest(model = model, data = data, title = title, ...)
        } else if (inherits(model, "glm")) {
            glmforest(model = model, data = data, title = title, ...)
        } else if (inherits(model, "lm")) {
            lmforest(model = model, data = data, title = title, ...)
        } else {
            stop(paste("Model class", model_class, 
                       "is not supported. Supported classes are: lm, glm, coxph, clogit"))
        }
    }
}
