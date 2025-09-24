#' Fit Multivariable Model
#'
#' A unified interface for fitting various types of regression models with 
#' consistent syntax across different model types (GLM, Cox, mixed effects, etc.)
#'
#' @param data A data.frame or data.table containing the analysis dataset.
#' @param outcome Character string specifying the outcome variable. For survival 
#'   analysis, use Surv() syntax, e.g., "Surv(time, status)".
#' @param predictors Character vector of predictor variable names.
#' @param model_type Character string specifying model type. Options: "glm", "lm", 
#'   "coxph", "clogit", "coxme", "glmer". Default is "glm".
#' @param family For GLM models, the error distribution family. Default is "binomial". 
#'   See \code{\link[stats]{family}} for options.
#' @param interaction_terms Character vector of interaction terms using colon syntax, 
#'   e.g., c("age:sex", "bmi:smoking"). Default is NULL.
#' @param strata For Cox/conditional logistic models, variable name for stratification. 
#'   Default is NULL.
#' @param cluster For Cox models, variable name for robust clustered standard errors. 
#'   Default is NULL.
#' @param weights Character string naming the weights variable in data. Default is NULL.
#' @param ... Additional arguments passed to the underlying model fitting function.
#'
#' @return A fitted model object with additional attributes:
#' \itemize{
#'   \item \code{formula_str}: Character. The formula used as a string
#'   \item \code{predictors}: Character vector. The predictor variables
#'   \item \code{model_type}: Character. The type of model fitted
#'   \item \code{interaction_terms}: Character vector. Interaction terms (if any)
#'   \item \code{data_name}: Character. Name of the data object used
#' }
#' 
#' @export
mmodel <- function(data, 
                   outcome, 
                   predictors,
                   model_type = "glm",
                   family = "binomial",
                   interaction_terms = NULL,
                   strata = NULL,
                   cluster = NULL,
                   weights = NULL,
                   ...) {
    
                                        # Don't create internal variables - work with 'data' directly
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
                                        # Store original data name for reference
    data_name <- deparse(substitute(data))
    
                                        # Build formula
    if (!is.null(interaction_terms)) {
        formula_str <- paste(outcome, "~", 
                             paste(c(predictors, interaction_terms), collapse = " + "))
    } else {
        formula_str <- paste(outcome, "~", paste(predictors, collapse = " + "))
    }
    
                                        # Add strata if provided (for Cox models)
    if (!is.null(strata) && model_type %in% c("coxph", "clogit")) {
        formula_str <- paste(formula_str, "+ strata(", strata, ")")
    }
    
    formula <- stats::as.formula(formula_str)
    
                                        # Fit model based on type - use 'data' directly
    if (model_type == "glm") {
        if (!is.null(weights)) {
            model <- stats::glm(formula, data = data, family = family, 
                                weights = data[[weights]], ...)
        } else {
            model <- stats::glm(formula, data = data, family = family, ...)
        }
        
    } else if (model_type == "lm") {
        if (!is.null(weights)) {
            model <- stats::lm(formula, data = data, weights = data[[weights]], ...)
        } else {
            model <- stats::lm(formula, data = data, ...)
        }
        
    } else if (model_type == "coxph") {
        if (!requireNamespace("survival", quietly = TRUE)) 
            stop("Package 'survival' required for Cox models")
        
        if (!is.null(cluster)) {
            model <- survival::coxph(formula, data = data, 
                                     cluster = data[[cluster]], ...)
        } else {
            model <- survival::coxph(formula, data = data, ...)
        }
        
    } else if (model_type == "clogit") {
        if (!requireNamespace("survival", quietly = TRUE)) 
            stop("Package 'survival' required for conditional logistic regression")
        model <- survival::clogit(formula, data = data, ...)
        
    } else if (model_type == "coxme") {
        if (!requireNamespace("coxme", quietly = TRUE))
            stop("Package 'coxme' required for mixed effects Cox models")
        model <- coxme::coxme(formula, data = data, ...)
        
    } else if (model_type == "glmer") {
        if (!requireNamespace("lme4", quietly = TRUE))
            stop("Package 'lme4' required for mixed effects models")
        model <- lme4::glmer(formula, data = data, family = family, ...)
        
    } else {
        stop("Unsupported model type: ", model_type)
    }
    
                                        # Store the data directly in the model
    model$data <- data
    
                                        # Attach metadata as attributes
    data.table::setattr(model, "formula_str", formula_str)
    data.table::setattr(model, "predictors", predictors)
    data.table::setattr(model, "model_type", model_type)
    data.table::setattr(model, "data_name", data_name)
    
    if (!is.null(interaction_terms)) {
        data.table::setattr(model, "interaction_terms", interaction_terms)
    }
    if (!is.null(strata)) {
        data.table::setattr(model, "strata", strata)
    }
    if (!is.null(cluster)) {
        data.table::setattr(model, "cluster", cluster)
    }
    if (!is.null(weights)) {
        data.table::setattr(model, "weights", weights)
    }
    
                                        # Add class for method dispatch
    class(model) <- c("mmodel", class(model))
    
    return(model)
}
