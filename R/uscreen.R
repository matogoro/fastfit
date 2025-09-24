#' Univariable Screening Function
#'
#' Performs univariable regression analysis for multiple predictors against a 
#' single outcome, returning standardized results in data.table format.
#'
#' @param data A data.frame or data.table containing the analysis dataset.
#' @param outcome Character string specifying the outcome variable. For survival 
#'   analysis, use Surv() syntax, e.g., "Surv(time, status)".
#' @param predictors Character vector of predictor variable names to screen.
#' @param model_type Character string specifying model type. Options: "glm", "lm", 
#'   "coxph", "clogit". Default is "glm".
#' @param family For GLM models, the error distribution family. Default is "binomial".
#' @param conf_level Numeric between 0 and 1 for confidence intervals. Default is 0.95.
#' @param keep_models Logical. Whether to store fitted model objects in results. 
#'   Default is FALSE to save memory.
#' @param keep_qc_stats Logical. Include model quality statistics. Default is TRUE.
#' @param sort_by Character string naming column to sort by, or NULL to maintain 
#'   original predictor order. Options include "p_value", "effect", etc. Default is NULL.
#' @param add_reference_rows Logical. Add reference category rows for factors. 
#'   Default is TRUE.
#' @param ... Additional arguments passed to model fitting functions.
#'
#' @return A data.table with one row per coefficient containing:
#' \itemize{
#'   \item All columns from \code{m2dt()} output
#'   \item \code{.var_order}: Integer. Original position in predictors vector
#'   \item \code{model}: List column with model objects (if keep_models = TRUE)
#' }
#'
#' @details
#' This function fits separate univariable models for each predictor and combines
#' results into a single table for easy comparison and selection of variables for
#' multivariable models.
#'
#' The output preserves the original order of predictors unless sort_by is specified.
#' This is useful for maintaining logical groupings of variables.
#'
#' @examples
#' # Basic univariable screening
#' predictors <- c("age", "sex", "bmi", "smoking")
#' results <- uscreen(data, "outcome", predictors)
#' 
#' # View significant predictors only
#' results[sig_binary == TRUE]
#' 
#' # Cox model screening
#' cox_results <- uscreen(lung, "Surv(time, status)", 
#'                        c("age", "sex", "ph.ecog"),
#'                        model_type = "coxph")
#' 
#' # Sort by p-value
#' results_sorted <- uscreen(data, "outcome", predictors, 
#'                           sort_by = "p_value")
#'
#' @seealso \code{\link{m2dt}}, \code{\link{mmodel}}, \code{\link{u2m}}
#' @export
uscreen <- function(data, outcome, predictors, 
                    model_type = "glm",
                    family = "binomial",
                    conf_level = 0.95, 
                    keep_models = FALSE,
                    keep_qc_stats = TRUE,
                    sort_by = NULL,
                    add_reference_rows = TRUE,
                    ...) {
    
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
                                        # Change here: Add fill=TRUE to rbindlist
    results <- data.table::rbindlist(lapply(seq_along(predictors), function(i) {
        var <- predictors[i]
        
                                        # Build formula
        formula <- stats::as.formula(paste(outcome, "~", var))
        
                                        # Fit model based on type
        if (model_type == "glm") {
            model <- stats::glm(formula, data = data, family = family, ...)
            
        } else if (model_type == "lm") {
            model <- stats::lm(formula, data = data, ...)
            
        } else if (model_type == "coxph") {
            if (!requireNamespace("survival", quietly = TRUE)) 
                stop("Package 'survival' required for Cox models")
            model <- survival::coxph(formula, data = data, ...)
            
        } else if (model_type == "clogit") {
            if (!requireNamespace("survival", quietly = TRUE)) 
                stop("Package 'survival' required for conditional logistic")
            model <- survival::clogit(formula, data = data, ...)
            
        } else {
            stop("Unsupported model type: ", model_type)
        }
        
                                        # Extract results using m2dt
        dt <- m2dt(model, 
                   conf_level = conf_level,
                   variable_name = var,
                   keep_qc_stats = keep_qc_stats,
                   add_reference_rows = add_reference_rows)
        
                                        # Add original order index
        dt[, .var_order := i]
        
                                        # Optionally keep model object
        if (keep_models) {
            dt[, model := list(list(model))]
        }
        
        return(dt)
    }), fill = TRUE)  # KEY CHANGE: Add fill=TRUE here
    
                                        # Sort by original variable order (default)
    data.table::setorder(results, .var_order)
    
                                        # Apply additional sorting if requested
    if (!is.null(sort_by) && sort_by %in% names(results)) {
        data.table::setorder(results, cols = sort_by)
    }
    
    return(results)
}
