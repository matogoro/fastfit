#' Compare Multiple Models
#'
#' Fits multiple multivariable models with different predictor sets and combines
#' results for easy comparison.
#'
#' @param data Data.frame or data.table containing the dataset.
#' @param outcome Character string specifying the outcome variable. For survival 
#'   analysis, use Surv() syntax, e.g., "Surv(time, status)".
#' @param model_list List of character vectors, each containing predictor names 
#'   for one model.
#' @param model_names Character vector of names for each model. If NULL, uses 
#'   "Model 1", "Model 2", etc.
#' @param model_type Character string: "glm", "lm", "coxph", "clogit". Default "glm".
#' @param family For GLM models, the error distribution family. Default "binomial".
#' @param conf_level Numeric confidence level. Default 0.95.
#' @param keep_qc_stats Logical. Include model quality statistics. Default TRUE.
#' @param add_reference_rows Logical. Add reference category rows. Default TRUE.
#' @param ... Additional arguments passed to model fitting functions.
#'
#' @return A data.table combining results from all models, with the variable 
#'   column indicating which model each row belongs to.
#'
#' @details
#' This function is useful for comparing nested models, different variable 
#' selection strategies, or alternative model specifications. All models use
#' the same outcome and model type for valid comparison.
#' 
#' For automated univariable-to-multivariable workflows with variable selection,
#' use \code{fastfit()} instead. This function is designed for comparing
#' pre-specified model specifications.
#'
#' @examples
#' # Compare nested models
#' models <- list(
#'   demographics = c("age", "sex"),
#'   clinical = c("age", "sex", "bmi", "bp"),
#'   full = c("age", "sex", "bmi", "bp", "smoking", "diabetes")
#' )
#' 
#' comparison <- cmodels(data, "outcome", models,
#'                          model_names = c("Demographics", 
#'                                         "Clinical", 
#'                                         "Full Model"))
#' 
#' # View AIC across models
#' comparison[, .(AIC = unique(AIC)), by = variable]
#' 
#' # Create forest plot comparing models
#' forest_plot(comparison, facet_by = "variable")
#' 
#' # For screening-based selection, use fastfit instead:
#' # auto_model <- fastfit(data, "outcome", all_predictors, method = "screen")
#'
#' @seealso \code{\link{mmodel}} for fitting individual models, 
#'   \code{\link{fastfit}} for automated variable selection workflows
#' @export
cmodels <- function(data, outcome, model_list,
                    model_names = NULL,
                       model_type = "glm",
                       family = "binomial",
                       conf_level = 0.95,
                       keep_qc_stats = TRUE,
                       add_reference_rows = TRUE,
                       ...) {
    
    if (is.null(model_names)) {
        model_names <- paste("Model", seq_along(model_list))
    }
    
                                        # Fit all models and combine results
    results <- data.table::rbindlist(lapply(seq_along(model_list), function(i) {
        
                                        # Fit model using mmodel
        model <- mmodel(
            data = data,
            outcome = outcome,
            predictors = model_list[[i]],
            model_type = model_type,
            family = family,
            ...
        )
        
                                        # Convert to data.table
        dt <- m2dt(model,
                   conf_level = conf_level,
                   variable_name = model_names[i],
                   keep_qc_stats = keep_qc_stats,
                   add_reference_rows = add_reference_rows)
        
        return(dt)
    }), fill = TRUE)
    
    return(results)
}
