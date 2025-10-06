#' Compare Multiple Regression Models
#'
#' Fits multiple regression models and provides a comprehensive comparison table
#' with model quality metrics, convergence diagnostics, and selection guidance.
#'
#' @param data A data.frame or data.table containing the dataset.
#' @param outcome Character string specifying the outcome variable. For survival
#'   analysis, use Surv() syntax (e.g., "Surv(time, status)").
#' @param model_list List of character vectors, each containing predictor names
#'   for one model. Can also be a single character vector to auto-generate nested models.
#' @param model_names Character vector of names for each model. If NULL, uses
#'   "Model 1", "Model 2", etc. Default is NULL.
#' @param model_type Character string specifying model type. If "auto", detects
#'   based on outcome. Options: "auto", "glm", "lm", "coxph", "clogit".
#'   Default is "auto".
#' @param family For GLM models, the error distribution family. Default is "binomial".
#' @param conf_level Numeric confidence level for intervals. Default is 0.95.
#' @param include_coefficients Logical. If TRUE, includes a second table with
#'   coefficient estimates. Default is FALSE.
#' @param scoring_weights Named list of scoring weights. Each weight should be
#'   between 0 and 1, and they should sum to 1. Available metrics depend on model
#'   type. If NULL, uses sensible defaults. See Details for available metrics.
#' @param var_labels Named character vector for custom variable labels. Default is NULL.
#' @param ... Additional arguments passed to model fitting functions.
#'
#' @return A data.table with class "compfit_result" containing:
#'   \item{Model}{Model name/identifier}
#'   \item{N}{Sample size}
#'   \item{Events}{Number of events (for survival/logistic)}
#'   \item{Predictors}{Number of predictors}
#'   \item{Converged}{Whether model converged properly}
#'   \item{AIC}{Akaike Information Criterion}
#'   \item{BIC}{Bayesian Information Criterion}
#'   \item{Pseudo-R^2}{McFadden's pseudo-R-squared (GLM)}
#'   \item{Concordance}{C-statistic (logistic/survival)}
#'   \item{Brier Score}{Brier accuracy score (logistic)}
#'   \item{Global p}{Overall model p-value}
#'   \item{Recommendation}{Scoring metric for model selection}
#'   
#'   Attributes include:
#'   \item{models}{List of fitted model objects}
#'   \item{coefficients}{Coefficient comparison table (if requested)}
#'   \item{best_model}{Name of recommended model}
#'
#' @details
#' The function fits all specified models and computes comprehensive quality
#' metrics for comparison. The composite score combines multiple metrics:
#' lower AIC/BIC (information criteria), higher concordance (discrimination),
#' and model convergence status.
#' 
#' For GLMs, McFadden's pseudo-R-squared is calculated as 1 - (logLik/logLik_null).
#' For survival models, the global p-value comes from the log-rank test.
#' 
#' Models that fail to converge are flagged and penalized in the composite score.
#' 
#' Scoring weights can be customized based on model type:
#' \itemize{
#'   \item GLM: "convergence", "aic", "concordance", "pseudo_r2", "brier" 
#'   \item Cox: "convergence", "aic", "concordance", "global_p"
#'   \item Linear: "convergence", "aic", "pseudo_r2", "rmse"
#' }
#' Default weights emphasize discrimination (concordance) and model fit (AIC).
#'
#' @examples
#' \dontrun{
#' # Compare nested models
#' models <- list(
#'   base = c("age", "sex"),
#'   clinical = c("age", "sex", "bmi", "smoking"),
#'   full = c("age", "sex", "bmi", "smoking", "diabetes", "hypertension")
#' )
#' 
#' comparison <- compfit(data = mydata,
#'                      outcome = "disease",
#'                      model_list = models,
#'                      model_names = c("Base", "Clinical", "Full"))
#' print(comparison)
#' 
#' # Auto-detect model type for survival outcome
#' surv_models <- list(
#'   simple = c("age", "stage"),
#'   complex = c("age", "stage", "grade", "treatment")
#' )
#' 
#' surv_comp <- compfit(lung,
#'                     outcome = "Surv(time, status)",
#'                     model_list = surv_models,
#'                     model_type = "auto")
#' 
#' # Include coefficient comparison
#' detailed <- compfit(data, outcome, models,
#'                    include_coefficients = TRUE)
#' coef_table <- attr(detailed, "coefficients")
#' 
#' # Export comparison table
#' tbl2pdf(comparison, "model_comparison.pdf")
#' }
#'
#' @seealso
#' \code{\link{fit}} for individual model fitting,
#' \code{\link{fastfit}} for automated variable selection,
#' \code{\link{tbl2pdf}} for exporting results
#'
#' @export
compfit <- function(data,
                    outcome,
                    model_list,
                    model_names = NULL,
                    model_type = "auto",
                    family = "binomial",
                    conf_level = 0.95,
                    include_coefficients = FALSE,
                    scoring_weights = NULL,
                    var_labels = NULL,
                    ...) {
    
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }

    ## Auto-detect model type if requested
    if (model_type == "auto") {
        if (grepl("^Surv\\(", outcome)) {
            model_type <- "coxph"
            message("Auto-detected survival outcome, using Cox regression")
        } else if (is.factor(data[[outcome]]) || length(unique(data[[outcome]])) == 2) {
            model_type <- "glm"
            family <- "binomial"
            message("Auto-detected binary outcome, using logistic regression")
        } else if (is.numeric(data[[outcome]])) {
            model_type <- "lm"
            message("Auto-detected continuous outcome, using linear regression")
        } else {
            stop("Cannot auto-detect model type for outcome: ", outcome)
        }
    }

    ## Set model names if not provided
    if (is.null(model_names)) {
        model_names <- paste("Model", seq_along(model_list))
    }

    ## Initialize results storage
    comparison <- data.table::data.table()
    models <- list()
    coef_results <- list()

    ## Fit each model
    for (i in seq_along(model_list)) {
        message(sprintf("Fitting %s with %d predictors...", 
                        model_names[i], length(model_list[[i]])))
        
        ## Fit model using fit() for consistency
        model_result <- tryCatch({
            fit(data = data,
                outcome = outcome,
                predictors = model_list[[i]],
                model_type = model_type,
                family = family,
                conf_level = conf_level,
                var_labels = var_labels,
                keep_qc_stats = TRUE,
                ...)
        }, error = function(e) {
            message("  Warning - model failed to fit: ", e$message)
            NULL
        })

        ## Build comparison row - conditional on model type
        if (is.null(model_result)) {
            ## Model failed to fit
            row <- data.table::data.table(
                                   Model = model_names[i],
                                   N = nrow(data),
                                   Events = NA_integer_,
                                   Predictors = length(model_list[[i]]),
                                   Converged = "Failed",
                                   AIC = NA_real_,
                                   BIC = NA_real_,
                                   `Pseudo-R^2` = NA_real_,
                                   Concordance = NA_real_,
                                   `Global p` = NA_real_
                               )
            ## Add Brier Score only for GLM
            if (model_type == "glm") {
                row$`Brier Score` <- NA_real_
            }
        } else {

            ## Extract model and raw data
            model <- attr(model_result, "model")
            raw_data <- attr(model_result, "raw_data")
            models[[model_names[i]]] <- model
            
            ## Store coefficients if requested
            if (include_coefficients) {
                coef_results[[model_names[i]]] <- model_result
            }
            
            ## Check convergence
            converged <- check_convergence(model)
            
            ## Extract metrics
            metrics <- extract_model_metrics(model, raw_data, model_type)

            ## Build comparison row
            row <- data.table::data.table(
                                   Model = model_names[i],
                                   N = metrics$n,
                                   Events = metrics$events,
                                   Predictors = length(model_list[[i]]),
                                   Converged = converged,
                                   AIC = if (!is.null(metrics$aic) && !is.na(metrics$aic)) round(metrics$aic, 1) else NA_real_,
                                   BIC = if (!is.null(metrics$bic) && !is.na(metrics$bic)) round(metrics$bic, 1) else NA_real_,
                                   `Pseudo-R^2` = if (!is.null(metrics$pseudo_r2) && !is.na(metrics$pseudo_r2)) round(metrics$pseudo_r2, 3) else NA_real_,
                                   Concordance = if (!is.null(metrics$concordance) && !is.na(metrics$concordance)) round(metrics$concordance, 3) else NA_real_,
                                   `Global p` = if (!is.null(metrics$global_p) && !is.na(metrics$global_p)) format_pvalue(metrics$global_p, 3) else NA_character_
                               )
            
            ## Add Brier Score only for GLM
            if (model_type == "glm") {
                row$`Brier Score` <- if (!is.null(metrics$brier_score)) round(metrics$brier_score, 3) else NA_real_
            }
        }
        
        comparison <- rbind(comparison, row, fill = TRUE)
    }

    ## Calculate scores and sort
    comparison <- calculate_model_scores(comparison, model_type, scoring_weights)

    ## Format for display (ensure proper column order)
    if (model_type == "glm") {
        setcolorder(comparison, c("Model", "N", "Events", "Predictors", "Converged",
                                  "AIC", "BIC", "Pseudo-R^2", "Concordance", "Brier Score",
                                  "Global p", "Composite Score"))
    } else {
        setcolorder(comparison, c("Model", "N", "Events", "Predictors", "Converged",
                                  "AIC", "BIC", "Pseudo-R^2", "Concordance",
                                  "Global p", "Composite Score"))
    }
    ## Attach attributes
    setattr(comparison, "models", models)
    setattr(comparison, "model_type", model_type)
    setattr(comparison, "outcome", outcome)

    if (include_coefficients && length(coef_results) > 0) {
        coef_table <- combine_coefficient_tables(coef_results, model_names)
        setattr(comparison, "coefficients", coef_table)
    }

    ## Identify best model (first one after sorting)
    if (nrow(comparison) > 0) {
        setattr(comparison, "best_model", comparison$Model[1])
    }

    class(comparison) <- c("compfit_result", class(comparison))

    return(comparison)
}
