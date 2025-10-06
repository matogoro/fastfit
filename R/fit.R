#' Unified Interface for Model Fitting with Formatted Output
#'
#' Provides a consistent interface for fitting various regression models with
#' automatic formatting of results for publication. Supports multiple model types
#' with a unified syntax.
#'
#' @param data A data.frame or data.table containing the analysis dataset.
#' @param outcome Character string specifying the outcome variable. For survival
#'   analysis, use Surv() syntax (e.g., "Surv(time, status)").
#' @param predictors Character vector of predictor variable names to include
#'   in the model.
#' @param model_type Character string specifying the regression model type:
#'   \itemize{
#'     \item "glm" - Generalized linear model
#'     \item "lm" - Linear model
#'     \item "coxph" - Cox proportional hazards
#'     \item "clogit" - Conditional logistic regression
#'     \item "coxme" - Mixed effects Cox model
#'     \item "glmer" - Generalized linear mixed effects model
#'   }
#'   Default is "glm".
#' @param family For GLM models, the error distribution family. Common options:
#'   "binomial" (logistic), "poisson" (count), "gaussian" (normal), "Gamma".
#'   Default is "binomial". See \code{\link[stats]{family}}.
#' @param interaction_terms Character vector of interaction terms using colon syntax.
#'   Example: c("age:sex", "treatment:stage"). Default is NULL.
#' @param strata For Cox/conditional logistic models, character string naming the
#'   stratification variable. Default is NULL.
#' @param cluster For Cox models, character string naming the variable for robust
#'   clustered standard errors. Default is NULL.
#' @param weights Character string naming the weights variable in data. The weights
#'   column should contain numeric values. Default is NULL.
#' @param conf_level Numeric confidence level for confidence intervals. Must be
#'   between 0 and 1. Default is 0.95 (95 percent CI).
#' @param keep_qc_stats Logical. If TRUE, includes model quality statistics
#'   (AIC, BIC, concordance, etc) in the raw data attribute. Default is TRUE.
#' @param add_reference_rows Logical. If TRUE, adds rows for reference categories
#'   of factor variables. Default is TRUE.
#' @param digits Integer specifying decimal places for effect estimates.
#'   Default is 2.
#' @param digits_p Integer specifying decimal places for p-values.
#'   Default is 3.
#' @param var_labels Named character vector for custom variable labels.
#'   Default is NULL.
#' @param ... Additional arguments passed to the underlying model function
#'   (e.g., subset, na.action, control parameters).
#'
#' @return A data.table with class "fit_result" containing formatted results:
#'   \item{Variable}{Predictor name with custom label if provided}
#'   \item{Group}{Factor level, interaction term, or statistic label}
#'   \item{n}{Sample size}
#'   \item{events}{Number of events (survival/logistic models)}
#'   \item{Multivariable aOR/aHR/aRR (95 percent CI)}{Adjusted effect with CI}
#'   \item{p-value}{Formatted p-value}
#'   
#'   Attributes include:
#'   \item{model}{The fitted model object}
#'   \item{raw_data}{Unformatted numeric results}
#'   \item{formula_str}{The model formula as a string}
#'   \item{model_scope}{Either "Univariable" or "Multivariable"}
#'   \item{model_type}{The regression type used}
#'
#' @details
#' The function automatically detects whether the model is univariable (single
#' predictor) or multivariable (multiple predictors) and labels the output
#' accordingly. For multivariable models, effect estimates are labeled as
#' adjusted (aOR, aHR, aRR).
#' 
#' Interaction terms are specified using colon notation and are included in
#' addition to main effects. For stratified analyses in survival models, the
#' strata variable creates separate baseline hazards for each stratum level.
#' 
#' The formatted output is publication-ready and can be exported directly using
#' tbl2pdf(), tbl2tex(), or tbl2html() functions.
#'
#' @examples
#' \dontrun{
#' # Multivariable logistic regression
#' model1 <- fit(data = mydata,
#'               outcome = "disease",
#'               predictors = c("age", "sex", "smoking"),
#'               model_type = "glm",
#'               family = "binomial")
#' print(model1)
#' 
#' # Cox model with stratification
#' library(survival)
#' cox_model <- fit(data = lung,
#'                  outcome = "Surv(time, status)",
#'                  predictors = c("age", "sex", "ph.ecog"),
#'                  model_type = "coxph",
#'                  strata = "inst")
#' 
#' # Model with interactions and custom labels
#' labels <- c(age = "Age (years)",
#'             bmi = "Body Mass Index",
#'             treatment = "Treatment Group")
#' interact_model <- fit(data = trial_data,
#'                       outcome = "response",
#'                       predictors = c("age", "bmi", "treatment"),
#'                       interaction_terms = c("age:treatment"),
#'                       var_labels = labels)
#' 
#' # Linear regression with weights
#' weighted_model <- fit(data = survey_data,
#'                       outcome = "income",
#'                       predictors = c("education", "experience"),
#'                       model_type = "lm",
#'                       weights = "survey_weight")
#' 
#' # Access the underlying model
#' raw_model <- attr(model1, "model")
#' summary(raw_model)
#' 
#' # Access raw numeric results
#' raw_data <- attr(model1, "raw_data")
#' 
#' # Export formatted results
#' tbl2pdf(model1, "regression_results.pdf")
#' }
#'
#' @seealso
#' \code{\link{uscreen}} for univariable screening,
#' \code{\link{fastfit}} for complete analysis workflow,
#' \code{\link{m2dt}} for model-to-data.table conversion,
#' \code{\link{tbl2pdf}} for PDF export
#'
#' @export
fit <- function(data, 
                outcome,
                predictors,
                model_type = "glm",
                family = "binomial",
                interaction_terms = NULL,
                strata = NULL,
                cluster = NULL,
                weights = NULL,
                conf_level = 0.95,
                keep_qc_stats = TRUE,
                add_reference_rows = TRUE,
                digits = 2,
                digits_p = 3,
                var_labels = NULL,
                ...) {
    
    ## Don't create internal variables - work with 'data' directly
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
    ## Store original data name for reference
    data_name <- deparse(substitute(data))
    
    ## Build formula
    if (!is.null(interaction_terms)) {
        formula_str <- paste(outcome, "~", 
                             paste(c(predictors, interaction_terms), collapse = " + "))
    } else {
        formula_str <- paste(outcome, "~", paste(predictors, collapse = " + "))
    }
    
    ## Add strata if provided (for Cox models)
    if (!is.null(strata) && model_type %in% c("coxph", "clogit")) {
        formula_str <- paste(formula_str, "+ strata(", strata, ")")
    }
    
    formula <- stats::as.formula(formula_str)
    
    ## Fit model based on type - use 'data' directly
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
    
    ## Store the data directly in the model
    model$data <- data
    
    ## Attach metadata as attributes
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
    
    ## Convert to readable format using m2dt
    raw_data <- m2dt(model, 
                     conf_level = conf_level,
                     keep_qc_stats = keep_qc_stats,
                     add_reference_rows = add_reference_rows)

    ## Format results for publication
    formatted <- format_model_table(raw_data,
                                    digits = digits,
                                    digits_p = digits_p,
                                    var_labels = var_labels)

    ## Attach metadata
    setattr(formatted, "model", model)
    setattr(formatted, "raw_data", raw_data)
    setattr(formatted, "outcome", outcome)
    setattr(formatted, "predictors", predictors)
    setattr(formatted, "formula_str", formula_str)
    setattr(formatted, "model_scope", raw_data$model_scope[1])
    setattr(formatted, "model_type", raw_data$model_type[1])
    
    ## Add metadata from model fitting
    if (!is.null(interaction_terms)) {
        setattr(result, "interaction_terms", interaction_terms)
    }
    if (!is.null(strata)) {
        setattr(result, "strata", strata)
    }
    if (!is.null(cluster)) {
        setattr(result, "cluster", cluster)
    }
    if (!is.null(weights)) {
        setattr(result, "weights", weights)
    }

    class(formatted) <- c("fit_result", class(formatted))
    
    return(formatted)
}

#' Print method for fit results
#' @export
print.fit_result <- function(x, ...) {
    cat("\n", attr(x, "model_scope"), " ", attr(x, "model_type"), " Model\n", sep = "")
    cat("Formula: ", attr(x, "formula_str"), "\n", sep = "")
    
    ## Get sample size from raw results
    raw <- attr(x, "raw_data")
    if (!is.null(raw)) {
        cat("N = ", raw$n[1], sep = "")
        if (!is.na(raw$events[1])) cat(", Events = ", raw$events[1], sep = "")
        cat("\n")
    }
    cat("\n")
    
    NextMethod("print", x)
    invisible(x)
}
