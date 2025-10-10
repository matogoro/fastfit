#' Univariable Screening with Formatted Output
#'
#' Performs univariable regression analyses for multiple predictors against a single
#' outcome, returning publication-ready formatted results. Each predictor is tested
#' independently in its own model.
#'
#' @param data A data.frame or data.table containing the analysis dataset.
#' @param outcome Character string specifying the outcome variable name. For survival
#'   analysis, use Surv() syntax (e.g., "Surv(time, status)").
#' @param predictors Character vector of predictor variable names to screen.
#' @param model_type Character string specifying the regression model type.
#'   Options: "glm" (generalized linear model), "lm" (linear model),
#'   "coxph" (Cox proportional hazards), "clogit" (conditional logistic).
#'   Default is "glm".
#' @param family For GLM models, the error distribution family. Options include
#'   "binomial" (logistic regression), "poisson" (count data), "gaussian" (normal),
#'   "Gamma", etc. See \code{\link[stats]{family}}. Default is "binomial".
#' @param p_threshold Numeric value between 0 and 1 for filtering results by p-value.
#'   Only predictors with p <= threshold are returned. Default is 1 (no filtering).
#' @param conf_level Numeric confidence level for confidence intervals. Must be
#'   between 0 and 1. Default is 0.95 (95 percent CI).
#' @param add_reference_rows Logical. If TRUE, adds rows for reference categories
#'   of factor variables with baseline values (OR/HR = 1). Default is TRUE.
#' @param show_n_events Character vector specifying which optional columns to display.
#'   Options: "n", "events" (or "Events"). Default is c("n", "events") for 
#'   survival/logistic models, "n" only for other models. Set to NULL to hide
#'   these columns entirely.
#' @param digits Integer specifying decimal places for effect estimates (OR, HR, etc).
#'   Default is 2.
#' @param digits_p Integer specifying decimal places for p-values. Values less than
#'   10^(-digits_p) display as "< 0.001" etc. Default is 3.
#' @param var_labels Named character vector for custom variable labels. Names should
#'   match predictor names, values are display labels. Default is NULL.
#' @param keep_models Logical. If TRUE, stores all fitted model objects in the output.
#'   This can consume significant memory for large datasets or many predictors.
#'   Models are accessible via attr(result, "models"). Default is FALSE.
#' @param exponentiate Logical. Whether to exponentiate coefficients. Default is NULL,
#'   which automatically displays exponentiated coefficients for logistic/Poisson/Cox
#'   regression models and raw coefficients for log-link or linear regression models.
#' @param ... Additional arguments passed to the underlying model fitting functions
#'   (e.g., weights, subset, na.action).
#'
#' @return A data.table with class "uscreen_result" containing formatted results:
#'   \item{Variable}{Predictor name (or custom label if provided)}
#'   \item{Group}{Factor level or statistic type for continuous variables}
#'   \item{n}{Sample size}
#'   \item{events}{Number of events (for survival/logistic models)}
#'   \item{Univariable OR/HR/RR (95 percent CI)}{Formatted effect size with confidence interval}
#'   \item{p-value}{Formatted p-value}
#'   
#'   The returned object includes attributes:
#'   \item{raw_data}{Unformatted numeric results for further analysis}
#'   \item{models}{List of fitted model objects (if keep_models = TRUE)}
#'   \item{outcome}{The outcome variable name}
#'   \item{model_type}{The type of regression model used}
#'
#' @details
#' The function iterates through each predictor, fitting a separate univariable
#' model for each. This is useful for:
#' \itemize{
#'   \item Initial variable screening before multivariable modeling
#'   \item Understanding crude (unadjusted) associations
#'   \item Identifying multicollinearity issues
#'   \item Variable selection for further analysis
#' }
#' 
#' For factor variables with add_reference_rows = TRUE, the reference category
#' is shown with OR/HR = 1.00 (Reference) and no p-value. P-values are displayed
#' only for non-reference categories.
#' 
#' The formatted output is ready for export via tbl2pdf(), tbl2tex(), or tbl2html().
#'
#' @examples
#' \dontrun{
#' # Basic logistic regression screening
#' data(mtcars)
#' results <- uscreen(mtcars, 
#'                    outcome = "am",
#'                    predictors = c("mpg", "cyl", "disp", "hp"),
#'                    model_type = "glm",
#'                    family = "binomial")
#' print(results)
#' 
#' # Cox regression with custom labels
#' library(survival)
#' data(lung)
#' labels <- c(age = "Age (years)", 
#'             sex = "Sex", 
#'             ph.ecog = "ECOG Score")
#' cox_screen <- uscreen(lung,
#'                       outcome = "Surv(time, status)",
#'                       predictors = c("age", "sex", "ph.ecog"),
#'                       model_type = "coxph",
#'                       var_labels = labels)
#' 
#' # Filter by p-value threshold
#' significant <- uscreen(mydata,
#'                        outcome = "disease",
#'                        predictors = c("var1", "var2", "var3"),
#'                        p_threshold = 0.05)
#' 
#' # Keep models for diagnostics
#' with_models <- uscreen(mydata,
#'                        outcome = "outcome",
#'                        predictors = vars,
#'                        keep_models = TRUE)
#' models <- attr(with_models, "models")
#' plot(models[["age"]])  # Diagnostic plots
#' 
#' # Export to PDF
#' tbl2pdf(results, "screening_results.pdf")
#' }
#'
#' @seealso 
#' \code{\link{fit}} for multivariable modeling,
#' \code{\link{fastfit}} for complete univariable-to-multivariable workflow,
#' \code{\link{tbl2pdf}} for exporting results
#'
#' @export
uscreen <- function(data,
                    outcome,
                    predictors,
                    model_type = "glm",
                    family = "binomial",
                    p_threshold = 1,
                    conf_level = 0.95,
                    add_reference_rows = TRUE,
                    show_n_events = c("n", "events"),
                    digits = 2,
                    digits_p = 3,
                    var_labels = NULL,
                    keep_models = FALSE,
                    exponentiate = NULL,
                    ...) {
    
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
    ## Store models only if requested
    if (keep_models) {
        models <- list()
    }
    raw_results <- list()
    
    ## Fit univariable model for each predictor
    for (pred in predictors) {
        ## Build formula
        formula_str <- paste(outcome, "~", pred)
        formula <- stats::as.formula(formula_str)
        
        ## Fit model based on type
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
                stop("Package 'survival' required for conditional logistic regression")
            model <- survival::clogit(formula, data = data, ...)
        } else {
            stop("Unsupported model type: ", model_type)
        }

        ## Store the data directly in the model
        model$data <- data
        
        ## Store the model only if requested
        if (keep_models) {
            models[[pred]] <- model
        }
        
        ## Get raw results using m2dt
        raw_result <- m2dt(model,
                           conf_level = conf_level,
                           keep_qc_stats = FALSE,
                           add_reference_rows = add_reference_rows)
        
        ## Add predictor name for tracking
        raw_result[, predictor := pred]
        raw_results[[pred]] <- raw_result
        
        ## Clear model from memory if not keeping
        if (!keep_models) {
            rm(model)
        }
    }
    
    ## Combine all raw results
    combined_raw <- rbindlist(raw_results, fill = TRUE)
    
    ## Filter by p-value if requested
    if (p_threshold < 1) {
        passing_predictors <- unique(combined_raw[p_value <= p_threshold]$predictor)
        combined_raw <- combined_raw[predictor %in% passing_predictors]
        
        ## If filtering and keeping models, remove non-passing models
        if (keep_models) {
            models <- models[names(models) %in% passing_predictors]
        }
    }
    
    ## Format the combined results
    formatted <- format_model_table(combined_raw,
                                    show_n_events = show_n_events,
                                    digits = digits,
                                    digits_p = digits_p,
                                    var_labels = var_labels,
                                    exponentiate = exponentiate
                                    )
    
    ## Attach raw data
    setattr(formatted, "raw_data", combined_raw)
    
    ## Attach models only if kept
    if (keep_models) {
        setattr(formatted, "models", models)
    }
    
    setattr(formatted, "outcome", outcome)
    setattr(formatted, "model_type", unique(combined_raw$model_type)[1])
    setattr(formatted, "model_scope", "Univariable")
    setattr(formatted, "screening_type", "univariable")
    
    ## Add class for methods
    class(formatted) <- c("uscreen_result", class(formatted))
    
    return(formatted)
}

#' Print method for uscreen results
#' @keywords internal
#' @export
print.uscreen_result <- function(x, ...) {
    cat("\nUnivariable Screening Results\n")
    cat("Outcome: ", attr(x, "outcome"), "\n", sep = "")
    cat("Model Type: ", attr(x, "model_type"), "\n", sep = "")
    
    ## Get unique predictors from the formatted table
    n_predictors <- length(unique(x$Variable[x$Variable != ""]))
    cat("Predictors Screened: ", n_predictors, "\n", sep = "")
    
    ## Count significant from raw data
    raw <- attr(x, "raw_data")
    if (!is.null(raw) && "p_value" %in% names(raw)) {
        sig_predictors <- unique(raw[p_value < 0.05]$predictor)
        n_sig <- length(sig_predictors)
        cat("Significant (p < 0.05): ", n_sig, "\n", sep = "")
    }
    
    ## Note if models are stored
    if (!is.null(attr(x, "models"))) {
        cat("Models stored: Yes (", length(attr(x, "models")), ")\n", sep = "")
    }
    
    cat("\n")
    
    NextMethod("print", x)
    invisible(x)
}
