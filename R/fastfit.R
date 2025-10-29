#' Complete Regression Analysis Workflow
#'
#' Executes a comprehensive regression analysis pipeline that combines univariable 
#' screening, automatic or manual variable selection, and multivariable modeling 
#' in a single function call. This function is designed to streamline the complete 
#' analytical workflow from initial exploration to final adjusted models, with 
#' publication-ready formatted output showing both univariable and multivariable 
#' results side-by-side.
#'
#' @param data Data.frame or data.table containing the analysis dataset. The 
#'   function automatically converts data.frames to data.tables for processing.
#'   
#' @param outcome Character string specifying the outcome variable name. For 
#'   survival analysis, use \code{Surv()} syntax (e.g., \code{"Surv(os_months, os_status)"}).
#'   
#' @param predictors Character vector of predictor variable names to analyze. 
#'   All predictors are tested in univariable models. The subset included in 
#'   the multivariable model depends on the \code{method} parameter.
#'   
#' @param method Character string specifying the variable selection strategy:
#'   \itemize{
#'     \item \code{"screen"} - Automatic selection based on univariable p-value 
#'       threshold. Only predictors with p ≤ \code{p_threshold} in univariable 
#'       analysis are included in the multivariable model [default]
#'     \item \code{"all"} - Include all predictors in both univariable and 
#'       multivariable analyses (no selection)
#'     \item \code{"custom"} - Manual selection. All predictors in univariable 
#'       analysis, but only those specified in \code{multi_predictors} are 
#'       included in multivariable model
#'   }
#'   
#' @param multi_predictors Character vector of predictors to include in the 
#'   multivariable model when \code{method = "custom"}. Required when using 
#'   custom selection. Ignored for other methods. Default is \code{NULL}.
#'   
#' @param p_threshold Numeric p-value threshold for automatic variable selection 
#'   when \code{method = "screen"}. Predictors with univariable p-value ≤ 
#'   threshold are included in multivariable model. Common values: 0.05 (strict), 
#'   0.10 (moderate), 0.20 (liberal screening). Default is 0.05. Ignored for 
#'   other methods.
#'   
#' @param columns Character string specifying which result columns to display:
#'   \itemize{
#'     \item \code{"both"} - Show both univariable and multivariable results 
#'       side-by-side [default]
#'     \item \code{"uni"} - Show only univariable results
#'     \item \code{"multi"} - Show only multivariable results
#'   }
#'   
#' @param model_type Character string specifying the regression model type:
#'   \itemize{
#'     \item \code{"glm"} - Generalized linear model [default]
#'     \item \code{"lm"} - Linear regression
#'     \item \code{"coxph"} - Cox proportional hazards
#'     \item \code{"clogit"} - Conditional logistic regression
#'   }
#'   
#' @param family For GLM models, the error distribution and link function. 
#'   Common options: \code{"binomial"} (logistic) [default], \code{"poisson"} 
#'   (count data), \code{"gaussian"} (linear), \code{"Gamma"}. See 
#'   \code{\link[stats]{family}}.
#'   
#' @param conf_level Numeric confidence level for confidence intervals. Must be 
#'   between 0 and 1. Default is 0.95 (95\% CI).
#'   
#' @param add_reference_rows Logical. If \code{TRUE}, adds rows for reference 
#'   categories of factor variables with baseline values. Default is \code{TRUE}.
#'   
#' @param show_n Logical. If \code{TRUE}, includes sample size columns. 
#'   Default is \code{TRUE}.
#'   
#' @param show_events Logical. If \code{TRUE}, includes events columns (survival 
#'   and logistic models). Default is \code{TRUE}.
#'   
#' @param digits Integer specifying decimal places for effect estimates. 
#'   Default is 2.
#'   
#' @param digits_p Integer specifying decimal places for p-values. Default is 3.
#'   
#' @param var_labels Named character vector or list for custom variable display 
#'   labels. Default is \code{NULL}.
#'   
#' @param metrics Character specification for which statistics to display:
#'   \itemize{
#'     \item \code{"both"} - Show effect estimates with CI and p-values [default]
#'     \item \code{"effect"} - Show only effect estimates with CI
#'     \item \code{"p"} - Show only p-values
#'   }
#'   Can also be a character vector: \code{c("effect", "p")} is equivalent to 
#'   \code{"both"}.
#'   
#' @param return_type Character string specifying what to return:
#'   \itemize{
#'     \item \code{"table"} - Return formatted results table only [default]
#'     \item \code{"model"} - Return multivariable model object only
#'     \item \code{"both"} - Return list with both table and model
#'   }
#'   
#' @param keep_models Logical. If \code{TRUE}, stores univariable model objects 
#'   in the output. Can consume significant memory for many predictors. 
#'   Default is \code{FALSE}.
#'   
#' @param exponentiate Logical. Whether to exponentiate coefficients. Default 
#'   is \code{NULL} for automatic selection based on model type.
#'   
#' @param ... Additional arguments passed to model fitting functions (e.g., 
#'   \code{weights}, \code{subset}, \code{na.action}).
#'
#' @return Depends on \code{return_type} parameter:
#'   
#'   When \code{return_type = "table"} (default): A data.table with S3 class 
#'   \code{"fastfit_result"} containing:
#'   \describe{
#'     \item{Variable}{Character. Predictor name or custom label}
#'     \item{Group}{Character. Category level for factors, empty for continuous}
#'     \item{n/n_group}{Integer. Sample sizes (if \code{show_n = TRUE})}
#'     \item{events/events_group}{Integer. Event counts (if \code{show_events = TRUE})}
#'     \item{Univariable OR/HR/RR (95\% CI)}{Character. Unadjusted effect 
#'       (if \code{columns} includes "uni" and \code{metrics} includes "effect")}
#'     \item{Uni p}{Character. Univariable p-value (if \code{columns} includes 
#'       "uni" and \code{metrics} includes "p")}
#'     \item{Multivariable aOR/aHR/aRR (95\% CI)}{Character. Adjusted effect 
#'       (if \code{columns} includes "multi" and \code{metrics} includes "effect")}
#'     \item{Multi p}{Character. Multivariable p-value (if \code{columns} 
#'       includes "multi" and \code{metrics} includes "p")}
#'   }
#'   
#'   When \code{return_type = "model"}: The fitted multivariable model object 
#'   (glm, lm, coxph, etc.).
#'   
#'   When \code{return_type = "both"}: A list with two elements:
#'   \describe{
#'     \item{table}{The formatted results data.table}
#'     \item{model}{The fitted multivariable model object}
#'   }
#'   
#'   The table includes the following attributes:
#'   \describe{
#'     \item{outcome}{Character. The outcome variable name}
#'     \item{model_type}{Character. The regression model type}
#'     \item{method}{Character. The variable selection method used}
#'     \item{columns}{Character. Which columns were displayed}
#'     \item{model}{The multivariable model object (if fitted)}
#'     \item{uni_results}{The complete univariable screening results}
#'     \item{n_multi}{Integer. Number of predictors in multivariable model}
#'   }
#'
#' @details
#' \strong{Analysis Workflow:}
#' 
#' The function implements a complete regression analysis pipeline:
#' \enumerate{
#'   \item \strong{Univariable screening}: Fits separate models for each 
#'     predictor (outcome ~ predictor). Each predictor is tested independently 
#'     to understand crude associations.
#'   \item \strong{Variable selection}: Based on the \code{method} parameter:
#'     \itemize{
#'       \item \code{"screen"}: Automatically selects predictors with univariable 
#'         p ≤ \code{p_threshold}
#'       \item \code{"all"}: Includes all predictors (no selection)
#'       \item \code{"custom"}: Uses predictors specified in \code{multi_predictors}
#'     }
#'   \item \strong{Multivariable modeling}: Fits a single model with selected 
#'     predictors (outcome ~ predictor1 + predictor2 + ...). Estimates are 
#'     adjusted for all other variables in the model.
#'   \item \strong{Output formatting}: Combines results into publication-ready 
#'     table with appropriate effect measures and formatting.
#' }
#' 
#' \strong{Variable Selection Strategies:}
#' 
#' \emph{Screening Method} (\code{method = "screen"}):
#' \itemize{
#'   \item Uses p-value threshold for automatic selection
#'   \item Liberal thresholds (e.g., 0.20) cast a wide net to avoid missing 
#'     important predictors
#'   \item Stricter thresholds (e.g., 0.05) focus on strongly associated predictors
#'   \item Helps reduce overfitting and multicollinearity
#'   \item Common in exploratory analyses and when sample size is limited
#' }
#' 
#' \emph{All Method} (\code{method = "all"}):
#' \itemize{
#'   \item No variable selection - includes all predictors
#'   \item Appropriate when all variables are theoretically important
#'   \item Risk of overfitting with many predictors relative to sample size
#'   \item Useful for confirmatory analyses with pre-specified models
#' }
#' 
#' \emph{Custom Method} (\code{method = "custom"}):
#' \itemize{
#'   \item Manual selection based on subject matter knowledge
#'   \item Runs univariable analysis for all predictors (for comparison)
#'   \item Includes only specified predictors in multivariable model
#'   \item Ideal for theory-driven model building
#'   \item Allows comparison of unadjusted vs adjusted effects for all variables
#' }
#' 
#' \strong{Interpreting Results:}
#' 
#' When \code{columns = "both"} (default), tables show:
#' \itemize{
#'   \item \strong{Univariable columns}: Crude associations, unadjusted for 
#'     other variables. Labeled as "Univariable OR/HR/RR (95\% CI)" and "Uni p"
#'   \item \strong{Multivariable columns}: Adjusted associations, accounting 
#'     for all other predictors in the model. Labeled as "Multivariable aOR/aHR/aRR 
#'     (95\% CI)" and "Multi p" ("a" = adjusted)
#'   \item Variables not meeting selection criteria show "-" in multivariable columns
#' }
#' 
#' Comparing univariable and multivariable results helps identify:
#' \itemize{
#'   \item \strong{Confounding}: Large changes in effect estimates
#'   \item \strong{Independent effects}: Similar univariable and multivariable 
#'     estimates
#'   \item \strong{Mediation}: Attenuated effects in multivariable model
#'   \item \strong{Suppression}: Effects that emerge only after adjustment
#' }
#' 
#' \strong{Sample Size Considerations:}
#' 
#' Rule of thumb for multivariable models:
#' \itemize{
#'   \item \strong{Logistic regression}: ≥10 events per predictor variable
#'   \item \strong{Cox regression}: ≥10 events per predictor variable  
#'   \item \strong{Linear regression}: ≥10-20 observations per predictor
#' }
#' 
#' Use screening methods to reduce predictor count when these ratios are not met.
#'
#' @seealso 
#' \code{\link{uscreen}} for univariable screening only,
#' \code{\link{fit}} for fitting a single multivariable model,
#' \code{\link{compfit}} for comparing multiple models,
#' \code{\link{desctbl}} for descriptive statistics
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' data(clintrial_labels)
#' library(survival)
#' 
#' # Example 1: Basic screening with p < 0.05 threshold
#' result1 <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "hypertension", 
#'                   "diabetes", "treatment", "stage"),
#'     method = "screen",
#'     p_threshold = 0.05,
#'     var_labels = clintrial_labels
#' )
#' print(result1)
#' # Shows both univariable and multivariable results
#' # Only significant univariable predictors in multivariable model
#' 
#' # Example 2: Liberal screening threshold (p < 0.20)
#' result2 <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "hypertension",
#'                   "diabetes", "ecog", "treatment", "stage", "grade"),
#'     method = "screen",
#'     p_threshold = 0.20,  # More liberal for screening
#'     var_labels = clintrial_labels
#' )
#' print(result2)
#' 
#' # Example 3: Include all predictors (no selection)
#' result3 <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     method = "all",
#'     var_labels = clintrial_labels
#' )
#' print(result3)
#' # All predictors in both analyses
#' 
#' # Example 4: Custom variable selection
#' result4 <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "treatment", "stage"),
#'     method = "custom",
#'     multi_predictors = c("age", "treatment", "stage"),  # Manual selection
#'     var_labels = clintrial_labels
#' )
#' print(result4)
#' # Univariable for all, multivariable for selected only
#' 
#' # Example 5: Cox regression with screening
#' cox_result <- fastfit(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "bmi", "treatment", "stage", "grade"),
#'     model_type = "coxph",
#'     method = "screen",
#'     p_threshold = 0.10,
#'     var_labels = clintrial_labels
#' )
#' print(cox_result)
#' # Returns hazard ratios (HR/aHR)
#' 
#' # Example 6: Show only multivariable results
#' multi_only <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     method = "all",
#'     columns = "multi",  # Multivariable results only
#'     var_labels = clintrial_labels
#' )
#' print(multi_only)
#' 
#' # Example 7: Show only univariable results
#' uni_only <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     columns = "uni",  # Univariable results only
#'     var_labels = clintrial_labels
#' )
#' print(uni_only)
#' 
#' # Example 8: Show only effect estimates (no p-values)
#' effects_only <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment"),
#'     metrics = "effect",  # Effect estimates only
#'     var_labels = clintrial_labels
#' )
#' print(effects_only)
#' 
#' # Example 9: Show only p-values (no effect estimates)
#' pvals_only <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment"),
#'     metrics = "p",  # P-values only
#'     var_labels = clintrial_labels
#' )
#' print(pvals_only)
#' 
#' # Example 10: Return both table and model object
#' both <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     method = "all",
#'     return_type = "both"
#' )
#' 
#' # Access the table
#' print(both$table)
#' 
#' # Access the model
#' summary(both$model)
#' 
#' # Model diagnostics
#' plot(both$model)
#' 
#' # Example 11: Return only the model object
#' model_only <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment"),
#'     return_type = "model"
#' )
#' 
#' # This returns a glm object directly
#' summary(model_only)
#' 
#' # Example 12: Linear regression
#' linear_result <- fastfit(
#'     data = clintrial,
#'     outcome = "bmi",
#'     predictors = c("age", "sex", "smoking", "creatinine"),
#'     model_type = "lm",
#'     method = "all",
#'     var_labels = clintrial_labels
#' )
#' print(linear_result)
#' 
#' # Example 13: Poisson regression
#' poisson_result <- fastfit(
#'     data = clintrial,
#'     outcome = "los_days",
#'     predictors = c("age", "treatment", "surgery", "stage"),
#'     model_type = "glm",
#'     family = "poisson",
#'     method = "all",
#'     var_labels = clintrial_labels
#' )
#' print(poisson_result)
#' # Returns rate ratios (RR/aRR)
#' 
#' # Example 14: Keep univariable models for diagnostics
#' with_models <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "bmi", "creatinine"),
#'     keep_models = TRUE
#' )
#' 
#' # Access univariable models
#' uni_results <- attr(with_models, "uni_results")
#' uni_models <- attr(uni_results, "models")
#' summary(uni_models[["age"]])
#' 
#' # Example 15: Check how many predictors made it to multivariable
#' result <- fastfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "bmi", "smoking", "hypertension",
#'                   "diabetes", "ecog", "treatment", "stage", "grade"),
#'     method = "screen",
#'     p_threshold = 0.10
#' )
#' 
#' n_multi <- attr(result, "n_multi")
#' cat("Predictors in multivariable model:", n_multi, "\n")
#' 
#' # Example 16: Complete publication workflow
#' final_table <- fastfit(
#'     data = clintrial,
#'     outcome = "Surv(os_months, os_status)",
#'     predictors = c("age", "sex", "race", "bmi", "smoking", 
#'                   "hypertension", "diabetes", "ecog",
#'                   "treatment", "stage", "grade"),
#'     model_type = "coxph",
#'     method = "screen",
#'     p_threshold = 0.10,
#'     columns = "both",
#'     metrics = "both",
#'     var_labels = clintrial_labels,
#'     digits = 2,
#'     digits_p = 3
#' )
#' print(final_table)
#' 
#' # Can export directly to PDF/LaTeX/HTML for publication
#' # tbl2pdf(final_table, "regression_results.pdf")
#'
#' @export
fastfit <- function(data,
                    outcome, 
                    predictors,
                    method = "screen",
                    multi_predictors = NULL,
                    p_threshold = 0.05,
                    columns = "both",
                    model_type = "glm",
                    family = "binomial",
                    conf_level = 0.95,
                    add_reference_rows = TRUE,
                    show_n = TRUE,
                    show_events = TRUE,
                    digits = 2,
                    digits_p = 3,
                    var_labels = NULL,
                    metrics = "both",
                    return_type = "table",
                    keep_models = FALSE,
                    exponentiate = NULL,
                    ...) {
    
    ## Input validation
    method <- match.arg(method, c("screen", "all", "custom"))
    columns <- match.arg(columns, c("both", "uni", "multi"))
    return_type <- match.arg(return_type, c("table", "model", "both"))

    if (method == "custom" && is.null(multi_predictors)) {
        stop("multi_predictors must be specified when method='custom'")
    }

    ## Convert metrics to standardized format
    if (length(metrics) == 1 && metrics == "both") {
        metrics <- c("effect", "p")
    }

    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }

    ## Step 1: Univariable analysis (if needed)
    uni_results <- NULL
    uni_raw <- NULL
    
    if (columns %in% c("both", "uni")) {
        message("Running univariable analysis...")
        uni_results <- uscreen(
            data = data,
            outcome = outcome,
            predictors = predictors,
            model_type = model_type,
            family = family,
            conf_level = conf_level,
            add_reference_rows = add_reference_rows,
            show_n = show_n,
            show_events = show_events,
            digits = digits,
            digits_p = digits_p,
            var_labels = var_labels,
            keep_models = keep_models,
            exponentiate = exponentiate,
            ...
        )
        ## Extract raw data for variable selection
        uni_raw <- attr(uni_results, "raw_data")
    }

    ## Step 2: Determine predictors for multivariable model
    multi_vars <- NULL
    multi_model <- NULL
    multi_results <- NULL
    multi_raw <- NULL

    if (columns %in% c("both", "multi")) {
        if (method == "screen") {
            ## Screen based on p-value threshold using raw data
            if (is.null(uni_raw)) {
                ## Need to run univariable if not already done
                uni_temp <- uscreen(data, outcome, predictors, model_type, 
                                    family, conf_level = conf_level,
                                    add_reference_rows = add_reference_rows, ...)
                uni_raw <- attr(uni_temp, "raw_data")
            }
            
            ## Use raw data for filtering
            multi_vars <- unique(uni_raw[p_value <= p_threshold]$predictor)
            
            if (length(multi_vars) == 0) {
                warning("No variables meet p <= ", p_threshold, " threshold")
                if (return_type == "model") return(NULL)
            }
            
        } else if (method == "all") {
            ## Use all predictors
            multi_vars <- predictors
            
        } else if (method == "custom") {
            ## Use specified predictors
            multi_vars <- multi_predictors
        }
        
        ## Fit multivariable model if we have predictors
        if (length(multi_vars) > 0) {
            message(sprintf("Fitting multivariable model with %d predictors...", 
                            length(multi_vars)))
            
            multi_results <- fit(
                data = data,
                outcome = outcome,
                predictors = multi_vars,
                model_type = model_type,
                family = family,
                conf_level = conf_level,
                show_n = show_n,
                show_events = show_events,
                digits = digits,
                digits_p = digits_p,
                var_labels = var_labels,
                keep_qc_stats = FALSE,  # Don't need QC stats for display
                add_reference_rows = add_reference_rows,
                exponentiate = exponentiate,
                ...
            )
            
            ## Extract model and raw data
            multi_model <- attr(multi_results, "model")
            multi_raw <- attr(multi_results, "raw_data")
        }
    }

    ## Step 3: Handle return types
    if (return_type == "model") {
        return(multi_model)
    }

    ## Step 4: Format combined output using the new formatted tables
    result <- format_fastfit_combined(
        uni_formatted = uni_results,
        multi_formatted = multi_results,
        uni_raw = uni_raw,
        multi_raw = multi_raw,
        predictors = predictors,
        columns = columns,
        metrics = metrics,
        show_n = show_n,
        show_events = show_events,
        var_labels = var_labels,
        exponentiate = exponentiate
    )

    ## Add attributes
    setattr(result, "outcome", outcome)
    setattr(result, "model_type", model_type)
    setattr(result, "method", method)
    setattr(result, "columns", columns)
    setattr(result, "uni_raw", uni_raw)
    setattr(result, "multi_raw", multi_raw)

    if (!is.null(multi_model)) {
        setattr(result, "model", multi_model)
    }

    if (columns != "uni" && length(multi_vars) > 0) {
        setattr(result, "n_multi", length(multi_vars))
    }

    class(result) <- c("fastfit_result", class(result))

    if (return_type == "both") {
        return(list(table = result, model = multi_model))
    } else {
        return(result)
    }
}

#' Format combined fastfit output from formatted tables
#' @keywords internal
format_fastfit_combined <- function(uni_formatted, multi_formatted, 
                                    uni_raw, multi_raw,
                                    predictors, columns, metrics, 
                                    show_n, show_events, var_labels,
                                    exponentiate = NULL) {
    
    ## Determine effect column name from the formatted tables
    effect_cols <- grep("\\(95% CI\\)$", names(uni_formatted %||% multi_formatted), value = TRUE)
    effect_type <- if (length(effect_cols) > 0) {
                       gsub(" \\(95% CI\\)", "", effect_cols[1])
                   } else {
                       "Effect"
                   }
    
    result <- data.table::data.table()
    
    ## Get unique variables from both tables
    all_vars <- unique(c(
        if (!is.null(uni_formatted)) uni_formatted$Variable[uni_formatted$Variable != ""],
        if (!is.null(multi_formatted)) multi_formatted$Variable[multi_formatted$Variable != ""]
    ))
    
    for (var in all_vars) {
        ## Get rows for this variable from formatted tables
        uni_var_rows <- if (!is.null(uni_formatted)) {
                            ## Find the variable and its subsequent rows
                            var_start <- which(uni_formatted$Variable == var)
                            if (length(var_start) > 0) {
                                var_end <- min(c(which(uni_formatted$Variable != "" & 
                                                       seq_len(nrow(uni_formatted)) > var_start[1]),
                                                 nrow(uni_formatted) + 1)) - 1
                                uni_formatted[var_start[1]:var_end]
                            } else {
                                NULL
                            }
                        } else {
                            NULL
                        }
        
        multi_var_rows <- if (!is.null(multi_formatted)) {
                              var_start <- which(multi_formatted$Variable == var)
                              if (length(var_start) > 0) {
                                  var_end <- min(c(which(multi_formatted$Variable != "" & 
                                                         seq_len(nrow(multi_formatted)) > var_start[1]),
                                                   nrow(multi_formatted) + 1)) - 1
                                  multi_formatted[var_start[1]:var_end]
                              } else {
                                  NULL
                              }
                          } else {
                              NULL
                          }
        
        ## Determine number of rows needed (max of uni and multi)
        n_rows <- max(
            if (!is.null(uni_var_rows)) nrow(uni_var_rows) else 0,
            if (!is.null(multi_var_rows)) nrow(multi_var_rows) else 0
        )
        
        if (n_rows == 0) next
        
        ## Build combined rows
        for (i in seq_len(n_rows)) {
            row <- data.table::data.table()
            
            ## Variable and Group columns from either source
            if (!is.null(uni_var_rows) && i <= nrow(uni_var_rows)) {
                row[, Variable := uni_var_rows$Variable[i]]
                row[, Group := uni_var_rows$Group[i]]
                if (show_n && "n" %in% names(uni_var_rows)) {
                    row[, n := uni_var_rows$n[i]]
                }
                if (show_events && "Events" %in% names(uni_var_rows)) {
                    row[, Events := uni_var_rows$Events[i]]
                }
            } else if (!is.null(multi_var_rows) && i <= nrow(multi_var_rows)) {
                row[, Variable := multi_var_rows$Variable[i]]
                row[, Group := multi_var_rows$Group[i]]
                if (show_n && "n" %in% names(multi_var_rows)) {
                    row[, n := multi_var_rows$n[i]]
                }
                if (show_events && "Events" %in% names(multi_var_rows)) {
                    row[, Events := multi_var_rows$Events[i]]
                }
            }
            
            ## Univariable columns
            if (columns %in% c("both", "uni") && !is.null(uni_var_rows) && i <= nrow(uni_var_rows)) {
                effect_col <- grep("\\(95% CI\\)$", names(uni_var_rows), value = TRUE)[1]
                if ("effect" %in% metrics && !is.na(effect_col)) {
                    row[, uni_effect := uni_var_rows[[effect_col]][i]]
                }
                if ("p" %in% metrics && "p-value" %in% names(uni_var_rows)) {
                    row[, uni_p := uni_var_rows[["p-value"]][i]]
                }
            } else if (columns %in% c("both", "uni")) {
                if ("effect" %in% metrics) row[, uni_effect := ""]
                if ("p" %in% metrics) row[, uni_p := ""]
            }
            
            ## Multivariable columns
            if (columns %in% c("both", "multi") && !is.null(multi_var_rows) && i <= nrow(multi_var_rows)) {
                effect_col <- grep("\\(95% CI\\)$", names(multi_var_rows), value = TRUE)[1]
                if ("effect" %in% metrics && !is.na(effect_col)) {
                    row[, multi_effect := multi_var_rows[[effect_col]][i]]
                }
                if ("p" %in% metrics && "p-value" %in% names(multi_var_rows)) {
                    row[, multi_p := multi_var_rows[["p-value"]][i]]
                }
            } else if (columns %in% c("both", "multi")) {
                ## Variable not in multivariable model
                if ("effect" %in% metrics) row[, multi_effect := "-"]
                if ("p" %in% metrics) row[, multi_p := "-"]
            }
            
            result <- rbind(result, row, fill = TRUE)
        }
    }
    
    ## Clean up column names for display
    if (columns == "both") {
        if ("uni_effect" %in% names(result)) {
            ## Determine the effect type from raw data and exponentiate parameter
            if (!is.null(exponentiate)) {
                if (exponentiate) {
                                        # User explicitly wants exponentiated values
                    if (!is.null(uni_raw)) {
                        if ("OR" %in% names(uni_raw) || "exp_coef" %in% names(uni_raw)) {
                            effect_type <- "OR"
                        } else if ("HR" %in% names(uni_raw) || "exp_coef" %in% names(uni_raw)) {
                            effect_type <- "HR"
                        } else if ("RR" %in% names(uni_raw) || "exp_coef" %in% names(uni_raw)) {
                            effect_type <- "RR"
                        } else {
                            effect_type <- "Exp(Coef)"
                        }
                    }
                } else {
                                        # User explicitly wants coefficients
                    effect_type <- "Coefficient"
                }
            } else {
                                        # Use existing auto-detection
                effect_type <- if (!is.null(uni_raw) && "OR" %in% names(uni_raw)) "OR"
                               else if (!is.null(uni_raw) && "HR" %in% names(uni_raw)) "HR"
                               else if (!is.null(uni_raw) && "RR" %in% names(uni_raw)) "RR"
                               else "Estimate"
            }
            
            setnames(result, "uni_effect", paste0("Univariable ", effect_type, " (95% CI)"))
        }
        
        if ("uni_p" %in% names(result)) {
            setnames(result, "uni_p", "Uni p")
        }
        
        if ("multi_effect" %in% names(result)) {
            ## Determine the adjusted effect type
            if (!is.null(exponentiate)) {
                if (exponentiate) {
                                        # User explicitly wants exponentiated values - use adjusted notation
                    if (!is.null(multi_raw)) {
                        if ("OR" %in% names(multi_raw) || "exp_coef" %in% names(multi_raw)) {
                            effect_type <- "aOR"  # Preserve adjusted notation
                        } else if ("HR" %in% names(multi_raw) || "exp_coef" %in% names(multi_raw)) {
                            effect_type <- "aHR"  # Preserve adjusted notation
                        } else if ("RR" %in% names(multi_raw) || "exp_coef" %in% names(multi_raw)) {
                            effect_type <- "aRR"  # Preserve adjusted notation
                        } else {
                            effect_type <- "Adj. Exp(Coef)"
                        }
                    }
                } else {
                                        # User explicitly wants coefficients
                    effect_type <- "Adj. Coefficient"
                }
            } else {
                                        # Use existing auto-detection with adjusted notation
                effect_type <- if (!is.null(multi_raw) && "OR" %in% names(multi_raw)) "aOR"
                               else if (!is.null(multi_raw) && "HR" %in% names(multi_raw)) "aHR"
                               else if (!is.null(multi_raw) && "RR" %in% names(multi_raw)) "aRR"
                               else "Estimate"
            }
            
            setnames(result, "multi_effect", paste0("Multivariable ", effect_type, " (95% CI)"))
        }
        
        if ("multi_p" %in% names(result)) {
            setnames(result, "multi_p", "Multi p")
        }
        
    } else if (columns == "uni") {
        ## Keep as univariable (same logic as above for univariable)
        if ("uni_effect" %in% names(result)) {
            if (!is.null(exponentiate)) {
                if (exponentiate) {
                    if (!is.null(uni_raw)) {
                        if ("OR" %in% names(uni_raw) || "exp_coef" %in% names(uni_raw)) {
                            effect_type <- "OR"
                        } else if ("HR" %in% names(uni_raw) || "exp_coef" %in% names(uni_raw)) {
                            effect_type <- "HR"
                        } else if ("RR" %in% names(uni_raw) || "exp_coef" %in% names(uni_raw)) {
                            effect_type <- "RR"
                        } else {
                            effect_type <- "Exp(Coef)"
                        }
                    }
                } else {
                    effect_type <- "Coefficient"
                }
            } else {
                effect_type <- if (!is.null(uni_raw) && "OR" %in% names(uni_raw)) "OR"
                               else if (!is.null(uni_raw) && "HR" %in% names(uni_raw)) "HR"
                               else if (!is.null(uni_raw) && "RR" %in% names(uni_raw)) "RR"
                               else "Estimate"
            }
            setnames(result, "uni_effect", paste0("Univariable ", effect_type, " (95% CI)"))
        }
        
    } else if (columns == "multi") {
        ## Use adjusted notation (same logic as multivariable above)
        if ("multi_effect" %in% names(result)) {
            if (!is.null(exponentiate)) {
                if (exponentiate) {
                    if (!is.null(multi_raw)) {
                        if ("OR" %in% names(multi_raw) || "exp_coef" %in% names(multi_raw)) {
                            effect_type <- "aOR"
                        } else if ("HR" %in% names(multi_raw) || "exp_coef" %in% names(multi_raw)) {
                            effect_type <- "aHR"
                        } else if ("RR" %in% names(multi_raw) || "exp_coef" %in% names(multi_raw)) {
                            effect_type <- "aRR"
                        } else {
                            effect_type <- "Adj. Exp(Coef)"
                        }
                    }
                } else {
                    effect_type <- "Adj. Coefficient"
                }
            } else {
                effect_type <- if (!is.null(multi_raw) && "OR" %in% names(multi_raw)) "aOR"
                               else if (!is.null(multi_raw) && "HR" %in% names(multi_raw)) "aHR"
                               else if (!is.null(multi_raw) && "RR" %in% names(multi_raw)) "aRR"
                               else "Estimate"
            }
            setnames(result, "multi_effect", paste0("Multivariable ", effect_type, " (95% CI)"))
        }
    }
    
    return(result)
}

#' Print method for fastfit results
#' @keywords internal
#' @export
print.fastfit_result <- function(x, ...) {
    cat("\nFastfit Analysis Results\n")
    cat("Outcome: ", attr(x, "outcome"), "\n", sep = "")
    cat("Model Type: ", attr(x, "model_type"), "\n", sep = "")
    cat("Method: ", attr(x, "method"), "\n", sep = "")
    
    if (!is.null(attr(x, "n_multi"))) {
        cat("Multivariable predictors: ", attr(x, "n_multi"), "\n", sep = "")
    }
    
    cat("\n")
    NextMethod("print", x)
    invisible(x)
}
