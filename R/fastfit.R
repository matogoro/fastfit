#' Perform Complete Statistical Analysis Pipeline
#'
#' Executes a comprehensive regression analysis workflow including univariable
#' screening, automatic or manual variable selection, and multivariable modeling.
#' Supports multiple model types with publication-ready output formatting.
#'
#' @param data Data.frame or data.table containing the analysis dataset.
#' @param outcome Character string specifying the outcome variable name.
#' @param predictors Character vector of predictor variable names to analyze.
#' @param method Character string specifying variable selection strategy:
#'   "screen" (univariable to multivariable based on p-value threshold),
#'   "all" (all predictors in both analyses), or
#'   "custom" (all in univariable, selected in multivariable). Default is "screen".
#' @param multi_predictors Character vector of predictors for multivariable
#'   model when method is "custom". Default is NULL.
#' @param p_threshold Numeric p-value threshold for variable selection when
#'   method is "screen". Default is 0.05.
#' @param columns Character string specifying which columns to display:
#'   "both", "uni", or "multi". Default is "both".
#' @param model_type Character string: "glm", "lm", "coxph", or "clogit".
#'   Default is "glm".
#' @param family For GLM models, the error distribution and link function.
#'   Default is "binomial".
#' @param conf_level Numeric confidence level for intervals. Default is 0.95.
#' @param add_reference_rows Logical. If TRUE, includes reference category rows.
#'   Default is TRUE.
#' @param show_n_events Character vector specifying which optional columns to display.
#'   Options: "n", "events" (or "Events"). Default is c("n", "events") for 
#'   survival/logistic models, "n" only for other models. Set to NULL to hide
#'   these columns entirely.
#' @param digits Integer specifying decimal places for effect estimates.
#'   Default is 2.
#' @param digits_p Integer specifying decimal places for p-values.
#'   Default is 3.
#' @param var_labels Named character vector for display labels. Names match
#'   predictors, values are labels. Default is NULL.
#' @param metrics Character specification for statistics: "effect" (OR/HR/estimate
#'   with CI), "p" (p-value only), or "both" (default).
#' @param return_type Character string: "table" (formatted results), "model"
#'   (multivariable model only), or "both" (list with table and model).
#'   Default is "table".
#' @param keep_models Logical. If TRUE, stores univariable model objects.
#'   Default is FALSE.
#' @param exponentiate Logical. Whether to exponentiate coefficients. Default is NULL,
#'   which automatically displays exponentiated coefficients for logistic/Poisson/Cox
#'   regression models and raw coefficients for log-link or linear regression models.
#' @param ... Additional arguments passed to model fitting functions.
#'
#' @return A data.table with formatted results (when return_type is "table"),
#'   a model object (when "model"), or a list with both (when "both").
#'   The table includes variable, group, sample sizes, and effect/p-value
#'   columns as requested. Includes attributes for outcome, model_type,
#'   method, columns, model (if fitted), and n_multi.
#'
#' @details
#' The function implements a complete regression workflow:
#' 
#' 1. Univariable screening: Each predictor is tested individually against
#'    the outcome using the specified model type.
#'    
#' 2. Variable selection: Based on the method parameter, variables are selected
#'    for multivariable analysis.
#'    
#' 3. Multivariable modeling: Selected variables are combined in a single model.
#' 
#' 4. Output formatting: Results are formatted for publication with appropriate
#'    effect measures (OR for logistic, HR for Cox, estimates for linear).
#'
#' The automatic p-value threshold screening (method = "screen") helps identify
#' potentially important predictors while reducing multicollinearity. The custom
#' method allows for theory-driven model building.
#'
#' @examples
#' \dontrun{
#' # Basic logistic regression with screening
#' result <- fastfit(data = mydata,
#'                   outcome = "disease",
#'                   predictors = c("age", "sex", "bmi", "smoking"),
#'                   method = "screen",
#'                   p_threshold = 0.2)
#' 
#' # Cox regression with all variables
#' cox_result <- fastfit(data = survival_data,
#'                       outcome = "survival",
#'                       predictors = c("age", "stage", "grade"),
#'                       model_type = "coxph",
#'                       method = "all")
#' 
#' # Custom selection with labels
#' labels <- c(age = "Age (years)", sex = "Sex", bmi = "Body Mass Index")
#' custom_result <- fastfit(data = mydata,
#'                          outcome = "outcome",
#'                          predictors = c("age", "sex", "bmi", "smoking"),
#'                          method = "custom",
#'                          multi_predictors = c("age", "sex"),
#'                          var_labels = labels)
#' 
#' # Return both table and model
#' both <- fastfit(data = mydata,
#'                 outcome = "disease",
#'                 predictors = vars,
#'                 return_type = "both")
#' both$table  # Access table
#' summary(both$model)  # Access model
#' 
#' # Export results
#' exporttbl(result, "regression_results.pdf")
#' }
#'
#' @seealso 
#' \code{\link{uscreen}} for univariable screening only,
#' \code{\link{mmodel}} for multivariable modeling only,
#' \code{\link{m2dt}} for converting models to tables,
#' \code{\link{exporttbl}} for exporting results to PDF/LaTeX/HTML,
#' \code{\link{desctbl}} for descriptive statistics
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
                    show_n_events = c("n", "events"),
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
            show_n_events = show_n_events,
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
                show_n_events = show_n_events,
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
        show_n_events = show_n_events,
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
                                    show_n_events, var_labels,
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
                if ("n" %in% show_n_events && "n" %in% names(uni_var_rows)) {
                    row[, n := uni_var_rows$n[i]]
                }
                if (("Events" %in% show_n_events || "events" %in% show_n_events) && "Events" %in% names(uni_var_rows)) {
                    row[, Events := uni_var_rows$Events[i]]
                }
            } else if (!is.null(multi_var_rows) && i <= nrow(multi_var_rows)) {
                row[, Variable := multi_var_rows$Variable[i]]
                row[, Group := multi_var_rows$Group[i]]
                if ("n" %in% show_n_events && "n" %in% names(multi_var_rows)) {
                    row[, n := multi_var_rows$n[i]]
                }
                if (("Events" %in% show_n_events || "events" %in% show_n_events) && "Events" %in% names(multi_var_rows)) {
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
