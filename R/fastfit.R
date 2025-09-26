#' Perform Complete Statistical Analysis with Flexible Options
#'
#' @param data Data.frame or data.table containing the dataset.
#' @param outcome Character string specifying the outcome variable.
#' @param predictors Character vector of predictor variable names for univariable screening.
#' @param method Character string specifying the analysis method:
#'   "screen" - Univariable to multivariable based on p-value threshold (default).
#'   "all" - All predictors in both univariable and multivariable.
#'   "custom" - Univariable all, multivariable with selected predictors only.
#' @param multi_predictors For method="custom", character vector of predictors for multivariable model.
#' @param p_threshold For method="screen", p-value threshold for selection. Default 0.05.
#' @param columns Character string specifying which results to include:
#'   "both" - Both univariable and multivariable columns (default).
#'   "uni" - Univariable results only.
#'   "multi" - Multivariable results only.
#' @param model_type Character string: "glm", "lm", "coxph", "clogit". Default "glm".
#' @param family For GLM models, the family. Default "binomial".
#' @param conf_level Confidence level. Default 0.95.
#' @param add_reference_rows Add reference category rows. Default TRUE.
#' @param var_labels Named character vector for custom variable labels. Names should
#'   match variable names in predictors, values are display labels.
#' @param metrics Character vector of metrics to include per column:
#'   "effect" - Effect size with CI (OR/HR/RR/Estimate).
#'   "p" - P-value.
#'   "both" - Both effect and p-value (default).
#' @param return_type What to return:
#'   "table" - Formatted data.table (default).
#'   "model" - Multivariable model object only.
#'   "both" - List with table and model.
#' @param keep_models Logical. Store model objects in results. Default FALSE.
#' @param ... Additional arguments passed to model fitting functions.
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
                    var_labels = NULL,  # Add this parameter
                    metrics = "both",
                    return_type = "table",
                    keep_models = FALSE,
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
    if (columns %in% c("both", "uni")) {
        message("Running univariable analysis...")
        uni_results <- uscreen(
            data = data,
            outcome = outcome,
            predictors = predictors,
            model_type = model_type,
            family = family,
            conf_level = conf_level,
            keep_models = keep_models,
            keep_qc_stats = TRUE,
            add_reference_rows = add_reference_rows,
            var_labels = var_labels,
            ...
        )
    }

    ## Step 2: Determine predictors for multivariable model
    multi_vars <- NULL
    multi_model <- NULL
    multi_results <- NULL

    if (columns %in% c("both", "multi")) {
        if (method == "screen") {
            ## Screen based on p-value threshold
            if (is.null(uni_results)) {
                ## Need to run univariable if not already done
                uni_temp <- uscreen(data, outcome, predictors, model_type, 
                                    family, conf_level, FALSE, TRUE, 
                                    add_reference_rows, ...)
                multi_vars <- uni_temp[sig_binary == TRUE | p_value <= p_threshold, 
                                       unique(variable)]
            } else {
                multi_vars <- uni_results[sig_binary == TRUE | p_value <= p_threshold, 
                                          unique(variable)]
            }
            
            if (length(multi_vars) == 0) {
                warning("No significant variables at p <= ", p_threshold)
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
            
            multi_model <- mmodel(
                data = data,
                outcome = outcome,
                predictors = multi_vars,
                model_type = model_type,
                family = family,
                ...
            )
            
            multi_results <- m2dt(
                multi_model,
                conf_level = conf_level,
                variable_name = "Multivariable",
                keep_qc_stats = TRUE,
                add_reference_rows = add_reference_rows
            )

            if (!is.null(var_labels)) {
                                        # Extract base variable name from term
                multi_results[, base_var := gsub("^([^0-9]+).*", "\\1", term)]
                
                                        # Add label column
                multi_results[, label := ifelse(base_var %in% names(var_labels),
                                                var_labels[base_var],
                                                base_var)]
            }
        }
    }

    ## Step 3: Handle return types
    if (return_type == "model") {
        return(multi_model)
    }

    ## Step 4: Format combined output
    result <- format_fastfit_table(
        uni_results = uni_results,
        multi_results = multi_results,
        predictors = predictors,
        columns = columns,
        metrics = metrics,
        model_type = model_type
    )

    ## Add attributes
    data.table::setattr(result, "outcome", outcome)
    data.table::setattr(result, "model_type", model_type)
    data.table::setattr(result, "method", method)
    data.table::setattr(result, "columns", columns)

    if (!is.null(multi_model)) {
        data.table::setattr(result, "model", multi_model)
    }

    if (columns != "uni") {
        data.table::setattr(result, "n_multi", length(multi_vars))
    }

    result[]  # Force finalization

    if (return_type == "both") {
        return(list(table = result, model = multi_model))
    } else {
        return(result)
    }
}

#' Format fastfit table output
#' 
#' @keywords internal
format_fastfit_table <- function(uni_results, multi_results, predictors,
                                 columns, metrics, model_type) {
    
                                        # Determine effect column name
    effect_col <- if (model_type == "coxph") "HR" 
                  else if (model_type == "glm") "OR" 
                  else "Estimate"
    
                                        # Build result table
    result <- data.table::data.table()
    
    for (pred in predictors) {
                                        # Get rows for this predictor from uni_results
        if (!is.null(uni_results)) {
            pred_rows <- uni_results[variable == pred]
            
                                        # Use label if available, otherwise use variable name
            display_name <- if ("label" %in% names(pred_rows) && nrow(pred_rows) > 0) {
                                pred_rows$label[1]
                            } else {
                                pred
                            }
        } else {
            display_name <- pred
            pred_rows <- data.table::data.table(
                                         term = pred,
                                         reference = ""
                                     )
        }
        
        for (i in seq_len(nrow(pred_rows))) {
            row <- data.table::data.table(
                                   variable = if (i == 1) display_name else "",  # Use display_name here
                                   level = gsub(paste0("^", pred), "", pred_rows$term[i])
                               )
            
                                        # Clean up level display
            if (row$level == "") row[, level := "(main effect)"]
            
                                        # Add sample size if available
            if (!is.null(uni_results) && "n" %in% names(uni_results)) {
                row[, n := pred_rows$n[i]]
                if ("events" %in% names(pred_rows)) {
                    row[, events := pred_rows$events[i]]
                }
            }
            
            ## Univariable columns
            if (columns %in% c("both", "uni") && !is.null(uni_results)) {
                if ("effect" %in% metrics) {
                    ## FIX: Check for NA and handle properly
                    if ("reference" %in% names(pred_rows) && 
                        !is.na(pred_rows$reference[i]) && 
                        pred_rows$reference[i] != "") {
                        row[, uni_effect := pred_rows$reference[i]]
                    } else if (!is.na(pred_rows[[effect_col]][i])) {
                        row[, uni_effect := sprintf("%.2f (%.2f-%.2f)",
                                                    pred_rows[[effect_col]][i],
                                                    pred_rows$CI_lower[i],
                                                    pred_rows$CI_upper[i])]
                    } else {
                        row[, uni_effect := ""]
                    }
                }
                
                if ("p" %in% metrics) {
                    if (!is.na(pred_rows$p_value[i])) {
                        row[, uni_p := format_p_value(pred_rows$p_value[i])]
                    } else {
                        row[, uni_p := ""]
                    }
                }
            }
            
            ## Multivariable columns
            if (columns %in% c("both", "multi") && !is.null(multi_results)) {
                multi_row <- multi_results[term == pred_rows$term[i]]
                
                if (nrow(multi_row) > 0) {
                    if ("effect" %in% metrics) {
                        ## FIX: Check for NA and handle properly
                        if ("reference" %in% names(multi_row) && 
                            !is.na(multi_row$reference[1]) && 
                            multi_row$reference[1] != "") {
                            row[, multi_effect := multi_row$reference[1]]
                        } else if (!is.na(multi_row[[effect_col]][1])) {
                            row[, multi_effect := sprintf("%.2f (%.2f-%.2f)",
                                                          multi_row[[effect_col]][1],
                                                          multi_row$CI_lower[1],
                                                          multi_row$CI_upper[1])]
                        } else {
                            row[, multi_effect := ""]
                        }
                    }
                    
                    if ("p" %in% metrics) {
                        if (!is.na(multi_row$p_value[1])) {
                            row[, multi_p := format_p_value(multi_row$p_value[1])]
                        } else {
                            row[, multi_p := ""]
                        }
                    }
                } else {
                    ## Not in multivariable model
                    if ("effect" %in% metrics) row[, multi_effect := "-"]
                    if ("p" %in% metrics) row[, multi_p := "-"]
                }
            }
            
            result <- rbind(result, row, fill = TRUE)
        }
    }
    
    ## Clean up column names for display
    if (columns == "both") {
        if ("uni_effect" %in% names(result)) {
            data.table::setnames(result, "uni_effect", 
                                 paste0("Univariable ", effect_col, " (95% CI)"))
        }
        if ("uni_p" %in% names(result)) {
            data.table::setnames(result, "uni_p", "Uni p")
        }
        if ("multi_effect" %in% names(result)) {
            data.table::setnames(result, "multi_effect", 
                                 paste0("Multivariable ", effect_col, " (95% CI)"))
        }
        if ("multi_p" %in% names(result)) {
            data.table::setnames(result, "multi_p", "Multi p")
        }
    } else if (columns == "uni") {
        if ("uni_effect" %in% names(result)) {
            data.table::setnames(result, "uni_effect", paste0(effect_col, " (95% CI)"))
        }
        if ("uni_p" %in% names(result)) {
            data.table::setnames(result, "uni_p", "p-value")
        }
    } else if (columns == "multi") {
        if ("multi_effect" %in% names(result)) {
            data.table::setnames(result, "multi_effect", paste0(effect_col, " (95% CI)"))
        }
        if ("multi_p" %in% names(result)) {
            data.table::setnames(result, "multi_p", "p-value")
        }
    }
    
    return(result)
}

#' Format p-values for display
#' @keywords internal
format_p_value <- function(p) {
    if (p < 0.001) "< 0.001"
    else sprintf("%.3f", p)
}
