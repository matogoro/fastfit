#' Format model results for publication-ready display
#' @keywords internal
format_model_table <- function(data, 
                               effect_col = NULL,
                               digits = 2, 
                               digits_p = 3, 
                               var_labels = NULL,
                               show_n = TRUE,
                               reference_label = "reference") {
    
    result <- data.table::copy(data)
    
    ## Standardize column names
    if ("variable" %in% names(result)) {
        setnames(result, "variable", "Variable")
    }
    if ("group" %in% names(result)) {
        setnames(result, "group", "Group")
    }
    
    ## Auto-detect effect column if not specified
    if (is.null(effect_col)) {
        effect_col <- intersect(c("OR", "HR", "RR", "Estimate"), names(result))[1]
        if (length(effect_col) == 0) {
            stop("No effect measure column found (OR, HR, RR, or Estimate)")
        }
    }
    
    ## Apply variable labels if provided
    if (!is.null(var_labels) && "Variable" %in% names(result)) {
        for (var_name in names(var_labels)) {
            result[Variable == var_name, Variable := var_labels[var_name]]
        }
    }
    
    ## Clean up Group display (handle empty groups for continuous vars)
    if ("Group" %in% names(result)) {
        result[Group == "", Group := "-"]
    }

                                        # Format sample size columns with group-specific counts
    if (show_n && "n" %in% names(result)) {
                                        # Check if we have group-specific counts
        if ("n_group" %in% names(result)) {
                                        # For rows with n_group, use that; otherwise use n
            result[, display_n := ifelse(!is.na(n_group), 
                                         as.character(n_group),
                                         as.character(n))]
        } else {
                                        # No group counts available, use total n
            result[, display_n := as.character(n)]
        }
        
                                        # Format with commas
        result[, display_n := ifelse(!is.na(display_n) & as.numeric(display_n) >= 1000,
                                     format(as.numeric(display_n), big.mark = ","),
                                     display_n)]
        
                                        # Replace the n column with display_n
        result[, n := display_n]
        result[, display_n := NULL]
    }

                                        # Similar for events
    if ("events" %in% names(result)) {
        if ("events_group" %in% names(result)) {
            result[, display_events := ifelse(!is.na(events_group),
                                              as.character(events_group),
                                              as.character(events))]
        } else {
            result[, display_events := as.character(events)]
        }
        
        result[, display_events := ifelse(!is.na(display_events) & as.numeric(display_events) >= 1000,
                                          format(as.numeric(display_events), big.mark = ","),
                                          display_events)]
        
        result[, events := display_events]
        result[, display_events := NULL]
    }
    
    ## Eliminate repeated variable names (only show on first row for each variable)
    if ("Variable" %in% names(result)) {
        current_var <- ""
        for (i in seq_len(nrow(result))) {
            if (result$Variable[i] == current_var) {
                result[i, Variable := ""]
            } else {
                current_var <- result$Variable[i]
            }
        }
    }
    
    ## Format events column if present
    if ("events" %in% names(result)) {
        ## Similar logic for events_group if available
        if ("events_group" %in% names(result)) {
            result[, display_events := ifelse(!is.na(events_group),
                                              as.character(events_group),
                                              as.character(events))]
        } else {
            result[, display_events := as.character(events)]
        }
        
        result[, display_events := ifelse(!is.na(display_events) & as.numeric(display_events) >= 1000,
                                          format(as.numeric(display_events), big.mark = ","),
                                          display_events)]
    }
    
    ## Create effect column label based on model scope
    model_scope <- if ("model_scope" %in% names(result)) {
                       unique(result$model_scope)[1]
                   } else {
                       "Effect"
                   }
    
    ## Create appropriate label
    if (model_scope == "Univariable") {
        effect_label <- paste0("Univariable ", effect_col, " (95% CI)")
    } else if (model_scope == "Multivariable") {
        adjusted_col <- if (effect_col == "OR") "aOR" 
                        else if (effect_col == "HR") "aHR"
                        else if (effect_col == "RR") "aRR"
                        else effect_col
        effect_label <- paste0("Multivariable ", adjusted_col, " (95% CI)")
    } else {
        effect_label <- paste0(effect_col, " (95% CI)")
    }
    
    ## Format effect sizes with CI
    if ("CI_lower" %in% names(result) && "CI_upper" %in% names(result)) {
        is_reference <- FALSE
        if ("reference" %in% names(result)) {
            is_reference <- !is.na(result$reference) & result$reference == reference_label
        }
        
        if (effect_col %in% c("OR", "HR", "RR")) {
            result[, (effect_label) := ifelse(
                         is_reference,
                         "-",
                                       ifelse(!is.na(get(effect_col)),
                                              sprintf("%.*f (%.*f-%.*f)",
                                                      digits, get(effect_col),
                                                      digits, CI_lower,
                                                      digits, CI_upper),
                                              "")
                     )]
        } else {
            result[, (effect_label) := ifelse(
                         is_reference,
                         "-",
                                       ifelse(!is.na(get(effect_col)),
                                              sprintf("%.*f (%.*f, %.*f)",
                                                      digits, get(effect_col),
                                                      digits, CI_lower,
                                                      digits, CI_upper),
                                              "")
                     )]
        }
    }
    
    ## Format p-values
    if ("p_value" %in% names(result)) {
        result[, `p-value` := format_pvalues_fit(p_value, digits_p)]
        
        if ("reference" %in% names(result)) {
            result[!is.na(reference) & reference == reference_label, `p-value` := ""]
        }
    }
    
    ## Select columns for final output
    display_cols <- character()
    
    if ("Variable" %in% names(result)) display_cols <- c(display_cols, "Variable")
    if ("Group" %in% names(result)) display_cols <- c(display_cols, "Group")

    if (show_n && "n" %in% names(result)) {
        display_cols <- c(display_cols, "n")
    }
    
    if ("events" %in% names(result)) {
        display_cols <- c(display_cols, "events")
    }
    
    if (effect_label %in% names(result)) display_cols <- c(display_cols, effect_label)
    if ("p-value" %in% names(result)) display_cols <- c(display_cols, "p-value")
    
    formatted <- result[, ..display_cols]
    
    return(formatted)
}

#' Format p-values for display
#' @keywords internal
format_pvalues_fit <- function(p, digits = 3) {
    ifelse(is.na(p), "",
    ifelse(p < 0.001, "< 0.001",
           sprintf(paste0("%.", digits, "f"), p)))
}
