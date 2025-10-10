#' Format model results for publication-ready display
#' @keywords internal
format_model_table <- function(data, 
                               effect_col = NULL,
                               digits = 2, 
                               digits_p = 3, 
                               var_labels = NULL,
                               show_n_events = c("n", "Events"),
                               reference_label = "reference",
                               exponentiate = NULL) {
    
    result <- data.table::copy(data)
    
    ## Standardize column names
    if ("variable" %in% names(result)) {
        setnames(result, "variable", "Variable")
    }
    if ("group" %in% names(result)) {
        setnames(result, "group", "Group")
    }

    ## Handle the exponentiate parameter to choose which columns to use
    if (!is.null(exponentiate)) {
        if (exponentiate && "exp_coef" %in% names(result)) {
            ## Check model type
            if ("OR" %in% names(result)) {
                effect_col <- "OR"
            } else if ("HR" %in% names(result)) {
                effect_col <- "HR"
            } else if ("RR" %in% names(result)) {
                effect_col <- "RR"
            } else {
                ## Generic model - create OR/RR based on model type
                model_type <- unique(result$model_type)[1]
                if (grepl("Logistic", model_type)) {
                    result[, `:=`(
                        OR = exp_coef,
                        CI_lower = exp_lower,
                        CI_upper = exp_upper
                    )]
                    effect_col <- "OR"
                } else if (grepl("Poisson", model_type)) {
                    result[, `:=`(
                        RR = exp_coef,
                        CI_lower = exp_lower,
                        CI_upper = exp_upper
                    )]
                    effect_col <- "RR"
                } else {
                    result[, `:=`(
                        Estimate = exp_coef,
                        CI_lower = exp_lower,
                        CI_upper = exp_upper
                    )]
                    effect_col <- "Estimate"
                }
            }
        } else if (!exponentiate && "coef" %in% names(result)) {
            result[, `:=`(
                Estimate = coef,
                CI_lower = coef_lower,
                CI_upper = coef_upper
            )]
            effect_col <- "Estimate"
        }
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

    ## Format sample size BEFORE eliminating repeated variable names
    ## This preserves the group-specific counts
    if ("n" %in% names(result)) {
        if ("n_group" %in% names(result)) {
            ## Use group-specific n where available
            result[, n := ifelse(!is.na(n_group), 
                                 as.character(n_group),
                                 as.character(n))]
        }
        ## Format with commas
        result[, n := ifelse(!is.na(n) & as.numeric(n) >= 1000,
                             format(as.numeric(n), big.mark = ","),
                             as.character(n))]
    }

    ## Similar for events
    if ("events" %in% names(result)) {
        if ("events_group" %in% names(result)) {
            result[, events := ifelse(!is.na(events_group),
                                      as.character(events_group),
                                      as.character(events))]
        }
        result[, events := ifelse(!is.na(events) & as.numeric(events) >= 1000,
                                  format(as.numeric(events), big.mark = ","),
                                  as.character(events))]
    }
    
    ## NOW eliminate repeated variable names (AFTER processing n and events)
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
            result[!is.na(reference) & reference == reference_label, `p-value` := "-"]
        }
    }

    ## Handle show_n_events parameter
    if (is.null(show_n_events)) {
        ## User explicitly wants no n or events columns
        show_n_events <- character(0)
    } else if (is.character(show_n_events)) {
        ## Validate the input
        valid_options <- c("n", "Events", "events")
        invalid <- setdiff(show_n_events, valid_options)
        if (length(invalid) > 0) {
            warning("Invalid show_n_events options ignored: ", paste(invalid, collapse = ", "))
            show_n_events <- intersect(show_n_events, valid_options)
        }
        ## Normalize "events" to "Events"
        show_n_events[show_n_events == "events"] <- "Events"
    }
    
    ## Select columns for final output
    display_cols <- character()
    
    if ("Variable" %in% names(result)) display_cols <- c(display_cols, "Variable")
    if ("Group" %in% names(result)) display_cols <- c(display_cols, "Group")

    if (("n" %in% show_n_events) && ("n" %in% names(result))) {
        display_cols <- c(display_cols, "n")
    }
    
    if (("Events" %in% show_n_events || "events" %in% show_n_events) && ("events" %in% names(result))) {
        display_cols <- c(display_cols, "events")
    }
    
    if (effect_label %in% names(result)) display_cols <- c(display_cols, effect_label)
    if ("p-value" %in% names(result)) display_cols <- c(display_cols, "p-value")
    
    formatted <- result[, ..display_cols]

    if ("events" %in% names(formatted)) {
        setnames(formatted, "events", "Events")
    }

    return(formatted)
}

#' Format p-values for display
#' @keywords internal
format_pvalues_fit <- function(p, digits = 3) {
    ifelse(is.na(p), "-",
    ifelse(p < 0.001, "< 0.001",
           sprintf(paste0("%.", digits, "f"), p)))
}
