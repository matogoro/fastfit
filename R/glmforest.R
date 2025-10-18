#' Create a Forest Plot for Generalized Linear Models
#'
#' Generates a publication-ready forest plot combining a data table and 
#' graphical representation of odds ratios (or other effect measures) from a GLM.
#' The plot includes variable names, group levels, sample sizes, effect estimates 
#' with confidence intervals, and p-values.
#'
#' @param model A generalized linear model object of class \code{glm}.
#' @param data A data frame or data.table containing the original data used to 
#'   fit the model. If NULL (default), the function attempts to extract data 
#'   from the model object.
#' @param title Character string specifying the plot title. Default is "Generalized Linear Model".
#' @param effect_label Character string for the effect measure label. Default is "Odds Ratio" 
#'   for logistic regression, "Risk Ratio" for log-link models.
#' @param digits Integer specifying the number of decimal places for effect estimates 
#'   and confidence intervals. Default is 2.
#' @param font_size Numeric value controlling the base font size. 
#'   Default is 1.0.
#' @param annot_size Numeric value controlling the relative font size for annotations. 
#'   Default is 3.88.
#' @param header_size Numeric value controlling the relative font size for headers. 
#'   Default is 5.82.
#' @param title_size Numeric value controlling the relative font size for the plot title. 
#'   Default is 23.28.
#' @param tbl_width Numeric value between 0 and 1 specifying the proportion of plot 
#'   width allocated to the data table. If NULL (default), automatically calculated
#'   based on content and plot dimensions.
#' @param plot_width Numeric value specifying the intended plot width in inches.
#'   Used to optimize tbl_width calculation. Default is NULL.
#' @param plot_height Numeric value specifying the intended plot height in inches.
#'   Default is NULL.
#' @param col_width_var Numeric value specifying the proportion of plot width 
#'   allocated to the variable column. Default is 0.01.
#' @param col_width_level Numeric value specifying the proportion of plot width 
#'   allocated to the group column. Default is 0.17.
#' @param col_width_n Numeric value specifying the proportion of plot width 
#'   allocated to the sample size column. Default is 0.17.
#' @param col_width_events Numeric value specifying the proportion of plot width 
#'   allocated to the sample size column. Default is 0.06.
#' @param col_width_or Numeric value specifying the proportion of plot width 
#'   allocated to the effect estimate column. Default is 0.05.
#' @param show_n_events Character vector specifying which columns to show. 
#'   Default is c("n", "Events"). Use NULL for neither.
#' @param indent_groups Logical. Whether to indent group levels under variables.
#'   Default is FALSE.
#' @param condense_table Logical. Whether to condense binary variables into single rows.
#'   Forces indent_groups = TRUE. Default is FALSE.
#' @param ref_label Character string to display for reference categories. 
#'   Default is "reference".
#' @param var_labels Named character vector for custom variable labels. Names should 
#'   match variable names in the model, values are the display labels.
#' @param color Character string specifying the color for effect estimate point estimates. 
#'   Default is "#3A7E8A" (teal).
#' @param exponentiate Logical. Whether to exponentiate coefficients. Default is TRUE 
#'   for logit and log links, FALSE otherwise.
#'
#' @return A ggplot object containing the forest plot. The object has an attribute
#'   "recommended_dims" with suggested width and height in inches.
#'
#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom stats AIC confint deviance nobs
#' @importFrom grDevices axisTicks
#'
#' @seealso 
#' \code{\link[stats]{glm}} for fitting GLM models,
#' \code{\link[ggplot2]{ggplot}} for the underlying plotting system
#'
glmforest <- function(model, data = NULL,
                      title = "Generalized Linear Model",
                      effect_label = NULL,
                      digits = 2,
                      font_size = 1.0,
                      annot_size = 3.88,
                      header_size = 5.82,
                      title_size = 23.28,
                      plot_width = NULL,
                      plot_height = NULL,
                      tbl_width = 0.6,
                      col_width_var = 0.01,
                      col_width_level = 0.17,
                      col_width_n = 0.17,
                      col_width_events = 0.06,
                      col_width_or = 0.05,
                      show_n_events = c("n", "Events"),
                      indent_groups = FALSE,
                      condense_table = FALSE,
                      ref_label = "reference",
                      var_labels = NULL,
                      color = "#3A7E8A",
                      exponentiate = NULL) {
    
    ## Check for required packages
    if (!requireNamespace("data.table", quietly = TRUE)) {
        stop("Package 'data.table' is required but not installed.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("Package 'ggplot2' is required but not installed.")
    }
    if (!requireNamespace("grid", quietly = TRUE)) {
        stop("Package 'grid' is required but not installed.")
    }
    
    stopifnot(inherits(model, "glm"))
    
    ## Determine if we should exponentiate based on link function
    if(is.null(exponentiate)) {
        link_func <- model$family$link
        exponentiate <- link_func %in% c("logit", "log", "cloglog")
    }
    
    ## Set effect label based on model family and link
    if(is.null(effect_label)) {
        if(model$family$family == "binomial" && model$family$link == "logit") {
            effect_label <- "Odds Ratio"
        } else if(model$family$link == "log") {
            effect_label <- "Risk Ratio"
        } else if(exponentiate) {
            effect_label <- "Exp(Coefficient)"
        } else {
            effect_label <- "Coefficient"
        }
    }
    
    ## Get data
    if(is.null(data)){
        warning("The `data` argument is not provided. Data will be extracted from model fit.")
        data <- model$data
        if (is.null(data))
            data <- model$model
        if (is.null(data))
            stop("The `data` argument should be provided either to glmforest or glm.")
    }
    
    ## Convert to data.table if not already
    if(!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
    terms <- attr(model$terms, "dataClasses")[-1]
    
    ## Extract coefficients
    coef_summary <- summary(model)$coefficients
    conf_int <- stats::confint(model)
    
    ## Remove intercept if present
    if("(Intercept)" %in% rownames(coef_summary)) {
        intercept_idx <- which(rownames(coef_summary) == "(Intercept)")
        coef_summary <- coef_summary[-intercept_idx, , drop = FALSE]
        conf_int <- conf_int[-intercept_idx, , drop = FALSE]
    }
    
    coef <- data.table::data.table(
                            term = rownames(coef_summary),
                            estimate = coef_summary[, "Estimate"],
                            std_error = coef_summary[, "Std. Error"],
                            statistic = coef_summary[, "z value"],
                            p_value = coef_summary[, "Pr(>|z|)"],
                            conf_low = conf_int[, 1],
                            conf_high = conf_int[, 2]
                        )
    
    ## Get model statistics
    gmodel <- list(
        nobs = nobs(model),
        null_deviance = model$null.deviance,
        residual_deviance = model$deviance,
        AIC = stats::AIC(model),
        family = paste0(model$family$family, " (", model$family$link, " link)")
    )
    
    ## Format observations and AIC with commas
    gmodel$nobs_formatted <- format(gmodel$nobs, big.mark = ",", scientific = FALSE)
    gmodel$AIC_formatted <- format(round(gmodel$AIC, 2), big.mark = ",", scientific = FALSE, nsmall = 2)
    
    ## Calculate pseudo R-squared (McFadden's)
    gmodel$pseudo_r2 <- 1 - (model$deviance / model$null.deviance)
    
    ## Extract statistics for every variable - preserving order
    allTerms <- lapply(seq_along(terms), function(i){
        var <- names(terms)[i]
        
        if (terms[i] %in% c("factor", "character")) {
            ## Get the factor levels from the model's xlevels (proper order)
            if(var %in% names(model$xlevels)) {
                factor_levels <- model$xlevels[[var]]
                
                ## Create data table with proper levels
                level_counts <- data[!is.na(get(var)), .N, by = var]
                data.table::setnames(level_counts, c("level", "Freq"))
                
                ## Ensure all levels are present
                all_levels_dt <- data.table::data.table(
                                                 level = factor_levels,
                                                 Freq = 0
                                             )
                
                ## Update with actual counts
                for(lev in factor_levels) {
                    if(lev %in% level_counts$level) {
                        all_levels_dt[level == lev, Freq := level_counts[level == lev, Freq]]
                    }
                }
                
                all_levels_dt[, var := var]
                all_levels_dt[, pos := .I]
                all_levels_dt[, var_order := i]
                all_levels_dt[, .(var, level, Freq, pos, var_order)]
            } else {
                ## Fallback for variables not in xlevels
                adf <- data[!is.na(get(var)), .N, by = var]
                data.table::setnames(adf, old = c(var, "N"), new = c("level", "Freq"))
                adf[, var := var]
                adf[, pos := .I]
                adf[, var_order := i]
                adf[, .(var, level, Freq, pos, var_order)]
            }
        }
        else if (terms[i] == "numeric") {
            data.table::data.table(var = var, level = "-", Freq = sum(!is.na(data[[var]])), 
                                   pos = 1, var_order = i)
        }
        else {
            vars = grep(paste0("^", var, "*."), coef$term, value=TRUE)
            data.table::data.table(var = vars, level = "", Freq = nrow(data),
                                   pos = seq_along(vars), var_order = i)
        }
    })
    
    allTermsDF <- data.table::rbindlist(allTerms)
    data.table::setnames(allTermsDF, c("var", "level", "N", "pos", "var_order"))

    ## Add events for binomial models
    if (model$family$family == "binomial") {
        outcome_var <- all.vars(model$formula)[1]
        
        if (!outcome_var %in% names(data)) {
            warning(paste("Outcome variable", outcome_var, "not found in data. Events column will not be created."))
            allTermsDF[, Events := NA_integer_]
        } else {
            outcome_data <- data[[outcome_var]]
            
            if (is.factor(outcome_data)) {
                outcome_binary <- as.numeric(outcome_data) == 2
            } else {
                outcome_binary <- outcome_data
            }
            
            allTermsDF[, Events := {
                if (level == "-") {
                    sum(outcome_binary, na.rm = TRUE)
                } else {
                    if (var %in% names(data)) {
                        sum(outcome_binary[data[[var]] == level], na.rm = TRUE)
                    } else {
                        NA_integer_
                    }
                }
            }, by = seq_len(nrow(allTermsDF))]
        }
    } else {
        allTermsDF[, Events := NA_integer_]
    }

    ## Apply condensing and indenting if requested
    if (condense_table || indent_groups) {
        if (condense_table) {
            indent_groups <- TRUE
        }
        
        allTermsDF[, inds := ifelse(level == "-", var, paste0(var, level))]
        orig_inds_map <- data.table::copy(allTermsDF[, .(var, level, inds, N, Events)])
        
        processed_rows <- list()
        row_counter <- 1
        unique_vars <- unique(allTermsDF[, var])
        
        for (v in unique_vars) {
            var_rows <- allTermsDF[var == v]
            
            if (nrow(var_rows) == 1) {
                processed_rows[[row_counter]] <- var_rows
                row_counter <- row_counter + 1
            } else {
                is_binary <- nrow(var_rows) == 2
                
                if (condense_table && is_binary) {
                    non_ref_row <- var_rows[2]
                    condensed_row <- data.table::copy(non_ref_row)
                    condensed_row[, var := paste0(v, " (", level, ")")]
                    condensed_row[, level := "-"]
                    processed_rows[[row_counter]] <- condensed_row
                    row_counter <- row_counter + 1
                } else {
                    if (indent_groups) {
                        header_row <- data.table::data.table(
                                                      var = v,
                                                      level = "-",
                                                      N = NA_integer_,
                                                      Events = NA_integer_,
                                                      pos = var_rows$pos[1],
                                                      var_order = var_rows$var_order[1],
                                                      inds = NA_character_
                                                  )
                        processed_rows[[row_counter]] <- header_row
                        row_counter <- row_counter + 1
                        
                        for (i in 1:nrow(var_rows)) {
                            group_row <- data.table::copy(var_rows[i])
                            group_row[, var := paste0("    ", level)]
                            group_row[, level := "-"]
                            processed_rows[[row_counter]] <- group_row
                            row_counter <- row_counter + 1
                        }
                    } else {
                        for (i in 1:nrow(var_rows)) {
                            processed_rows[[row_counter]] <- var_rows[i]
                            row_counter <- row_counter + 1
                        }
                    }
                }
            }
        }
        
        allTermsDF <- data.table::rbindlist(processed_rows, fill = TRUE)
        
        for (i in 1:nrow(allTermsDF)) {
            current_var <- allTermsDF$var[i]
            
            if (is.na(allTermsDF$inds[i]) && !grepl("^    ", current_var) && !grepl("\\(", current_var)) {
                next
            }
            
            if (grepl("\\(", current_var)) {
                orig_var <- gsub(" \\(.*\\)", "", current_var)
                orig_level <- gsub(".*\\((.*)\\)", "\\1", current_var)
                
                matching <- orig_inds_map[var == orig_var & level == orig_level]
                if (nrow(matching) > 0) {
                    allTermsDF[i, `:=`(inds = matching$inds[1],
                                       N = matching$N[1],
                                       Events = matching$Events[1])]
                }
            }
            else if (grepl("^    ", current_var)) {
                clean_level <- gsub("^    ", "", current_var)
                
                parent_var <- NA_character_
                for (j in (i-1):1) {
                    if (!grepl("^    ", allTermsDF$var[j]) && allTermsDF$var[j] != "") {
                        parent_var <- gsub(" \\(.*\\)", "", allTermsDF$var[j])
                        break
                    }
                }
                
                if (!is.na(parent_var)) {
                    matching <- orig_inds_map[var == parent_var & level == clean_level]
                    if (nrow(matching) > 0) {
                        allTermsDF[i, `:=`(inds = matching$inds[1],
                                           N = matching$N[1],
                                           Events = matching$Events[1])]
                    }
                }
            }
        }
    } else {
        allTermsDF[, inds := ifelse(level == "-", var, paste0(var, level))]
    }
    
    ## Process coefficients
    coef[, term := gsub(term, pattern = "`", replacement = "")]
    coef[, inds := term]
    
    ## Merge data
    toShow <- merge(allTermsDF, coef, by.x = "inds", by.y = "inds", all.x = TRUE, sort = FALSE)
    
    ## Sort by variable order first, then position within variable
    data.table::setorder(toShow, var_order, pos)
    
    ## Add variable-based shading indicator
    toShow[, shade_group := var_order %% 2]
    
    ## Select columns
    toShow <- toShow[, .(var, level, N, Events, p_value, estimate, conf_low, conf_high, pos, var_order, shade_group)]
    
    ## Format the values based on exponentiate setting
    toShowExpClean <- data.table::copy(toShow)
    
    ## Create formatted columns for display
    if(exponentiate) {
        toShowExpClean[, effect := ifelse(is.na(estimate), 
                                          NA_real_,
                                          exp(estimate))]
        toShowExpClean[, effect_formatted := ifelse(is.na(N) & is.na(estimate),
                                                    "",
                                             ifelse(is.na(estimate), 
                                                    ref_label,
                                                    format(round(exp(estimate), digits), nsmall = digits)))]
        toShowExpClean[, conf_low_formatted := ifelse(is.na(conf_low), 
                                                      NA_character_,
                                                      format(round(exp(conf_low), digits), nsmall = digits))]
        toShowExpClean[, conf_high_formatted := ifelse(is.na(conf_high), 
                                                       NA_character_,
                                                       format(round(exp(conf_high), digits), nsmall = digits))]
    } else {
        toShowExpClean[, effect := estimate]
        toShowExpClean[, effect_formatted := ifelse(is.na(N) & is.na(estimate),
                                                    "",
                                             ifelse(is.na(estimate), 
                                                    ref_label,
                                                    format(round(estimate, digits), nsmall = digits)))]
        toShowExpClean[, conf_low_formatted := ifelse(is.na(conf_low), 
                                                      NA_character_,
                                                      format(round(conf_low, digits), nsmall = digits))]
        toShowExpClean[, conf_high_formatted := ifelse(is.na(conf_high), 
                                                       NA_character_,
                                                       format(round(conf_high, digits), nsmall = digits))]
    }
    
    ## Format p-values
    toShowExpClean[, p_formatted := ifelse(is.na(p_value), 
                                           NA_character_,
                                    ifelse(p_value < 0.001, 
                                           "< 0.001",
                                           format(round(p_value, 3), nsmall = 3)))]
    
    ## Create the combined effect string with expression for italic p
    effect_abbrev <- if(effect_label == "Odds Ratio") "aOR" else if(effect_label == "Risk Ratio") "aRR" else "Coef"

    toShowExpClean[, effect_string_expr := ifelse(
                         is.na(N) & is.na(estimate),
                         "''",
                         fcase(
                             is.na(estimate), paste0("'", ref_label, "'"),
                             
                             p_value < 0.001 & !exponentiate & (conf_low < 0 | conf_high < 0),
                             paste0("'", effect_formatted, " (", conf_low_formatted, ", ", 
                                    conf_high_formatted, "); '*~italic(p)~'< 0.001'"),
                             
                             p_value < 0.001,
                             paste0("'", effect_formatted, " (", conf_low_formatted, "-", 
                                    conf_high_formatted, "); '*~italic(p)~'< 0.001'"),
                             
                             !exponentiate & (conf_low < 0 | conf_high < 0),
                             paste0("'", effect_formatted, " (", conf_low_formatted, ", ", 
                                    conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'"),
                             
                             default = paste0("'", effect_formatted, " (", conf_low_formatted, "-", 
                                              conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'")
                         )
                     )]
    
    ## Format N and events with thousands separator
    toShowExpClean[, n_formatted := ifelse(is.na(N), "", format(N, big.mark = ",", scientific = FALSE))]
    toShowExpClean[, events_formatted := ifelse(is.na(Events), "", format(Events, big.mark = ",", scientific = FALSE))]
    
    ## Clean up variable names for display
    toShowExpClean[, var_display := as.character(var)]
    
    if (indent_groups || condense_table) {
        toShowExpClean[, var_display := var]
        
        for (v in unique(toShowExpClean$var)) {
            if (v != "" && !grepl("^    ", v)) {
                clean_v <- gsub(" \\(.*\\)", "", v)
                
                label <- if (!is.null(var_labels) && clean_v %in% names(var_labels)) {
                             var_labels[clean_v]
                         } else if (!is.null(attr(data[[clean_v]], "label"))) {
                             attr(data[[clean_v]], "label")
                         } else {
                             NULL
                         }
                
                if (!is.null(label)) {
                    if (grepl("\\(", v)) {
                        category <- gsub(".*\\((.*)\\)", "\\1", v)
                        if (category %in% c("1", "Yes", "yes")) {
                            toShowExpClean[var == v, var_display := label]
                        } else {
                            toShowExpClean[var == v, var_display := paste0(label, " (", category, ")")]
                        }
                    } else {
                        toShowExpClean[var == v, var_display := label]
                    }
                }
            }
        }
        
        toShowExpClean[, level := ""]
    } else {
        for(v in unique(toShowExpClean$var)) {
            if(v %in% toShowExpClean$var) {
                if(!is.null(var_labels) && v %in% names(var_labels)) {
                    toShowExpClean[var == v, var_display := var_labels[v]]
                }
                else if(!is.null(attr(data[[v]], "label"))) {
                    toShowExpClean[var == v, var_display := attr(data[[v]], "label")]
                }
            }
        }
    }
    
    toShowExpClean[duplicated(var), var_display := ""]
    
    ## Handle missing estimates for plotting
    if(exponentiate) {
        toShowExpClean[is.na(estimate), estimate := 0]
    } else {
        toShowExpClean[is.na(estimate), estimate := 0]
    }
    
    ## Reorder (flip) - but maintain the variable grouping
    toShowExpClean <- toShowExpClean[order(nrow(toShowExpClean):1)]
    toShowExpClean[, x_pos := .I]
    
    ## Add shade_group based on var_order for alternating by variable
    toShowExpClean[, shade_group := var_order %% 2]
    
    ## Number of rows for later adjustments
    n_rows <- nrow(toShowExpClean)
    
    ## Calculate plot ranges with better handling of extreme cases
    if(exponentiate) {
        rangeb <- range(toShow$conf_low, toShow$conf_high, na.rm = TRUE)
        breaks <- grDevices::axisTicks(rangeb/2, log = TRUE, nint = 7)
        
        min_ci <- min(toShow$conf_low, na.rm = TRUE)
        max_ci <- max(toShow$conf_high, na.rm = TRUE)
        
        is_one_sided <- (min_ci > 0) || (max_ci < 0)
        
        if(!1 %in% breaks) {
            breaks <- sort(unique(c(breaks, 1)))
        }
        
        if(min_ci > 0) {
            rangeb[1] <- log(0.9)
        } else if(max_ci < 0) {
            rangeb[2] <- log(1.1)
        }
        
        breaks <- breaks[breaks >= exp(rangeb[1]) & breaks <= exp(rangeb[2])]
        reference_value <- 1
    } else {
        rangeb <- range(toShow$conf_low, toShow$conf_high, na.rm = TRUE)
        breaks <- pretty(rangeb, n = 7)
        breaks <- breaks[breaks >= rangeb[1] & breaks <= rangeb[2]]
        
        ticks_df <- data.frame(
            x = -0.5,
            xend = -0.7,
            y = breaks,
            yend = breaks
        )
        
        is_one_sided <- (min(rangeb) > 0) || (max(rangeb) < 0)
        
        if(!0 %in% breaks) {
            breaks <- sort(unique(c(breaks, 0)))
        }
        
        if(min(rangeb) > 0) {
            rangeb[1] <- -0.1 * abs(max(rangeb))
        } else if(max(rangeb) < 0) {
            rangeb[2] <- 0.1 * abs(min(rangeb))
        }
        
        reference_value <- 0
    }
    
    ## Calculate dynamic column widths based on content
    max_var_len <- max(nchar(toShowExpClean$var_display) * 0.75, nchar("Variable"), na.rm = TRUE)
    max_n_len <- max(nchar(toShowExpClean$n_formatted), nchar("n"), na.rm = TRUE)
    
    ## For effect column, calculate without expressions
    effect_display_len <- ifelse(toShowExpClean$effect_string_expr == paste0("'", ref_label, "'"),
                                 nchar(ref_label),
                                 nchar(paste0(toShowExpClean$effect_formatted, " (", 
                                              toShowExpClean$conf_low_formatted, "-",
                                              toShowExpClean$conf_high_formatted, "); p = ",
                                              toShowExpClean$p_formatted)))
    max_effect_len <- max(effect_display_len, nchar(paste0(effect_abbrev, " (95% CI); p-value")), na.rm = TRUE) * 1.1
    
    ## Calculate total character width needed
    if (indent_groups || condense_table) {
        total_text_chars <- max_var_len + max_n_len + max_effect_len + 4
    } else {
        max_level_len <- max(nchar(toShowExpClean$level) * 0.75, nchar("Group"), na.rm = TRUE)
        total_text_chars <- max_var_len + max_level_len + max_n_len + max_effect_len + 4
    }
    
    ## Calculate optimal tbl_width if not provided
    if(is.null(tbl_width)) {
        char_to_inch <- 0.08 * font_size
        text_width_needed <- total_text_chars * char_to_inch
        
        if(exponentiate) {
            min_ci_value <- min(toShow$conf_low, na.rm = TRUE)
            max_ci_value <- max(toShow$conf_high, na.rm = TRUE)
        } else {
            min_ci_value <- min(toShow$conf_low, na.rm = TRUE)
            max_ci_value <- max(toShow$conf_high, na.rm = TRUE)
        }
        
        forest_data_range <- max_ci_value - min_ci_value
        forest_width_min <- max(5, forest_data_range * 1.2 + 2)
        
        if(!is.null(plot_width)) {
            available_for_text <- plot_width - forest_width_min
            if(available_for_text > text_width_needed) {
                tbl_width <- text_width_needed / plot_width
            } else {
                tbl_width <- min(0.5, available_for_text / plot_width)
                warning(paste0("Specified plot width (", plot_width, 
                               ") may be too narrow. Consider width of at least ", 
                               round(text_width_needed + forest_width_min, 1), " inches."))
            }
        } else {
            text_proportion_needed <- text_width_needed / (text_width_needed + forest_width_min * 1.2)
            
            if(n_rows < 10) {
                tbl_width <- text_proportion_needed * 0.95
            } else if(n_rows > 30) {
                tbl_width <- text_proportion_needed * 0.85
            } else {
                tbl_width <- text_proportion_needed * 0.9
            }
        }
        
        tbl_width <- max(0.25, min(0.55, tbl_width))
        
        if(is_one_sided) {
            tbl_width <- tbl_width * 0.85
        }
    }
    
    ## Calculate recommended plot dimensions
    rec_height <- max(5, min(20, 3 + n_rows * 0.25))
    
    if(!is.null(plot_width)) {
        rec_width <- plot_width
        if(!is.null(plot_height)) {
            rec_height <- plot_height
        }
    } else {
        char_to_inch <- 0.08 * font_size
        text_width_needed <- total_text_chars * char_to_inch
        
        if(exponentiate) {
            min_ci_value <- min(toShow$conf_low, na.rm = TRUE)
            max_ci_value <- max(toShow$conf_high, na.rm = TRUE)
        } else {
            min_ci_value <- min(toShow$conf_low, na.rm = TRUE) 
            max_ci_value <- max(toShow$conf_high, na.rm = TRUE)
        }
        forest_data_range <- max_ci_value - min_ci_value
        forest_width_needed <- max(6, forest_data_range * 1.2 + 3)
        
        rec_width <- max((text_width_needed / tbl_width) * 1.5,  
                         text_width_needed + forest_width_needed)
        rec_width <- max(12, min(30, rec_width))
    }
    
    ## Provide dimension recommendations
    if(is.null(plot_width) || is.null(plot_height)) {
        message(paste0("Recommended plot dimensions: width = ", round(rec_width, 1), 
                       " inches, height = ", round(rec_height, 1), " inches"))
    }
    
    ## Calculate the range for the plot
    rangeplot <- rangeb
    range_width <- diff(rangeb)
    
    if(is.null(tbl_width)) {
        table_space_factor <- 1.5
    } else {
        table_space_factor <- tbl_width / (1 - tbl_width)
    }
    
    rangeplot[1] <- rangeb[1] - (range_width * table_space_factor * 1.3)
    rangeplot[2] <- rangeb[2] + (range_width * 0.05)
    
    width <- diff(rangeplot)

    ## Handle show_n_events parameter
    if (is.null(show_n_events)) {
        show_n_events <- character(0)
    } else {
        show_n_events[show_n_events == "events"] <- "Events"
    }
    
    ## Start positions from the left edge of the extended range
    y_variable <- rangeplot[1] + col_width_var * width
    
    if (indent_groups || condense_table) {
        y_level <- y_variable + (col_width_var + 0.05) * width
        next_base <- y_variable
    } else {
        y_level <- y_variable + col_width_level * width
        next_base <- y_level
    }
    
    if ("n" %in% show_n_events) {
        y_n <- y_level + col_width_n * width
        next_pos <- y_n
    } else {
        next_pos <- y_level
    }
    
    if ("Events" %in% show_n_events) {
        y_events <- next_pos + col_width_events * width
        y_or <- y_events + col_width_or * width
    } else {
        y_or <- next_pos + col_width_or * width
    }

    forest_start <- rangeb[1] - 0.05 * range_width
    
    annot_font <- font_size * annot_size
    header_font <- font_size * header_size
    
    ## Custom ticks data
    ticks_df <- data.frame(
        x = -0.5,
        xend = -0.7,
        y = breaks,
        yend = breaks
    )
    
    ## Format deviance values
    null_dev_formatted <- format(round(gmodel$null_deviance, 2), nsmall = 2)
    resid_dev_formatted <- format(round(gmodel$residual_deviance, 2), nsmall = 2)
    pseudo_r2_formatted <- format(round(gmodel$pseudo_r2, 3), nsmall = 3)
    
    ## Create the plot
    if(exponentiate) {
        p <- ggplot2::ggplot(toShowExpClean, ggplot2::aes(x_pos, exp(estimate))) +
            
            ggplot2::geom_rect(ggplot2::aes(xmin = x_pos - .5, xmax = x_pos + .5,
                                            ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                                            fill = ordered(shade_group + 1))) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
            ggplot2::scale_size_continuous(range = c(1, 6), guide = "none") +
            ggplot2::scale_fill_manual(values = c("#FFFFFF", "#EEEEEE"), guide = "none") +
            
            ggplot2::geom_point(ggplot2::aes(size = N), pch = 22, color = "#000000", fill = color) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = exp(conf_low), ymax = exp(conf_high)), width = 0.15) +
            
            ggplot2::annotate(geom = "segment",
                              x = -0.5, xend = -0.5,
                              y = exp(rangeb[1]),
                              yend = exp(rangeb[2]),
                              color = "#000000", linewidth = 1) +
            
            ggplot2::annotate(geom = "segment", 
                              x = -0.5, xend = max(toShowExpClean$x_pos) + 0.5, 
                              y = reference_value, yend = reference_value, 
                              linetype = "longdash") +
            
            ggplot2::geom_segment(data = ticks_df,
                                  ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                                  inherit.aes = FALSE,
                                  color = "#000000",
                                  linewidth = 1) +
            ggplot2::geom_text(data = ticks_df,
                               ggplot2::aes(x = xend - 0.05, y = y, label = sprintf("%g", y)),
                               inherit.aes = FALSE,
                               hjust = 0.5,
                               vjust = 1.3,
                               size = annot_font * 1.5) +
            
            ggplot2::coord_flip(ylim = exp(rangeplot)) +
            ggplot2::ggtitle(title) +
            ggplot2::scale_y_log10(name = effect_label,
                                   labels = sprintf("%g", breaks),
                                   expand = c(0.02, 0.02),
                                   breaks = breaks) +
            ggplot2::theme_light() +
            ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0),
                           panel.grid.minor.y = ggplot2::element_blank(),
                           panel.grid.minor.x = ggplot2::element_blank(),
                           panel.grid.major.y = ggplot2::element_blank(),
                           panel.grid.major.x = ggplot2::element_blank(),
                           legend.position = "none",
                           panel.border = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank(),
                           axis.ticks.y = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(),
                           plot.title = ggplot2::element_text(size = font_size * title_size, face = "bold", hjust = 0.5)) +
            ggplot2::xlab("") +
            
            ## Variable column
            ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.5, y = exp(y_variable),
                              label = "Variable", fontface = "bold", hjust = 0,
                              size = header_font) +
            
            {if (indent_groups || condense_table) {
                 ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = exp(y_variable),
                                   label = toShowExpClean$var_display, 
                                   fontface = ifelse(grepl("^    ", toShowExpClean$var_display), "plain", "bold"), 
                                   hjust = 0,
                                   size = annot_font)
             } else {
                 ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = exp(y_variable),
                                   label = toShowExpClean$var_display, fontface = "bold", hjust = 0,
                                   size = annot_font)
             }} +
            
            ## Group/level column
            {if (!(indent_groups || condense_table)) {
                 list(
                     ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.5, y = exp(y_level),
                                       label = "Group", fontface = "bold", hjust = 0,
                                       size = header_font),
                     ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = exp(y_level),
                                       label = toShowExpClean$level, hjust = 0,
                                       size = annot_font)
                 )
             }} +
            
            ## N column (conditional)
            {if ("n" %in% show_n_events) {
                 list(
                     ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.5, y = exp(y_n),
                                       label = "n", fontface = "bold.italic", hjust = 0.5,
                                       size = header_font),
                     ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = exp(y_n),
                                       label = toShowExpClean$n_formatted, hjust = 0.5,
                                       size = annot_font)
                 )
             }} +
            
            ## Events column (conditional)
            {if ("Events" %in% show_n_events) {
                 list(
                     ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.5, y = exp(y_events),
                                       label = "Events", fontface = "bold", hjust = 0.5,
                                       size = header_font),
                     ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = exp(y_events),
                                       label = toShowExpClean$events_formatted, hjust = 0.5,
                                       size = annot_font)
                 )
             }} +
            
            ## Effect column
            ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.4, y = exp(y_or),
                              label = paste0("bold('", effect_abbrev, " (95% CI); '*bolditalic(p)*'-value')"),
                              hjust = 0, size = header_font, parse = TRUE) +
            
            ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = exp(y_or),
                              label = toShowExpClean$effect_string_expr, hjust = 0,
                              size = annot_font, parse = TRUE) +
            
            ## X-axis label
            ggplot2::annotate(geom = "text", x = -1.5, y = reference_value,
                              label = effect_label, fontface = "bold",
                              hjust = 0.5, vjust = 2, size = annot_font * 1.5) +
            
            ## Model statistics footer
            ggplot2::annotate(geom = "text", x = 0.5, y = exp(y_variable),
                              label = paste0("Observations: ", gmodel$nobs_formatted,
                                             "\nModel: ", gmodel$family,
                                             "\nNull (Residual) Deviance: ", null_dev_formatted, " (", resid_dev_formatted, ")",
                                             "\nPseudo R²: ", pseudo_r2_formatted,
                                             "\nAIC: ", gmodel$AIC_formatted),
                              size = annot_font * 0.8, hjust = 0, vjust = 1.2, fontface = "italic")
        
    } else {
        ## Non-exponentiated plot (linear scale)
        p <- ggplot2::ggplot(toShowExpClean, ggplot2::aes(x_pos, estimate)) +
            
            ggplot2::geom_rect(ggplot2::aes(xmin = x_pos - .5, xmax = x_pos + .5,
                                            ymin = rangeplot[1], ymax = rangeplot[2],
                                            fill = ordered(shade_group + 1))) +
            ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
            ggplot2::scale_size_continuous(range = c(1, 6), guide = "none") +
            ggplot2::scale_fill_manual(values = c("#FFFFFF", "#EEEEEE"), guide = "none") +
            
            ggplot2::geom_point(ggplot2::aes(size = N), pch = 22, color = "#000000", fill = color) +
            ggplot2::geom_errorbar(ggplot2::aes(ymin = conf_low, ymax = conf_high), width = 0.15) +
            
            ggplot2::annotate(geom = "segment",
                              x = -0.5, xend = -0.5,
                              y = rangeb[1],
                              yend = rangeb[2],
                              color = "#000000", linewidth = 1) +
            
            ggplot2::annotate(geom = "segment", 
                              x = -0.5, xend = max(toShowExpClean$x_pos) + 0.5, 
                              y = reference_value, yend = reference_value, 
                              linetype = "longdash") +
            
            ggplot2::geom_segment(data = ticks_df,
                                  ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
                                  inherit.aes = FALSE,
                                  color = "#000000",
                                  linewidth = 1) +
            ggplot2::geom_text(data = ticks_df,
                               ggplot2::aes(x = xend - 0.05, y = y, label = sprintf("%g", y)),
                               inherit.aes = FALSE,
                               hjust = 0.5,
                               vjust = 1.3,
                               size = annot_font * 1.5) +
            
            ggplot2::coord_flip(ylim = rangeplot) +
            ggplot2::ggtitle(title) +
            ggplot2::scale_y_continuous(name = effect_label,
                                        labels = sprintf("%g", breaks),
                                        expand = c(0.02, 0.02),
                                        breaks = breaks) +
            ggplot2::theme_light() +
            ggplot2::theme(plot.margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0),
                           panel.grid.minor.y = ggplot2::element_blank(),
                           panel.grid.minor.x = ggplot2::element_blank(),
                           panel.grid.major.y = ggplot2::element_blank(),
                           panel.grid.major.x = ggplot2::element_blank(),
                           legend.position = "none",
                           panel.border = ggplot2::element_blank(),
                           axis.title.y = ggplot2::element_blank(),
                           axis.text.y = ggplot2::element_blank(),
                           axis.ticks.y = ggplot2::element_blank(),
                           axis.title.x = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank(),
                           plot.title = ggplot2::element_text(size = font_size * title_size, face = "bold", hjust = 0.5)) +
            ggplot2::xlab("") +
            
            ## Variable column
            ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.5, y = y_variable,
                              label = "Variable", fontface = "bold", hjust = 0,
                              size = header_font) +
            
            {if (indent_groups || condense_table) {
                 ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = y_variable,
                                   label = toShowExpClean$var_display, 
                                   fontface = ifelse(grepl("^    ", toShowExpClean$var_display), "plain", "bold"), 
                                   hjust = 0,
                                   size = annot_font)
             } else {
                 ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = y_variable,
                                   label = toShowExpClean$var_display, fontface = "bold", hjust = 0,
                                   size = annot_font)
             }} +
            
            ## Group/level column
            {if (!(indent_groups || condense_table)) {
                 list(
                     ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.5, y = y_level,
                                       label = "Group", fontface = "bold", hjust = 0,
                                       size = header_font),
                     ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = y_level,
                                       label = toShowExpClean$level, hjust = 0,
                                       size = annot_font)
                 )
             }} +
            
            ## N column (conditional)
            {if ("n" %in% show_n_events) {
                 list(
                     ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.5, y = y_n,
                                       label = "n", fontface = "bold.italic", hjust = 0.5,
                                       size = header_font),
                     ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = y_n,
                                       label = toShowExpClean$n_formatted, hjust = 0.5,
                                       size = annot_font)
                 )
             }} +
            
            ## Events column (conditional)
            {if ("Events" %in% show_n_events) {
                 list(
                     ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.5, y = y_events,
                                       label = "Events", fontface = "bold", hjust = 0.5,
                                       size = header_font),
                     ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = y_events,
                                       label = toShowExpClean$events_formatted, hjust = 0.5,
                                       size = annot_font)
                 )
             }} +
            
            ## Effect column
            ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.4, y = y_or,
                              label = paste0("bold('", effect_abbrev, " (95% CI); '*bolditalic(p)*'-value')"),
                              hjust = 0, size = header_font, parse = TRUE) +
            
            ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = y_or,
                              label = toShowExpClean$effect_string_expr, hjust = 0,
                              size = annot_font, parse = TRUE) +
            
            ## X-axis label
            ggplot2::annotate(geom = "text", x = -1.5, y = reference_value,
                              label = effect_label, fontface = "bold",
                              hjust = 0.5, vjust = 2, size = annot_font * 1.5) +
            
            ## Model statistics at bottom
            ggplot2::annotate(geom = "text", x = 0.5, y = y_variable,
                              label = paste0("Observations: ", gmodel$nobs_formatted,
                                             "\nModel: ", gmodel$family,
                                             "\nNull (Residual) Deviance: ", null_dev_formatted, " (", resid_dev_formatted, ")",
                                             "\nPseudo R²: ", pseudo_r2_formatted,
                                             "\nAIC: ", gmodel$AIC_formatted),
                              size = annot_font * 0.8, hjust = 0, vjust = 1.2, fontface = "italic")
    }
    
    ## Add recommended dimensions as an attribute
    attr(p, "recommended_dims") <- list(width = rec_width, height = rec_height)
    
    ## Return the plot
    return(p)
}
