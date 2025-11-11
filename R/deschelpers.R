#' Process variable wrapper - returns both formatted and raw
#' @keywords internal
process_variable <- function(data, var, group_var = NULL, 
                             stats_continuous, stats_categorical,
                             digits, na_include, na_label,
                             test, test_continuous, test_categorical,
                             total, total_label, var_labels, na_percent, ...) {
    
    ## Get variable label
    var_label <- if (!is.null(var_labels) && var %chin% names(var_labels)) {
                     var_labels[var]
                 } else {
                     var
                 }
    
    ## Determine variable type and process accordingly
    if (grepl("^Surv\\(", var)) {
        return(process_survival(
            data = data,
            var = var,
            var_label = var_label,
            group_var = group_var,
            digits = digits,
            na_include = na_include,
            na_label = na_label,
            test = test,
            total = total,
            total_label = total_label,
            ...
        ))
    } else if (is.numeric(data[[var]])) {
        return(process_continuous(
            data = data,
            var = var,
            var_label = var_label,
            group_var = group_var,
            stats = stats_continuous,
            digits = digits,
            na_include = na_include,
            na_label = na_label,
            test = test,
            test_type = test_continuous,
            total = total,
            total_label = total_label,
            ...
        ))
    } else {
        return(process_categorical(
            data = data,
            var = var,
            var_label = var_label,
            group_var = group_var,
            stats = stats_categorical,
            na_include = na_include,
            na_label = na_label,
            test = test,
            test_type = test_categorical,
            total = total,
            total_label = total_label,
            na_percent = na_percent,
            ...
        ))
    }
}

#' Process continuous variable
#' @keywords internal
process_continuous <- function(data, var, var_label, group_var, stats, digits,
                               na_include, na_label, test, test_type,
                               total, total_label, ...) {
    
    ## Pre-allocate lists for better performance
    n_stats <- length(stats)
    formatted_list <- vector("list", n_stats + if(na_include && any(is.na(data[[var]]))) 1 else 0)
    raw_list <- vector("list", n_stats + if(na_include && any(is.na(data[[var]]))) 1 else 0)
    
    if (!is.null(group_var)) {
        if (is.factor(data[[group_var]])) {
            groups <- levels(data[[group_var]])
        } else {
            groups <- unique(data[[group_var]])
            groups <- groups[!is.na(groups)]
        }
        
        first_stat <- TRUE
        p_values <- list()
        
        ## Calculate p-values first if testing
        if (test) {
            for (stat_type in stats) {
                if (stat_type == "range") {
                    p_values[[stat_type]] <- NULL
                } else {
                    p_values[[stat_type]] <- perform_continuous_test(
                        data, var, group_var, test_type, stat_type, ...
                    )
                }
            }
        }
        
        for (i in seq_along(stats)) {
            stat_type <- stats[i]
            
            ## Initialize rows
            formatted_row <- list(
                variable = if (first_stat) var_label else "",
                level = get_stat_label(stat_type)
            )
            
            raw_row <- list(
                variable = if (first_stat) var_label else "",
                level = stat_type,
                stat_type = stat_type
            )
            
            ## Add total column
            if (!isFALSE(total)) {
                values <- data[[var]][!is.na(data[[var]])]
                
                ## Formatted
                formatted_row[[total_label]] <- calc_continuous_stat(values, stat_type, digits)
                
                ## Raw - main value plus supplementary columns
                if (stat_type == "mean_sd") {
                    raw_row[[total_label]] <- mean(values)
                    raw_row[[paste0(total_label, "_sd")]] <- sd(values)
                    raw_row[[paste0(total_label, "_n")]] <- length(values)
                } else if (stat_type == "median_iqr") {
                    raw_row[[total_label]] <- median(values)
                    raw_row[[paste0(total_label, "_q1")]] <- quantile(values, 0.25)
                    raw_row[[paste0(total_label, "_q3")]] <- quantile(values, 0.75)
                    raw_row[[paste0(total_label, "_n")]] <- length(values)
                } else if (stat_type == "median_range") {
                    raw_row[[total_label]] <- median(values)
                    raw_row[[paste0(total_label, "_min")]] <- min(values)
                    raw_row[[paste0(total_label, "_max")]] <- max(values)
                    raw_row[[paste0(total_label, "_n")]] <- length(values)
                } else if (stat_type == "range") {
                    raw_row[[total_label]] <- min(values)  # Store min as main value
                    raw_row[[paste0(total_label, "_max")]] <- max(values)
                    raw_row[[paste0(total_label, "_n")]] <- length(values)
                }
            }
            
            ## Add group columns
            for (grp in groups) {
                grp_data <- data[get(group_var) == grp, get(var)]
                grp_data <- grp_data[!is.na(grp_data)]
                grp_col <- as.character(grp)
                
                ## Formatted
                formatted_row[[grp_col]] <- calc_continuous_stat(grp_data, stat_type, digits)
                
                ## Raw with supplementary columns
                if (length(grp_data) > 0) {
                    if (stat_type == "mean_sd") {
                        raw_row[[grp_col]] <- mean(grp_data)
                        raw_row[[paste0(grp_col, "_sd")]] <- sd(grp_data)
                        raw_row[[paste0(grp_col, "_n")]] <- length(grp_data)
                    } else if (stat_type == "median_iqr") {
                        raw_row[[grp_col]] <- median(grp_data)
                        raw_row[[paste0(grp_col, "_q1")]] <- quantile(grp_data, 0.25)
                        raw_row[[paste0(grp_col, "_q3")]] <- quantile(grp_data, 0.75)
                        raw_row[[paste0(grp_col, "_n")]] <- length(grp_data)
                    } else if (stat_type == "median_range") {
                        raw_row[[grp_col]] <- median(grp_data)
                        raw_row[[paste0(grp_col, "_min")]] <- min(grp_data)
                        raw_row[[paste0(grp_col, "_max")]] <- max(grp_data)
                        raw_row[[paste0(grp_col, "_n")]] <- length(grp_data)
                    } else if (stat_type == "range") {
                        raw_row[[grp_col]] <- min(grp_data)
                        raw_row[[paste0(grp_col, "_max")]] <- max(grp_data)
                        raw_row[[paste0(grp_col, "_n")]] <- length(grp_data)
                    }
                } else {
                    raw_row[[grp_col]] <- NA_real_
                    raw_row[[paste0(grp_col, "_n")]] <- 0
                }
            }
            
            ## Add p_value
            if (test && stat_type != "range" && !is.null(p_values[[stat_type]])) {
                formatted_row[["p_value"]] <- p_values[[stat_type]]
                raw_row[["p_value"]] <- p_values[[stat_type]]
            }
            
            formatted_list[[i]] <- data.table::as.data.table(formatted_row)
            raw_list[[i]] <- data.table::as.data.table(raw_row)
            first_stat <- FALSE
        }

        ## Missing count handling for grouped data
        if (na_include && any(is.na(data[[var]]))) {
            miss_formatted <- list(
                variable = "",
                level = na_label
            )
            
            miss_raw <- list(
                variable = "",
                level = na_label,
                stat_type = "missing"
            )
            
            ## Add total column
            if (!isFALSE(total)) {
                n_miss <- sum(is.na(data[[var]]))
                
                ## Just show count, not percentage for continuous
                miss_formatted[[total_label]] <- if (n_miss >= 1000) {
                                                     format(n_miss, big.mark = ",")
                                                 } else {
                                                     as.character(n_miss)
                                                 }
                miss_raw[[total_label]] <- n_miss
            }
            
            ## Add group columns
            for (grp in groups) {
                grp_col <- as.character(grp)
                n_miss <- sum(is.na(data[[var]]) & data[[group_var]] == grp, na.rm = TRUE)
                
                miss_formatted[[grp_col]] <- if (n_miss >= 1000) {
                                                 format(n_miss, big.mark = ",")
                                             } else {
                                                 as.character(n_miss)
                                             }
                miss_raw[[grp_col]] <- n_miss
            }
            
            formatted_list[[n_stats + 1]] <- data.table::as.data.table(miss_formatted)
            raw_list[[n_stats + 1]] <- data.table::as.data.table(miss_raw)
        }
        
        ## Combine results using rbindlist
        formatted_result <- data.table::rbindlist(formatted_list, fill = TRUE)
        raw_result <- data.table::rbindlist(raw_list, fill = TRUE)
        
    } else {
        ## Ungrouped version (no group_var)
        first_stat <- TRUE
        
        for (i in seq_along(stats)) {
            stat_type <- stats[i]
            
            formatted_row <- list(
                variable = if (first_stat) var_label else "",
                level = get_stat_label(stat_type)
            )
            
            raw_row <- list(
                variable = if (first_stat) var_label else "",
                level = stat_type,
                stat_type = stat_type
            )
            
            ## Add total column
            if (!isFALSE(total)) {
                values <- data[[var]][!is.na(data[[var]])]
                formatted_row[[total_label]] <- calc_continuous_stat(values, stat_type, digits)
                
                ## Raw data
                if (stat_type == "mean_sd") {
                    raw_row[[total_label]] <- mean(values)
                    raw_row[[paste0(total_label, "_sd")]] <- sd(values)
                    raw_row[[paste0(total_label, "_n")]] <- length(values)
                } else if (stat_type == "median_iqr") {
                    raw_row[[total_label]] <- median(values)
                    raw_row[[paste0(total_label, "_q1")]] <- quantile(values, 0.25)
                    raw_row[[paste0(total_label, "_q3")]] <- quantile(values, 0.75)
                    raw_row[[paste0(total_label, "_n")]] <- length(values)
                } else if (stat_type == "median_range") {
                    raw_row[[total_label]] <- median(values)
                    raw_row[[paste0(total_label, "_min")]] <- min(values)
                    raw_row[[paste0(total_label, "_max")]] <- max(values)
                    raw_row[[paste0(total_label, "_n")]] <- length(values)
                } else if (stat_type == "range") {
                    raw_row[[total_label]] <- min(values)
                    raw_row[[paste0(total_label, "_max")]] <- max(values)
                    raw_row[[paste0(total_label, "_n")]] <- length(values)
                }
            }
            
            formatted_list[[i]] <- data.table::as.data.table(formatted_row)
            raw_list[[i]] <- data.table::as.data.table(raw_row)
            first_stat <- FALSE
        }
        
        ## Missing count handling for ungrouped data
        if (na_include && any(is.na(data[[var]]))) {
            n_miss <- sum(is.na(data[[var]]))
            
            miss_formatted <- list(
                variable = "",
                level = na_label
            )
            
            miss_raw <- list(
                variable = "",
                level = na_label,
                stat_type = "missing"
            )
            
            if (!isFALSE(total)) {
                miss_formatted[[total_label]] <- if (n_miss >= 1000) {
                                                     format(n_miss, big.mark = ",")
                                                 } else {
                                                     as.character(n_miss)
                                                 }
                miss_raw[[total_label]] <- n_miss
            }
            
            formatted_list[[n_stats + 1]] <- data.table::as.data.table(miss_formatted)
            raw_list[[n_stats + 1]] <- data.table::as.data.table(miss_raw)
        }
        
        formatted_result <- data.table::rbindlist(formatted_list, fill = TRUE)
        raw_result <- data.table::rbindlist(raw_list, fill = TRUE)
    }
    
    return(list(formatted = formatted_result, raw = raw_result))
}

#' Process categorical variable
#' @keywords internal
process_categorical <- function(data, var, var_label, group_var, stats,
                                na_include, na_label, test, test_type,
                                total, total_label, na_percent, ...) {
    
    ## Get levels
    if (is.factor(data[[var]])) {
        levels_to_show <- levels(data[[var]])
    } else {
        levels_to_show <- unique(data[[var]])
        levels_to_show <- levels_to_show[!is.na(levels_to_show)]
    }
    
    ## Add NA level if requested
    if (na_include && any(is.na(data[[var]]))) {
        levels_to_show <- c(levels_to_show, NA)
    }
    
    ## Pre-allocate lists
    n_levels <- length(levels_to_show)
    formatted_list <- vector("list", n_levels)  # No +1, no separate header row
    raw_list <- vector("list", n_levels)
    
    if (!is.null(group_var)) {
        if (is.factor(data[[group_var]])) {
            groups <- levels(data[[group_var]])
        } else {
            groups <- unique(data[[group_var]])
            groups <- groups[!is.na(groups)]
        }
        
        ## Calculate p-value first
        p_value <- if (test) {
                       perform_categorical_test(data, var, group_var, test_type, ...)
                   } else {
                       NULL
                   }
        
        ## Level rows - first level gets the variable name
        for (i in seq_along(levels_to_show)) {
            lvl <- levels_to_show[i]
            
            level_formatted <- list(
                variable = if (i == 1) var_label else "",  # Variable name on first level only
                level = if (is.na(lvl)) na_label else as.character(lvl)
            )
            
            level_raw <- list(
                variable = if (i == 1) var_label else "",
                level = if (is.na(lvl)) na_label else as.character(lvl),
                stat_type = "category"
            )
            
            ## Calculate total denominator
            if (na_percent || is.na(lvl)) {
                total_denom <- nrow(data)
            } else {
                total_denom <- sum(!is.na(data[[var]]))
            }
            
            ## Add total column
            if (!isFALSE(total)) {
                if (is.na(lvl)) {
                    n <- sum(is.na(data[[var]]))
                } else {
                    n <- sum(data[[var]] == lvl, na.rm = TRUE)
                }
                
                level_formatted[[total_label]] <- format_categorical_stat(n, total_denom, stats)
                level_raw[[total_label]] <- n
                level_raw[[paste0(total_label, "_total")]] <- total_denom
            }
            
            ## Add group columns
            for (grp in groups) {
                grp_col <- as.character(grp)
                
                ## Calculate group denominator
                if (na_percent || is.na(lvl)) {
                    grp_denom <- sum(data[[group_var]] == grp, na.rm = TRUE)
                } else {
                    grp_denom <- sum(!is.na(data[[var]]) & data[[group_var]] == grp, na.rm = TRUE)
                }
                
                if (is.na(lvl)) {
                    n <- sum(is.na(data[[var]]) & data[[group_var]] == grp, na.rm = TRUE)
                } else {
                    n <- sum(data[[var]] == lvl & data[[group_var]] == grp, na.rm = TRUE)
                }
                
                level_formatted[[grp_col]] <- format_categorical_stat(n, grp_denom, stats)
                level_raw[[grp_col]] <- n
                level_raw[[paste0(grp_col, "_total")]] <- grp_denom
            }
            
            ## Add p-value to first level only
            if (i == 1 && test && !is.null(p_value)) {
                level_formatted[["p_value"]] <- p_value
                level_raw[["p_value"]] <- p_value
            }
            
            formatted_list[[i]] <- data.table::as.data.table(level_formatted)
            raw_list[[i]] <- data.table::as.data.table(level_raw)
        }
        
        formatted_result <- data.table::rbindlist(formatted_list, fill = TRUE)
        raw_result <- data.table::rbindlist(raw_list, fill = TRUE)
        
    } else {
        ## Ungrouped version - first level gets variable name
        ## Calculate denominator
        if (na_percent) {
            total_denom <- nrow(data)
        } else {
            total_denom <- sum(!is.na(data[[var]]))
        }
        
        for (i in seq_along(levels_to_show)) {
            lvl <- levels_to_show[i]
            
            level_formatted <- list(
                variable = if (i == 1) var_label else "",  # Variable name on first level only
                level = if (is.na(lvl)) na_label else as.character(lvl)
            )
            
            level_raw <- list(
                variable = if (i == 1) var_label else "",
                level = if (is.na(lvl)) na_label else as.character(lvl),
                stat_type = "category"
            )
            
            if (!isFALSE(total)) {
                if (is.na(lvl)) {
                    n <- sum(is.na(data[[var]]))
                } else {
                    n <- sum(data[[var]] == lvl, na.rm = TRUE)
                }
                
                level_formatted[[total_label]] <- format_categorical_stat(n, total_denom, stats)
                level_raw[[total_label]] <- n
                level_raw[[paste0(total_label, "_total")]] <- total_denom
            }
            
            formatted_list[[i]] <- data.table::as.data.table(level_formatted)
            raw_list[[i]] <- data.table::as.data.table(level_raw)
        }
        
        formatted_result <- data.table::rbindlist(formatted_list, fill = TRUE)
        raw_result <- data.table::rbindlist(raw_list, fill = TRUE)
    }
    
    return(list(formatted = formatted_result, raw = raw_result))
}

#' Process survival variable
#' @keywords internal
process_survival <- function(data, var, var_label, group_var, digits,
                             na_include, na_label, test, total, total_label, ...) {
    
    ## Parse Surv() expression
    surv_match <- regexec("Surv\\(([^,]+),\\s*([^)]+)\\)", var)
    surv_parts <- regmatches(var, surv_match)[[1]]
    
    if (length(surv_parts) < 3) {
        stop("Invalid Surv() syntax: ", var)
    }
    
    time_var <- trimws(surv_parts[2])
    status_var <- trimws(surv_parts[3])
    
    ## Check if survival package available
    if (!requireNamespace("survival", quietly = TRUE)) {
        stop("Package 'survival' required for survival analysis")
    }
    
    formatted_result <- data.table::data.table()
    raw_result <- data.table::data.table()
    
    if (!is.null(group_var)) {
        if (is.factor(data[[group_var]])) {
            groups <- levels(data[[group_var]])
        } else {
            groups <- unique(data[[group_var]])
            groups <- groups[!is.na(groups)]
        }
        
        ## Calculate p-value (log-rank test)
        p_value <- if (test) {
                       clean_data <- data[!is.na(get(time_var)) & !is.na(get(status_var)) & !is.na(get(group_var))]
                       if (nrow(clean_data) > 0) {
                           formula <- survival::Surv(get(time_var), get(status_var)) ~ get(group_var)
                           tryCatch({
                               survival::survdiff(formula, data = clean_data)$pvalue
                           }, error = function(e) NA_real_)
                       } else {
                           NA_real_
                       }
                   } else {
                       NULL
                   }
        
        ## Variable name and median row
        formatted_row <- list(
            variable = var_label,
            level = "Median (95% CI)"
        )
        
        raw_row <- list(
            variable = var_label,
            level = "median",
            stat_type = "survival"
        )
        
        ## Add total column
        if (!isFALSE(total)) {
            surv_obj <- survival::Surv(data[[time_var]], data[[status_var]])
            fit <- survival::survfit(surv_obj ~ 1)
            
            median_surv <- summary(fit)$table["median"]
            ci_lower <- summary(fit)$table["0.95LCL"]
            ci_upper <- summary(fit)$table["0.95UCL"]
            
            formatted_row[[total_label]] <- sprintf(
                paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
                median_surv, ci_lower, ci_upper
            )
            
            raw_row[[total_label]] <- median_surv
            raw_row[[paste0(total_label, "_ci_lower")]] <- ci_lower
            raw_row[[paste0(total_label, "_ci_upper")]] <- ci_upper
        }
        
        ## Add group columns
        for (grp in groups) {
            grp_col <- as.character(grp)
            grp_data <- data[get(group_var) == grp]
            
            surv_obj <- survival::Surv(grp_data[[time_var]], grp_data[[status_var]])
            fit <- survival::survfit(surv_obj ~ 1)
            
            median_surv <- summary(fit)$table["median"]
            ci_lower <- summary(fit)$table["0.95LCL"]
            ci_upper <- summary(fit)$table["0.95UCL"]
            
            formatted_row[[grp_col]] <- sprintf(
                paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
                median_surv, ci_lower, ci_upper
            )
            
            raw_row[[grp_col]] <- median_surv
            raw_row[[paste0(grp_col, "_ci_lower")]] <- ci_lower
            raw_row[[paste0(grp_col, "_ci_upper")]] <- ci_upper
        }
        
        ## Add p-value
        if (test && !is.null(p_value)) {
            formatted_row[["p_value"]] <- p_value
            raw_row[["p_value"]] <- p_value
        }
        
        formatted_result <- data.table::as.data.table(formatted_row)
        raw_result <- data.table::as.data.table(raw_row)
        
    } else {
        ## Ungrouped version
        formatted_row <- list(
            variable = var_label,
            level = "Median (95% CI)"
        )
        
        raw_row <- list(
            variable = var_label,
            level = "median",
            stat_type = "survival"
        )
        
        if (!isFALSE(total)) {
            surv_obj <- survival::Surv(data[[time_var]], data[[status_var]])
            fit <- survival::survfit(surv_obj ~ 1)
            
            median_surv <- summary(fit)$table["median"]
            ci_lower <- summary(fit)$table["0.95LCL"]
            ci_upper <- summary(fit)$table["0.95UCL"]
            
            formatted_row[[total_label]] <- sprintf(
                paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
                median_surv, ci_lower, ci_upper
            )
            
            raw_row[[total_label]] <- median_surv
            raw_row[[paste0(total_label, "_ci_lower")]] <- ci_lower
            raw_row[[paste0(total_label, "_ci_upper")]] <- ci_upper
        }
        
        formatted_result <- data.table::as.data.table(formatted_row)
        raw_result <- data.table::as.data.table(raw_row)
    }
    
    return(list(formatted = formatted_result, raw = raw_result))
}

#' Calculate continuous statistic
#' @keywords internal
calc_continuous_stat <- function(x, stat_type, digits) {
    if (length(x) == 0) return("")
    
    switch(stat_type,
           "mean_sd" = {
               mean_val <- mean(x)
               sd_val <- sd(x)
               mean_str <- if (abs(mean_val) >= 1000) {
                               format(round(mean_val, digits), big.mark = ",")
                           } else {
                               sprintf(paste0("%.", digits, "f"), mean_val)
                           }
               sd_str <- if (abs(sd_val) >= 1000) {
                             format(round(sd_val, digits), big.mark = ",")
                         } else {
                             sprintf(paste0("%.", digits, "f"), sd_val)
                         }
               paste0(mean_str, " +/- ", sd_str)
           },
           "median_iqr" = {
               med_val <- median(x)
               q1_val <- quantile(x, 0.25)
               q3_val <- quantile(x, 0.75)
               med_str <- if (abs(med_val) >= 1000) {
                              format(round(med_val, digits), big.mark = ",")
                          } else {
                              sprintf(paste0("%.", digits, "f"), med_val)
                          }
               q1_str <- if (abs(q1_val) >= 1000) {
                             format(round(q1_val, digits), big.mark = ",")
                         } else {
                             sprintf(paste0("%.", digits, "f"), q1_val)
                         }
               q3_str <- if (abs(q3_val) >= 1000) {
                             format(round(q3_val, digits), big.mark = ",")
                         } else {
                             sprintf(paste0("%.", digits, "f"), q3_val)
                         }
               paste0(med_str, " [", q1_str, "-", q3_str, "]")
           },
           "median_range" = {
               med_val <- median(x)
               min_val <- min(x)
               max_val <- max(x)
               med_str <- if (abs(med_val) >= 1000) {
                              format(round(med_val, digits), big.mark = ",")
                          } else {
                              sprintf(paste0("%.", digits, "f"), med_val)
                          }
               min_str <- if (abs(min_val) >= 1000) {
                              format(round(min_val, digits), big.mark = ",")
                          } else {
                              sprintf(paste0("%.", digits, "f"), min_val)
                          }
               max_str <- if (abs(max_val) >= 1000) {
                              format(round(max_val, digits), big.mark = ",")
                          } else {
                              sprintf(paste0("%.", digits, "f"), max_val)
                          }
               ## Use comma separator for range with negative numbers
               if (min_val < 0 || max_val < 0) {
                   paste0(med_str, " (", min_str, ", ", max_str, ")")
               } else {
                   paste0(med_str, " (", min_str, "-", max_str, ")")
               }
           },
           "range" = {
               min_val <- min(x)
               max_val <- max(x)
               min_str <- if (abs(min_val) >= 1000) {
                              format(round(min_val, digits), big.mark = ",")
                          } else {
                              sprintf(paste0("%.", digits, "f"), min_val)
                          }
               max_str <- if (abs(max_val) >= 1000) {
                              format(round(max_val, digits), big.mark = ",")
                          } else {
                              sprintf(paste0("%.", digits, "f"), max_val)
                          }
               ## Use "to" for ranges with negative numbers
               if (min_val < 0 || max_val < 0) {
                   paste0(min_str, " to ", max_str)
               } else {
                   paste0(min_str, "-", max_str)
               }
           },
           sprintf(paste0("%.", digits, "f"), mean(x))  # Default
           )
}

#' Get statistic label
#' @keywords internal
get_stat_label <- function(stat_type) {
    switch(stat_type,
           "mean_sd" = "Mean +/- SD",
           "median_iqr" = "Median [IQR]",
           "median_range" = "Median (Range)",
           "range" = "Range",
           "n_miss" = "Missing",
           stat_type
           )
}

#' Format categorical statistic
#' @keywords internal
format_categorical_stat <- function(n, total, stat_type) {
    if (length(stat_type) > 1) stat_type <- stat_type[1]
    
    ## Format n with commas if >= 1000
    n_formatted <- if (n >= 1000) base::format(n, big.mark = ",") else as.character(n)
    
    switch(stat_type,
           "n" = n_formatted,
           "percent" = sprintf("%.1f%%", 100 * n / total),
           "n_percent" = sprintf("%s (%.1f%%)", n_formatted, 100 * n / total),
           n_formatted
           )
}

#' Perform continuous variable test
#' @keywords internal
perform_continuous_test <- function(data, var, group_var, test_type, stat_type = NULL, ...) {
    groups <- unique(data[[group_var]])
    groups <- groups[!is.na(groups)]
    n_groups <- length(groups)
    
    if (n_groups < 2) return(NA_real_)
    
    ## Remove NA values
    clean_data <- data[!is.na(get(var)) & !is.na(get(group_var))]
    
    ## Auto-select test based on statistic type and number of groups
    if (test_type == "auto") {
        if (!is.null(stat_type)) {
            if (grepl("mean", stat_type)) {
                ## For means, use parametric tests
                test_type <- if (n_groups == 2) "t" else "aov"
            } else if (grepl("median", stat_type)) {
                ## For medians, use non-parametric tests  
                test_type <- if (n_groups == 2) "wrs" else "kwt"
            } else {
                ## Default to non-parametric
                test_type <- if (n_groups == 2) "wrs" else "kwt"
            }
        } else {
            ## Original default behavior
            test_type <- if (n_groups == 2) "wrs" else "kwt"
        }
    }
    
    ## Perform test
    formula <- stats::as.formula(paste(var, "~", group_var))
    
    result <- tryCatch({
        switch(test_type,
               "t" = stats::t.test(formula, data = clean_data, ...)$p.value,
               "wrs" = stats::wilcox.test(formula, data = clean_data, ...)$p.value,
               "aov" = summary(stats::aov(formula, data = clean_data))[[1]][["Pr(>F)"]][1],
               "kwt" = stats::kruskal.test(formula, data = clean_data, ...)$p.value,
               NA_real_
               )
    }, error = function(e) NA_real_)
    
    return(result)
}

#' Perform categorical variable test
#' @keywords internal
perform_categorical_test <- function(data, var, group_var, test_type, ...) {
    
    ## Ensure we only use complete cases for both variables
    valid_idx <- !is.na(data[[var]]) & !is.na(data[[group_var]])
    
    ## Create contingency table with matched vectors
    tab <- table(data[[var]][valid_idx], data[[group_var]][valid_idx], useNA = "no")
    
    if (any(dim(tab) < 2)) return(NA_real_)
    
    ## Auto-select test
    if (test_type == "auto") {
        ## Suppress warning only for the expected value calculation
        expected <- suppressWarnings(stats::chisq.test(tab, simulate.p.value = FALSE)$expected)
        if (any(expected < 5)) {
            test_type <- "fisher.test"
        } else {
            test_type <- "chisq.test"
        }
    }
    
    ## Perform test
    result <- tryCatch({
        switch(test_type,
               "fisher.test" = stats::fisher.test(tab, ...)$p.value,
               "chisq.test" = stats::chisq.test(tab, ...)$p.value,
               NA_real_
               )
    }, error = function(e) NA_real_)
    
    return(result)
}

#' Format p-values in result table
#' @keywords internal
format_pvalues_desctbl <- function(result, digits_p) {
    if ("p_value" %chin% names(result)) {
        ## Handle character p_values first, then numeric
        result[, p_formatted := character(.N)]
        
        ## Use vectorized operations with proper conditionals
        result[!is.na(Variable) & Variable == "N", p_formatted := ""]
        result[is.na(p_value) & (is.na(Variable) | Variable != "N"), p_formatted := ""]
        
        ## Handle character p_values (already formatted)
        if (is.character(result$p_value)) {
            result[!is.na(p_value) & !is.na(Variable) & Variable != "N", 
                   p_formatted := p_value]
        } else {
            ## Handle numeric p_values
            result[!is.na(p_value) & p_value < 0.001 & (is.na(Variable) | Variable != "N"), 
                   p_formatted := "< 0.001"]
            result[!is.na(p_value) & p_value >= 0.001 & (is.na(Variable) | Variable != "N"), 
                   p_formatted := sprintf(paste0("%.", digits_p, "f"), p_value)]
        }
        
        ## Move p-value to last column
        cols <- names(result)
        cols <- c(cols[!cols %chin% c("p_value", "p_formatted")], "p_formatted")
        result <- result[, ..cols]
        setnames(result, "p_formatted", "p-value")
    }
    return(result)
}

#' Reorder columns to position total column
#' @keywords internal
reorder_total_column <- function(result, total, total_label) {
    if (total_label %chin% names(result)) {
        cols <- names(result)
        
        ## Get base columns and group columns
        base_cols <- c("Variable", "Group")
        group_cols <- cols[!cols %chin% c(base_cols, total_label, "p-value")]
        p_col <- if ("p-value" %chin% cols) "p-value" else NULL
        
        ## total = TRUE puts it first, total = "last" puts it last
        if (isTRUE(total) || (is.character(total) && total == "first")) {
            new_order <- c(base_cols, total_label, group_cols, p_col)
        } else if (is.character(total) && total == "last") {
            new_order <- c(base_cols, group_cols, total_label, p_col)
        } else {
            ## Default to first position for any other truthy value
            new_order <- c(base_cols, total_label, group_cols, p_col)
        }
        
        ## Keep only columns that exist
        new_order <- new_order[!is.na(new_order) & new_order %chin% cols]
        result <- result[, ..new_order]
    }
    return(result)
}
