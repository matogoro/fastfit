#' Process variable wrapper - returns both formatted and raw
#' @keywords internal
process_variable <- function(data, var, group_var = NULL, 
                             stats_continuous, stats_categorical,
                             digits, na_include, na_label,
                             test, test_continuous, test_categorical,
                             total, total_label, var_labels, na_percent, ...) {
    
    ## Get variable label
    var_label <- if (!is.null(var_labels) && var %in% names(var_labels)) {
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
    
    formatted_result <- data.table::data.table()
    raw_result <- data.table::data.table()
    
    if (!is.null(group_var)) {
        groups <- unique(data[[group_var]])
        groups <- groups[!is.na(groups)]
        
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
            
            formatted_result <- rbind(formatted_result, data.table::as.data.table(formatted_row), fill = TRUE)
            raw_result <- rbind(raw_result, data.table::as.data.table(raw_row), fill = TRUE)
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
                miss_raw[[paste0(total_label, "_n")]] <- nrow(data)
            }
            
            ## Add group columns
            for (grp in groups) {
                grp_data <- data[get(group_var) == grp]
                grp_col <- as.character(grp)
                n_miss <- sum(is.na(grp_data[[var]]))
                
                ## Just show count, not percentage for continuous
                miss_formatted[[grp_col]] <- if (n_miss >= 1000) {
                                                 format(n_miss, big.mark = ",")
                                             } else {
                                                 as.character(n_miss)
                                             }
                miss_raw[[grp_col]] <- n_miss
                miss_raw[[paste0(grp_col, "_n")]] <- nrow(grp_data)
            }
            
            formatted_result <- rbind(formatted_result, 
                                      data.table::as.data.table(miss_formatted), 
                                      fill = TRUE)
            raw_result <- rbind(raw_result, 
                                data.table::as.data.table(miss_raw), 
                                fill = TRUE)
        }
    } else {
        
        ## No grouping variable case
        first_stat <- TRUE
        
        for (stat_type in stats) {
            formatted_row <- list(
                variable = if (first_stat) var_label else "",
                level = get_stat_label(stat_type)
            )
            formatted_row[[total_label]] <- calc_continuous_stat(data[[var]], stat_type, digits)
            
            raw_row <- list(
                variable = if (first_stat) var_label else "",
                level = stat_type,
                stat_type = stat_type
            )
            raw_row[[total_label]] <- NA_real_  # Would need appropriate calculation here
            
            formatted_result <- rbind(formatted_result, 
                                      data.table::as.data.table(formatted_row), 
                                      fill = TRUE)
            raw_result <- rbind(raw_result, 
                                data.table::as.data.table(raw_row), 
                                fill = TRUE)
            first_stat <- FALSE
        }
        
        ## Missing for ungrouped case
        if (na_include && any(is.na(data[[var]]))) {
            n_miss <- sum(is.na(data[[var]]))
            miss_formatted <- list(
                variable = "",
                level = na_label
            )
            miss_formatted[[total_label]] <- if (n_miss >= 1000) format(n_miss, big.mark = ",") else as.character(n_miss)
            
            miss_raw <- list(
                variable = "",
                level = na_label,
                stat_type = "missing"
            )
            miss_raw[[total_label]] <- n_miss
            
            formatted_result <- rbind(formatted_result, 
                                      data.table::as.data.table(miss_formatted), 
                                      fill = TRUE)
            raw_result <- rbind(raw_result, 
                                data.table::as.data.table(miss_raw), 
                                fill = TRUE)
        }
    }
    
    return(list(formatted = formatted_result, raw = raw_result))
}

#' Process categorical variable
#' @keywords internal
process_categorical <- function(data, var, var_label, group_var, stats,
                                na_include, na_label, test, test_type,
                                total, total_label, na_percent, ...) {
    
    formatted_result <- data.table::data.table()
    raw_result <- data.table::data.table()
    
    ## Get levels
    if (is.factor(data[[var]])) {
        levels <- levels(data[[var]])
    } else {
        levels <- sort(unique(data[[var]]))
        levels <- levels[!is.na(levels)]
    }
    
    if (na_include && any(is.na(data[[var]]))) {
        levels <- c(levels, na_label)
    }
    
    ## Calculate test p-value
    p_val <- NA_real_
    if (test && !is.null(group_var)) {
        p_val <- perform_categorical_test(data, var, group_var, test_type, ...)
    }
    
    first_row <- TRUE
    
    for (i in seq_along(levels)) {
        level <- levels[i]
        
        formatted_row <- list(
            variable = if (first_row) var_label else "",
            level = as.character(level)
        )
        
        raw_row <- list(
            variable = if (first_row) var_label else "",
            level = as.character(level),
            stat_type = "categorical"
        )
        
        if (!is.null(group_var)) {
            groups <- unique(data[[group_var]])
            groups <- groups[!is.na(groups)]
            
            ## Add total column
            if (!isFALSE(total)) {
                if (level == na_label) {
                    n <- sum(is.na(data[[var]]))
                    total_n <- nrow(data)
                    
                    ## Only show percentage if na_percent = TRUE
                    if (na_percent) {
                        formatted_row[[total_label]] <- format_categorical_stat(n, total_n, stats)
                    } else {
                        
                        ## Just show count for NA when na_percent = FALSE
                        formatted_row[[total_label]] <- if (n >= 1000) {
                                                            format(n, big.mark = ",")
                                                        } else {
                                                            as.character(n)
                                                        }
                    }
                } else {
                    n <- sum(data[[var]] == level, na.rm = TRUE)
                    
                    ## Denominator depends on na_percent
                    total_n <- if (na_percent) nrow(data) else sum(!is.na(data[[var]]))
                    formatted_row[[total_label]] <- format_categorical_stat(n, total_n, stats)
                }
                
                ## Raw columns remain the same
                raw_row[[total_label]] <- n
                raw_row[[paste0(total_label, "_pct")]] <- 100 * n / total_n
                raw_row[[paste0(total_label, "_total")]] <- total_n
            }
            
            ## Add group columns
            for (grp in groups) {
                grp_data <- data[get(group_var) == grp]
                grp_col <- as.character(grp)
                
                if (level == na_label) {
                    n <- sum(is.na(grp_data[[var]]))
                    total_n <- nrow(grp_data)
                    
                    ## Only show percentage if na_percent = TRUE
                    if (na_percent) {
                        formatted_row[[grp_col]] <- format_categorical_stat(n, total_n, stats)
                    } else {
                        
                        ## Just show count for NA when na_percent = FALSE
                        formatted_row[[grp_col]] <- if (n >= 1000) {
                                                        format(n, big.mark = ",")
                                                    } else {
                                                        as.character(n)
                                                    }
                    }
                } else {
                    n <- sum(grp_data[[var]] == level, na.rm = TRUE)
                    
                    ## Denominator depends on na_percent
                    total_n <- if (na_percent) nrow(grp_data) else sum(!is.na(grp_data[[var]]))
                    formatted_row[[grp_col]] <- format_categorical_stat(n, total_n, stats)
                }
                
                ## Raw columns
                raw_row[[grp_col]] <- n
                raw_row[[paste0(grp_col, "_pct")]] <- 100 * n / total_n
                raw_row[[paste0(grp_col, "_total")]] <- total_n
            }
            
            ## Add p_value to first row only
            if (first_row && test && !is.na(p_val)) {
                formatted_row[["p_value"]] <- p_val
                raw_row[["p_value"]] <- p_val
            }
        } else {
            
            ## No grouping case - apply same logic
            if (level == na_label) {
                n <- sum(is.na(data[[var]]))
                total_n <- nrow(data)
                
                if (na_percent) {
                    formatted_row[[total_label]] <- format_categorical_stat(n, total_n, stats)
                } else {
                    formatted_row[[total_label]] <- if (n >= 1000) {
                                                        format(n, big.mark = ",")
                                                    } else {
                                                        as.character(n)
                                                    }
                }
            } else {
                n <- sum(data[[var]] == level, na.rm = TRUE)
                total_n <- if (na_percent) nrow(data) else sum(!is.na(data[[var]]))
                formatted_row[[total_label]] <- format_categorical_stat(n, total_n, stats)
            }
            
            raw_row[[total_label]] <- n
        }
        
        formatted_result <- rbind(formatted_result, data.table::as.data.table(formatted_row), fill = TRUE)
        raw_result <- rbind(raw_result, data.table::as.data.table(raw_row), fill = TRUE)
        first_row <- FALSE
    }
    
    return(list(formatted = formatted_result, raw = raw_result))
}

#' Process survival variable
#' @keywords internal
process_survival <- function(data, var, var_label, group_var, digits,
                             test, total, total_label, na_include = FALSE, 
                             na_label = "Unknown", ...) {
    
    if (!requireNamespace("survival", quietly = TRUE)) {
        stop("Package 'survival' required for survival statistics")
    }
    
    formatted_result <- data.table::data.table()
    raw_result <- data.table::data.table()
    
    ## Parse the Surv() expression
    surv_expr <- gsub("Surv\\(|\\)", "", var)
    surv_parts <- trimws(strsplit(surv_expr, ",")[[1]])
    time_var <- surv_parts[1]
    event_var <- surv_parts[2]
    
    ## Build rows with same structure
    formatted_row <- list(
        variable = var_label,
        level = "Median (95% CI)"
    )
    
    raw_row <- list(
        variable = var_label,
        level = "survival_median",
        stat_type = "survival"
    )
    
    if (!is.null(group_var)) {
        groups <- unique(data[[group_var]])
        groups <- groups[!is.na(groups)]
        
        ## Calculate overall median if requested
        if (!isFALSE(total)) {
            overall_fit <- survival::survfit(
                                         survival::Surv(data[[time_var]], data[[event_var]]) ~ 1
                                     )
            overall_median <- summary(overall_fit)$table["median"]
            overall_lower <- summary(overall_fit)$table["0.95LCL"]
            overall_upper <- summary(overall_fit)$table["0.95UCL"]
            
            ## Formatted
            formatted_row[[total_label]] <- format_survival_stat(
                overall_median, overall_lower, overall_upper, digits
            )
            
            ## Raw with CI columns
            raw_row[[total_label]] <- overall_median
            raw_row[[paste0(total_label, "_lower95")]] <- overall_lower
            raw_row[[paste0(total_label, "_upper95")]] <- overall_upper
            raw_row[[paste0(total_label, "_n")]] <- overall_fit$n
        }
        
        ## Calculate for each group
        for (grp in groups) {
            grp_data <- data[get(group_var) == grp]
            grp_col <- as.character(grp)
            grp_fit <- survival::survfit(
                                     survival::Surv(grp_data[[time_var]], grp_data[[event_var]]) ~ 1
                                 )
            grp_median <- summary(grp_fit)$table["median"]
            grp_lower <- summary(grp_fit)$table["0.95LCL"]
            grp_upper <- summary(grp_fit)$table["0.95UCL"]
            
            ## Formatted
            formatted_row[[grp_col]] <- format_survival_stat(
                grp_median, grp_lower, grp_upper, digits
            )
            
            ## Raw with CI columns
            raw_row[[grp_col]] <- grp_median
            raw_row[[paste0(grp_col, "_lower95")]] <- grp_lower
            raw_row[[paste0(grp_col, "_upper95")]] <- grp_upper
            raw_row[[paste0(grp_col, "_n")]] <- grp_fit$n
        }
        
        ## Add p_value if test is TRUE
        if (test) {
            formula <- as.formula(paste0("survival::Surv(", time_var, ", ", 
                                         event_var, ") ~ ", group_var))
            logrank <- survival::survdiff(formula, data = data)
            p_val <- 1 - stats::pchisq(logrank$chisq, length(logrank$n) - 1)
            
            formatted_row[["p_value"]] <- p_val
            raw_row[["p_value"]] <- p_val
        }
    } else {
        
        ## No grouping variable case
        overall_fit <- survival::survfit(
                                     survival::Surv(data[[time_var]], data[[event_var]]) ~ 1
                                 )
        overall_median <- summary(overall_fit)$table["median"]
        overall_lower <- summary(overall_fit)$table["0.95LCL"]
        overall_upper <- summary(overall_fit)$table["0.95UCL"]
        
        formatted_row[[total_label]] <- format_survival_stat(
            overall_median, overall_lower, overall_upper, digits
        )
        raw_row[[total_label]] <- overall_median
    }
    
    formatted_result <- rbind(formatted_result, data.table::as.data.table(formatted_row), fill = TRUE)
    raw_result <- rbind(raw_result, data.table::as.data.table(raw_row), fill = TRUE)
    
    ## Add missing row if requested
    if (na_include) {
        
        ## Check for missing in either time or event variable
        n_miss_time <- sum(is.na(data[[time_var]]))
        n_miss_event <- sum(is.na(data[[event_var]]))
        n_miss <- max(n_miss_time, n_miss_event)  # Use the maximum
        
        if (n_miss > 0) {
            miss_formatted <- list(
                variable = "",
                level = na_label
            )
            
            miss_raw <- list(
                variable = "",
                level = na_label,
                stat_type = "missing"
            )
            
            if (!is.null(group_var)) {
                
                ## Add total column
                if (!isFALSE(total)) {
                    miss_formatted[[total_label]] <- if (n_miss >= 1000) {
                                                         format(n_miss, big.mark = ",")
                                                     } else {
                                                         as.character(n_miss)
                                                     }
                    miss_raw[[total_label]] <- n_miss
                    miss_raw[[paste0(total_label, "_n")]] <- nrow(data)
                }
                
                ## Add group columns
                for (grp in groups) {
                    grp_data <- data[get(group_var) == grp]
                    grp_miss_time <- sum(is.na(grp_data[[time_var]]))
                    grp_miss_event <- sum(is.na(grp_data[[event_var]]))
                    grp_miss <- max(grp_miss_time, grp_miss_event)
                    grp_col <- as.character(grp)
                    
                    miss_formatted[[grp_col]] <- if (grp_miss >= 1000) {
                                                     format(grp_miss, big.mark = ",")
                                                 } else {
                                                     as.character(grp_miss)
                                                 }
                    miss_raw[[grp_col]] <- grp_miss
                    miss_raw[[paste0(grp_col, "_n")]] <- nrow(grp_data)
                }
            } else {
                
                ## No grouping case
                miss_formatted[[total_label]] <- if (n_miss >= 1000) {
                                                     format(n_miss, big.mark = ",")
                                                 } else {
                                                     as.character(n_miss)
                                                 }
                miss_raw[[total_label]] <- n_miss
            }
            
            formatted_result <- rbind(formatted_result, 
                                      data.table::as.data.table(miss_formatted), 
                                      fill = TRUE)
            raw_result <- rbind(raw_result, 
                                data.table::as.data.table(miss_raw), 
                                fill = TRUE)
        }
    }
    
    return(list(formatted = formatted_result, raw = raw_result))
}

#' Format survival statistic
#' @keywords internal
format_survival_stat <- function(median, lower, upper, digits) {
    if (is.na(median)) {
        return("Not reached")
    }
    sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
            median, lower, upper)
}

#' Calculate continuous statistic
#' @keywords internal
calc_continuous_stat <- function(x, stat_type, digits) {
    x <- x[!is.na(x)]
    
    if (length(x) == 0) return("—")
    
    switch(stat_type,
           "mean_sd" = {
               mean_val <- mean(x)
               sd_val <- sd(x)
               
               ## Format with commas if >= 1000
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
               paste0(mean_str, " ± ", sd_str)
           },
           "median_iqr" = {
               med <- median(x)
               q1 <- quantile(x, 0.25)
               q3 <- quantile(x, 0.75)
               
               ## Format with commas if >= 1000
               med_str <- if (abs(med) >= 1000) {
                              format(round(med, digits), big.mark = ",")
                          } else {
                              sprintf(paste0("%.", digits, "f"), med)
                          }
               q1_str <- if (abs(q1) >= 1000) {
                             format(round(q1, digits), big.mark = ",")
                         } else {
                             sprintf(paste0("%.", digits, "f"), q1)
                         }
               q3_str <- if (abs(q3) >= 1000) {
                             format(round(q3, digits), big.mark = ",")
                         } else {
                             sprintf(paste0("%.", digits, "f"), q3)
                         }
               ## Use comma separator for IQR with negative numbers
               if (q1 < 0 || q3 < 0) {
                   paste0(med_str, " [", q1_str, ", ", q3_str, "]")
               } else {
                   paste0(med_str, " [", q1_str, "-", q3_str, "]")
               }
           },
           "median_range" = {
               med <- median(x)
               min_val <- min(x)
               max_val <- max(x)
               med_str <- if (abs(med) >= 1000) {
                              format(round(med, digits), big.mark = ",")
                          } else {
                              sprintf(paste0("%.", digits, "f"), med)
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
           "mean_sd" = "Mean ± SD",
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
    
    ## Create contingency table
    tab <- table(data[[var]], data[[group_var]], useNA = "no")
    
    if (any(dim(tab) < 2)) return(NA_real_)
    
    ## Auto-select test
    if (test_type == "auto") {
        expected <- stats::chisq.test(tab, simulate.p.value = FALSE)$expected
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
    if ("p_value" %in% names(result)) {
                                        # Process row by row to avoid vectorization issues
        result[, p_formatted := character(.N)]  # Initialize as character vector
        
        for (i in seq_len(nrow(result))) {
            if (!is.na(result$Variable[i]) && result$Variable[i] == "N") {
                result$p_formatted[i] <- ""
            } else if (is.na(result$p_value[i])) {
                result$p_formatted[i] <- ""
            } else if (is.character(result$p_value[i])) {
                result$p_formatted[i] <- result$p_value[i]
            } else if (result$p_value[i] < 0.001) {
                result$p_formatted[i] <- "< 0.001"
            } else {
                result$p_formatted[i] <- sprintf(paste0("%.", digits_p, "f"), result$p_value[i])
            }
        }
        
                                        # Move p-value to last column
        cols <- names(result)
        cols <- c(cols[!cols %in% c("p_value", "p_formatted")], "p_formatted")
        result <- result[, ..cols]
        data.table::setnames(result, "p_formatted", "p-value")
    }
    return(result)
}

#' Reorder columns to position total column
#' @keywords internal
reorder_total_column <- function(result, total, total_label) {
    if (total_label %in% names(result)) {
        cols <- names(result)
        
        ## Get base columns and group columns
        base_cols <- c("variable", "level")
        group_cols <- cols[!cols %in% c(base_cols, total_label, "p-value")]
        p_col <- if ("p-value" %in% cols) "p-value" else NULL
        
        ## Modified logic: total = TRUE puts it first, total = "last" puts it last
        if (isTRUE(total) || (is.character(total) && total == "first")) {
            new_order <- c(base_cols, total_label, group_cols, p_col)
        } else if (is.character(total) && total == "last") {
            new_order <- c(base_cols, group_cols, total_label, p_col)
        } else {
            ## Default to first position for any other truthy value
            new_order <- c(base_cols, total_label, group_cols, p_col)
        }
        
        ## Keep only columns that exist
        new_order <- new_order[!is.na(new_order) & new_order %in% cols]
        result <- result[, ..new_order]
    }
    return(result)
}
