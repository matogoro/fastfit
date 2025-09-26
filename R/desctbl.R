#' Create Descriptive Statistics Table
#'
#' Generates descriptive statistics tables with group comparisons and statistical tests.
#'
#' @param data A data.frame or data.table.
#' @param by Character string for grouping variable. Default NULL.
#' @param variables Character vector of variable names to summarize.
#' @param stats_continuous Character vector for continuous stats. Default c("mean_sd", "median_iqr", "range").
#' @param stats_categorical Character string for categorical format. Default "n_percent".
#' @param digits Integer for decimal places. Default 1.
#' @param digits_p Integer for p-value decimal places. Default 3.
#' @param na_include Logical to include missing values. Default FALSE.
#' @param na_label Character string for missing label. Default "Unknown".
#' @param na_percent Logical for missing percentages. Default FALSE.
#' @param test Logical to perform tests. Default TRUE.
#' @param test_continuous Character for continuous tests. Default "auto".
#' @param test_categorical Character for categorical tests. Default "auto".
#' @param total Logical or character for total column. Default TRUE.
#' @param total_label Character string for total label. Default "Total".
#' @param var_labels Named character vector for labels.
#' @param ... Additional arguments.
#'
#' @return A data.table with descriptive statistics.
#'
#' @export
desctbl <- function(data,
                    by = NULL,
                    variables,
                    stats_continuous = c("mean_sd", "median_iqr", "range"),
                    stats_categorical = "n_percent",
                    digits = 1,
                    digits_p = 3,
                    na_include = FALSE,
                    na_label = "Unknown",
                    na_percent = FALSE,
                    test = TRUE,
                    test_continuous = "auto",
                    test_categorical = "auto",
                    total = TRUE,
                    total_label = "Total",
                    var_labels = NULL,
                    ...) {
    
    if (!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
    ## Set group_var from 'by' parameter
    group_var <- by
    group_var_label <- NULL
    
    ## Apply label to group variable if provided
    if (!is.null(group_var) && !is.null(var_labels) && group_var %in% names(var_labels)) {
        group_var_label <- var_labels[group_var]
    } else if (!is.null(group_var)) {
        group_var_label <- group_var
    }
    
    ## Variables are already provided as a vector
    vars <- variables
    
    ## Initialize result table
    result <- data.table::data.table()
    
    ## Process each variable
    for (var in vars) {
        var_data <- process_variable(
            data = data,
            var = var,
            group_var = group_var,
            stats_continuous = stats_continuous,
            stats_categorical = stats_categorical,
            digits = digits,
            na_include = na_include,
            na_label = na_label,
            test = test,
            test_continuous = test_continuous,
            test_categorical = test_categorical,
            total = total,
            total_label = total_label,
            var_labels = var_labels,
            na_percent = na_percent,
            ...
        )
        
        result <- rbind(result, var_data, fill = TRUE)
    }
    
    ## Add p-value column if tests requested
    if (test && !is.null(group_var)) {
        result <- format_p_values(result, digits_p)
    }
    
    ## Reorder columns if total position specified
    if (!isFALSE(total) && !is.null(group_var)) {
        result <- reorder_total_column(result, total, total_label)
    }
    
    result[]  # Force finalization
    return(result)
}

#' Process single variable for descriptive statistics
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
    
    ## Determine variable type
    if (grepl("^Surv\\(", var)) {
        ## Survival object
        result <- process_survival(
            data = data,
            var = var,
            var_label = var_label,
            group_var = group_var,
            digits = digits,
            test = test,
            total = total,
            total_label = total_label,
            ...
        )
    } else if (is.numeric(data[[var]])) {
        ## Continuous variable
        result <- process_continuous(
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
        )
    } else {
        ## Categorical variable
        result <- process_categorical(
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
        )
    }
    
    return(result)
}

#' Process continuous variable
#' @keywords internal
process_continuous <- function(data, var, var_label, group_var, stats, digits,
                               na_include, na_label, test, test_type,
                               total, total_label, ...) {
    
    result <- data.table::data.table()
    
    if (!is.null(group_var)) {
        groups <- unique(data[[group_var]])
        groups <- groups[!is.na(groups)]
        
        first_stat <- TRUE
        p_values <- list()
        
        ## Calculate all p-values first if testing (but not for range)
        if (test) {
            for (stat_type in stats) {
                
                ## Skip p-value calculation for raw range
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
            
            ## Build row as a list first
            row_list <- list(
                variable = if (first_stat) var_label else "",
                level = get_stat_label(stat_type)
            )
            
            ## Add total column if requested
            if (!isFALSE(total)) {
                total_stat <- calc_continuous_stat(data[[var]], stat_type, digits)
                row_list[[total_label]] <- total_stat
            }
            
            ## Add group columns
            for (grp in groups) {
                grp_data <- data[get(group_var) == grp, get(var)]
                grp_stat <- calc_continuous_stat(grp_data, stat_type, digits)
                row_list[[as.character(grp)]] <- grp_stat
            }
            
            ## Add p_value if test is TRUE and we have a value (exclude for raw range)
            if (test && stat_type != "range" && !is.null(p_values[[stat_type]])) {
                row_list[["p_value"]] <- p_values[[stat_type]]
            }
            
            ## Convert list to data.table
            row <- data.table::as.data.table(row_list)
            
            first_stat <- FALSE
            result <- rbind(result, row, fill = TRUE)
        }
        
    } else {
        
        ## No grouping variable
        first_stat <- TRUE
        for (stat_type in stats) {
            row_list <- list(
                variable = if (first_stat) var_label else "",
                level = get_stat_label(stat_type),
                value = calc_continuous_stat(data[[var]], stat_type, digits)
            )
            row <- data.table::as.data.table(row_list)
            first_stat <- FALSE
            result <- rbind(result, row, fill = TRUE)
        }
    }
    
    ## Missing count handling
    if (na_include) {
        n_miss <- sum(is.na(data[[var]]))
        if (n_miss > 0) {
            miss_list <- list(
                variable = "",
                level = na_label
            )
            
            if (!is.null(group_var)) {
                if (!isFALSE(total)) {
                    miss_list[[total_label]] <- paste0(n_miss, " (", 
                                                       sprintf("%.1f", 100 * n_miss / nrow(data)), "%)")
                }
                groups <- unique(data[[group_var]])
                groups <- groups[!is.na(groups)]
                for (grp in groups) {
                    grp_miss <- sum(is.na(data[get(group_var) == grp, get(var)]))
                    grp_n <- nrow(data[get(group_var) == grp])
                    miss_list[[as.character(grp)]] <- paste0(grp_miss, " (",
                                                             sprintf("%.1f", 100 * grp_miss / grp_n), "%)")
                }
            } else {
                miss_list[["value"]] <- paste0(n_miss, " (", 
                                               sprintf("%.1f", 100 * n_miss / nrow(data)), "%)")
            }
            
            miss_row <- data.table::as.data.table(miss_list)
            result <- rbind(result, miss_row, fill = TRUE)
        }
    }
    
    return(result)
}

#' Process categorical variable
#' @keywords internal
process_categorical <- function(data, var, var_label, group_var, stats,
                                na_include, na_label, test, test_type,
                                total, total_label, na_percent, ...) {
    
    result <- data.table::data.table()
    
    ## Get factor levels in original order
    if (is.factor(data[[var]])) {
        levels <- levels(data[[var]])
    } else {
        levels <- sort(unique(data[[var]]))
        levels <- levels[!is.na(levels)]
    }
    
    ## Add NA as a level if requested
    if (na_include && any(is.na(data[[var]]))) {
        levels <- c(levels, na_label)
    }
    
    ## Calculate test p-value first if requested (excluding missing values)
    p_val <- NA_real_
    if (test && !is.null(group_var)) {
        p_val <- perform_categorical_test(data, var, group_var, test_type, ...)
    }
    
    ## First row with variable name
    first_row <- TRUE
    
    for (i in seq_along(levels)) {
        level <- levels[i]
        
        ## Build row as a list first
        row_list <- list(
            variable = if (first_row) var_label else "",
            level = as.character(level)
        )
        
        if (!is.null(group_var)) {
            groups <- unique(data[[group_var]])
            groups <- groups[!is.na(groups)]
            
            ## Add total column if requested
            if (!isFALSE(total)) {
                if (level == na_label) {
                    ## For missing category
                    n <- sum(is.na(data[[var]]))
                    if (na_percent) {
                        ## Show n (%) with total N as denominator
                        total_n <- nrow(data)
                        row_list[[total_label]] <- format_categorical_stat(n, total_n, "n_percent")
                    } else {
                        ## Show only n
                        row_list[[total_label]] <- as.character(n)
                    }
                } else {
                    ## For non-missing categories
                    n <- sum(data[[var]] == level, na.rm = TRUE)
                    ## Denominator depends on na_percent setting
                    if (na_percent) {
                        ## Use total N so all percentages (including missing) sum to 100%
                        total_n <- nrow(data)
                    } else {
                        ## Use non-missing N so only valid categories sum to 100%
                        total_n <- sum(!is.na(data[[var]]))
                    }
                    row_list[[total_label]] <- format_categorical_stat(n, total_n, stats)
                }
            }
            
            ## Add group columns
            for (grp in groups) {
                grp_data <- data[get(group_var) == grp]
                if (level == na_label) {
                    ## For missing category
                    n <- sum(is.na(grp_data[[var]]))
                    if (na_percent) {
                        ## Show n (%) with group total as denominator
                        total_n <- nrow(grp_data)
                        row_list[[as.character(grp)]] <- format_categorical_stat(n, total_n, "n_percent")
                    } else {
                        ## Show only n
                        row_list[[as.character(grp)]] <- as.character(n)
                    }
                } else {
                    ## For non-missing categories
                    n <- sum(grp_data[[var]] == level, na.rm = TRUE)
                    ## Denominator depends on na_percent setting
                    if (na_percent) {
                        ## Use total group N so all percentages sum to 100%
                        total_n <- nrow(grp_data)
                    } else {
                        ## Use non-missing N so only valid categories sum to 100%
                        total_n <- sum(!is.na(grp_data[[var]]))
                    }
                    row_list[[as.character(grp)]] <- format_categorical_stat(n, total_n, stats)
                }
            }
            
            ## Add p_value to first row only if test is TRUE
            if (first_row && test && !is.na(p_val)) {
                row_list[["p_value"]] <- p_val
            }
        } else {
            ## No grouping
            if (level == na_label) {
                ## For missing category
                n <- sum(is.na(data[[var]]))
                if (na_percent) {
                    total_n <- nrow(data)
                    row_list[["value"]] <- format_categorical_stat(n, total_n, "n_percent")
                } else {
                    row_list[["value"]] <- as.character(n)
                }
            } else {
                ## For non-missing categories
                n <- sum(data[[var]] == level, na.rm = TRUE)
                ## Denominator depends on na_percent setting
                if (na_percent) {
                    ## Use total N so all percentages sum to 100%
                    total_n <- nrow(data)
                } else {
                    ## Use non-missing N so only valid categories sum to 100%
                    total_n <- sum(!is.na(data[[var]]))
                }
                row_list[["value"]] <- format_categorical_stat(n, total_n, stats)
            }
        }
        
        ## Convert list to data.table
        row <- data.table::as.data.table(row_list)
        
        first_row <- FALSE
        result <- rbind(result, row, fill = TRUE)
    }
    
    return(result)
}

#' Process survival variable
#' @keywords internal
process_survival <- function(data, var, var_label, group_var, digits,
                             test, total, total_label, ...) {
    
    if (!requireNamespace("survival", quietly = TRUE)) {
        stop("Package 'survival' required for survival statistics")
    }
    
    result <- data.table::data.table()
    
    ## Parse the Surv() expression
    surv_expr <- gsub("Surv\\(|\\)", "", var)
    surv_parts <- trimws(strsplit(surv_expr, ",")[[1]])
    time_var <- surv_parts[1]
    event_var <- surv_parts[2]
    
    if (!is.null(group_var)) {
        groups <- unique(data[[group_var]])
        groups <- groups[!is.na(groups)]
        
        ## Build row as a list
        row_list <- list(
            variable = var_label,
            level = "Median (95% CI)"
        )
        
        ## Calculate overall median if requested
        if (!isFALSE(total)) {
            overall_fit <- survival::survfit(
                                         survival::Surv(data[[time_var]], data[[event_var]]) ~ 1
                                     )
            overall_median <- summary(overall_fit)$table["median"]
            overall_lower <- summary(overall_fit)$table["0.95LCL"]
            overall_upper <- summary(overall_fit)$table["0.95UCL"]
            
            row_list[[total_label]] <- format_survival_stat(
                overall_median, overall_lower, overall_upper, digits
            )
        }
        
        ## Calculate for each group
        for (grp in groups) {
            grp_data <- data[get(group_var) == grp]
            grp_fit <- survival::survfit(
                                     survival::Surv(grp_data[[time_var]], grp_data[[event_var]]) ~ 1
                                 )
            grp_median <- summary(grp_fit)$table["median"]
            grp_lower <- summary(grp_fit)$table["0.95LCL"]
            grp_upper <- summary(grp_fit)$table["0.95UCL"]
            
            row_list[[as.character(grp)]] <- format_survival_stat(
                grp_median, grp_lower, grp_upper, digits
            )
        }
        
        ## Add p_value if test is TRUE
        if (test) {
            formula <- as.formula(paste0("survival::Surv(", time_var, ", ", 
                                         event_var, ") ~ ", group_var))
            logrank <- survival::survdiff(formula, data = data)
            p_val <- 1 - stats::pchisq(logrank$chisq, length(logrank$n) - 1)
            row_list[["p_value"]] <- p_val
        }
        
        ## Convert list to data.table
        row <- data.table::as.data.table(row_list)
        
    } else {
        ## No grouping variable
        overall_fit <- survival::survfit(
                                     survival::Surv(data[[time_var]], data[[event_var]]) ~ 1
                                 )
        overall_median <- summary(overall_fit)$table["median"]
        overall_lower <- summary(overall_fit)$table["0.95LCL"]
        overall_upper <- summary(overall_fit)$table["0.95UCL"]
        
        row_list <- list(
            variable = var_label,
            level = "Median survival (95% CI)",
            value = format_survival_stat(overall_median, overall_lower, 
                                         overall_upper, digits)
        )
        row <- data.table::as.data.table(row_list)
    }
    
    result <- rbind(result, row, fill = TRUE)
    return(result)
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
           "mean_sd" = sprintf(paste0("%.", digits, "f ± %.", digits, "f"), 
                               mean(x), sd(x)),
           "median_iqr" = sprintf(paste0("%.", digits, "f [%.", digits, "f-%.", digits, "f]"),
                                  median(x), quantile(x, 0.25), quantile(x, 0.75)),
           "median_range" = sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"),
                                    median(x), min(x), max(x)),
           "range" = sprintf(paste0("%.", digits, "f-%.", digits, "f"),
                             min(x), max(x)),
           "n_miss" = as.character(sum(is.na(x))),
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
    
    switch(stat_type,
           "n" = as.character(n),
           "percent" = sprintf("%.1f%%", 100 * n / total),
           "n_percent" = sprintf("%d (%.1f%%)", n, 100 * n / total),
           as.character(n)
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
format_p_values <- function(result, digits_p) {
    if ("p_value" %in% names(result)) {
        result[, p_formatted := ifelse(
                     is.na(p_value), "",
                                ifelse(p_value < 0.001, "< 0.001",
                                       sprintf(paste0("%.", digits_p, "f"), p_value))
                 )]
        
        ## Move p-value to last column
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
