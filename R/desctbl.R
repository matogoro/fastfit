#' Create Descriptive Statistics Table with Group Comparisons
#'
#' Generates publication-ready descriptive statistics tables with automatic 
#' variable type detection, group comparisons, and statistical testing. 
#' Handles continuous, categorical, and survival variables with appropriate
#' summary statistics and tests.
#'
#' @param data A data.frame or data.table containing the dataset to summarize.
#' @param by Character string specifying the grouping variable for stratified 
#'   analysis. When NULL (default), produces overall summaries only.
#' @param variables Character vector of variable names to summarize. Supports
#'   standard column names and Surv() expressions (e.g., "Surv(time, status)").
#' @param stats_continuous Character vector specifying statistics to compute for
#'   continuous variables. Options include "mean_sd", "median_iqr", 
#'   "median_range", and "range". Default is c("mean_sd", "median_iqr", "range").
#' @param stats_categorical Character string specifying format for categorical
#'   summaries: "n" (count only), "percent" (percentage only), or "n_percent"
#'   (count with percentage). Default is "n_percent".
#' @param digits Integer specifying decimal places for continuous statistics.
#'   Default is 1.
#' @param digits_p Integer specifying decimal places for p-values. Values less
#'   than 10^(-digits_p) display as "< 0.001" etc. Default is 3.
#' @param na_include Logical. If TRUE, missing values are displayed as a 
#'   separate category. Default is FALSE.
#' @param na_label Character string used to label missing values when 
#'   na_include = TRUE. Default is "Unknown".
#' @param na_percent Logical. Controls percentage calculation when na_include = TRUE.
#'   Categorical values only. If TRUE, percentages include missing/NAs in denominator
#'   (all sum to 100%). If FALSE, percentages exclude missing/NAs from denominator.
#'   Default is FALSE.
#' @param test Logical. If TRUE, performs appropriate statistical tests for
#'   group comparisons. Default is TRUE.
#' @param test_continuous Character string specifying test for continuous variables:
#'   "auto" (automatic selection based on statistic type), "t" (t-test),
#'   "aov" (ANOVA), "wrs" (Wilcoxon rank-sum), or "kwt" (Kruskal-Wallis).
#'   Default is "auto".
#' @param test_categorical Character string specifying test for categorical variables:
#'   "auto" (automatic selection based on expected frequencies), "fisher"
#'   (Fisher's exact test), or "chisq" (chi-squared test). Default is "auto".
#' @param total Logical or character string. Controls total column placement:
#'   TRUE or "first" places total column first, "last" places it after groups,
#'   FALSE excludes it. Default is TRUE.
#' @param total_label Character string for the total column header. 
#'   Default is "Total".
#' @param var_labels Named character vector for custom variable labels. Names
#'   should match variable names, values are display labels.
#' @param ... Additional arguments passed to statistical test functions.
#'
#' @return A data.table with class "desctbl" containing formatted descriptive
#'   statistics. The table has columns for variable, group, and statistics by
#'   group. Numeric values >= 1000 are formatted with commas.
#'   
#'   The returned object includes attributes:
#'   - raw_data: Numeric version of the table for further analysis
#'   - by_variable: The grouping variable used
#'   - variables: Vector of analyzed variables
#'
#' @details
#' The function automatically detects variable types and applies appropriate
#' summaries. For continuous variables with test = "auto", t-tests or ANOVA
#' are used for means, while Wilcoxon or Kruskal-Wallis tests are used for
#' medians. For categorical variables with test = "auto", Fisher's exact test
#' is used when expected frequencies are less than 5, otherwise chi-squared.
#' 
#' Survival variables should be specified using Surv() syntax from the survival
#' package. These display median survival with 95\% confidence intervals and
#' use log-rank tests for group comparisons.
#'
#' @examples
#' \dontrun{
#' # Basic descriptive table
#' desctbl(mtcars, variables = c("mpg", "cyl", "am"))
#' 
#' # Grouped comparison with custom labels
#' labels <- c(mpg = "Miles per Gallon", cyl = "Cylinders", am = "Transmission")
#' desctbl(mtcars, by = "vs", variables = c("mpg", "cyl", "am"),
#'         var_labels = labels)
#' 
#' # Including missing values
#' desctbl(data, by = "treatment", variables = c("age", "sex"),
#'         na_include = TRUE, na_label = "Missing")
#' 
#' # Survival analysis
#' desctbl(lung, by = "sex", 
#'         variables = c("age", "ph.ecog", "Surv(time, status)"))
#' 
#' # Access raw numeric data
#' result <- desctbl(data, by = "group", variables = vars)
#' raw_data <- attr(result, "raw_data")
#' }
#'
#' @export
#' @seealso \code{\link{exporttbl}} for exporting to PDF/LaTeX/HTML,
#'   \code{\link{fastfit}} for regression analysis
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
    
    ## Initialize both result tables
    result <- data.table::data.table()
    raw_result <- data.table::data.table()
    
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
        
        result <- rbind(result, var_data$formatted, fill = TRUE)
        raw_result <- rbind(raw_result, var_data$raw, fill = TRUE)
    }

    ## Standardize column names
    if ("variable" %in% names(result)) {
        data.table::setnames(result, "variable", "Variable")
        data.table::setnames(raw_result, "variable", "Variable")
    }
    if ("level" %in% names(result)) {
        data.table::setnames(result, "level", "Group")
        data.table::setnames(raw_result, "level", "Group")
    }

    ## Add p-value column if tests requested (after standardization)
    if (test && !is.null(group_var)) {
        result <- format_pvalues_desctbl(result, digits_p)
    }

    ## Add N row as first row if grouped (after all other processing)
    if (!is.null(group_var)) {
        
        ## Get the group values in the correct order
        if (is.factor(data[[group_var]])) {
            ## Use factor levels for proper ordering
            groups <- levels(data[[group_var]])
        } else {
            ## Fall back to unique values for non-factors
            groups <- unique(data[[group_var]])
            groups <- groups[!is.na(groups)]
        }
        
        ## Create N row
        n_row <- data.table::data.table(
                                 Variable = "N",
                                 Group = ""
                             )
        
        ## Calculate and add total if present
        if (total_label %in% names(result)) {
            n_total <- nrow(data)
            n_row[[total_label]] <- format(n_total, big.mark = ",")
        }
        
        ## Calculate for each group in the correct order
        for (grp in groups) {
            grp_col <- as.character(grp)
            if (grp_col %in% names(result)) {
                n_group <- sum(data[[group_var]] == grp, na.rm = TRUE)
                n_row[[grp_col]] <- format(n_group, big.mark = ",")
            }
        }
        
        ## Add empty p-value column if it exists
        if ("p-value" %in% names(result)) {
            n_row[["p-value"]] <- ""
        }
        
        ## Prepend N row
        result <- rbind(n_row, result, fill = TRUE)
    }

    ## Reorder columns if total position specified
    if (!isFALSE(total) && !is.null(group_var)) {
        result <- reorder_total_column(result, total, total_label)
    }

    ## Attach raw data and metadata as attributes
    data.table::setattr(result, "raw_data", raw_result)
    data.table::setattr(result, "by_variable", group_var)
    data.table::setattr(result, "variables", variables)

    result[]
    return(result)
}
