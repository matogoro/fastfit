#' Create a Forest Plot for Cox Proportional Hazards Models
#'
#' Generates a publication-ready forest plot combining a data table and 
#' graphical representation of hazard ratios from a Cox proportional hazards model.
#' The plot includes variable names, group levels, sample sizes, hazard ratios 
#' with confidence intervals, and p-values.
#'
#' @param model A Cox proportional hazards model object of class \code{coxph} 
#'   from the \code{survival} package.
#' @param data A data frame or data.table containing the original data used to 
#'   fit the model. If NULL (default), the function attempts to extract data 
#'   from the model object.
#' @param title Character string specifying the plot title. Default is "Cox Proportional Hazards Model".
#' @param digits Integer specifying the number of decimal places for hazard ratios 
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
#'   allocated to the events column. Default is 0.06.
#' @param col_width_hr Numeric value specifying the proportion of plot width 
#'   allocated to the hazard ratio column. Default is 0.05.
#' @param ref_label Character string to display for reference categories. 
#'   Default is "reference".
#' @param var_labels Named character vector for custom variable labels. Names should 
#'   match variable names in the model, values are the display labels.
#' @param color Character string specifying the color for hazard ratio point estimates. 
#'   Default is "#8A61D8" (purple).
#'
#' @return A ggplot object containing the forest plot. The object has an attribute
#'   "recommended_dims" with suggested width and height in inches. If \code{ggpubr} 
#'   is installed, returns a ggplot object via \code{ggpubr::as_ggplot()}. 
#'   Otherwise, draws the plot directly and returns it invisibly.
#'
#' @details
#' The function creates a forest plot with the following features:
#' \itemize{
#'   \item Alternating row shading by variable groups
#'   \item Point sizes scaled by sample size
#'   \item Logarithmic scale for hazard ratios
#'   \item Model statistics displayed at the bottom (events, global p-value, 
#'     concordance index with 95% CI, AIC)
#'   \item Italic formatting for p-values in the display
#'   \item Automatic calculation of optimal dimensions
#'   \item Customizable column widths for the data table portion
#' }
#'
#' The function will suggest optimal plot dimensions based on the model complexity.
#' Access these via \code{attr(plot, "recommended_dims")} or note the message output.
#'
#' @examples
#' \dontrun{
#' library(survival)
#' library(ggplot2)
#' 
#' # Fit a Cox model
#' cox_model <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data = lung)
#' 
#' # Create forest plot and save with recommended dimensions
#' p <- coxforest(cox_model)
#' ggsave("forest.pdf", p, 
#'        width = attr(p, "recommended_dims")$width,
#'        height = attr(p, "recommended_dims")$height,
#'        units = "in")
#' 
#' # Create forest plot with specified dimensions and custom column widths
#' p <- coxforest(cox_model, 
#'                title = "Survival Analysis Results",
#'                plot_width = 12,
#'                plot_height = 8,
#'                col_width_var = 0.03,
#'                col_width_level = 0.20,
#'                var_labels = c("age" = "Age (years)",
#'                              "sex" = "Sex",
#'                              "ph.ecog" = "ECOG Performance Status"))
#' ggsave("forest.pdf", p, width = 12, height = 8, units = "in")
#' }
#'
#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom stats AIC confint
#' @importFrom grDevices axisTicks
#'
#' @seealso 
#' \code{\link[survival]{coxph}} for fitting Cox models,
#' \code{\link[ggplot2]{ggplot}} for the underlying plotting system
#'
#' @author Paul H. McClelland <PaulHMcClelland@protonmail.com>
#'
coxforest <- function(model, data = NULL,
                      title = "Cox Proportional Hazards Model",
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
                      col_width_hr = 0.05,
                      show_n_events = c("n", "Events"),
                      ref_label = "reference",
                      var_labels = NULL,
                      color = "#8A61D8") {

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
    
    stopifnot(inherits(model, "coxph"))
    
    ## Get data
    if(is.null(data)){
        warning("The `data` argument is not provided. Data will be extracted from model fit.")
        
        ## First try to get data from model$data (stored by mmodel)
        if (!is.null(model$data)) {
            data <- model$data
        } else {
            ## Try to evaluate the call
            tryCatch({
                data <- eval(model$call$data)
            }, error = function(e) {
                data <- NULL
            })
            
            ## Try model$model as fallback
            if (is.null(data))
                data <- model$model
        }
        
        if (is.null(data))
            stop("The `data` argument should be provided either to coxforest or stored in the model.")
    }
    
    ## Convert to data.table if not already
    if(!data.table::is.data.table(data)) {
        data <- data.table::as.data.table(data)
    }
    
    terms <- attr(model$terms, "dataClasses")[-1]

    ## Extract coefficients
    coef_summary <- summary(model)$coefficients
    conf_int <- stats::confint(model)
    
    coef <- data.table::data.table(
                            term = rownames(coef_summary),
                            estimate = coef_summary[, "coef"],
                            std_error = coef_summary[, "se(coef)"],
                            statistic = coef_summary[, "z"],
                            p_value = coef_summary[, "Pr(>|z|)"],
                            conf_low = conf_int[, 1],
                            conf_high = conf_int[, 2]
                        )
    
    ## Get model statistics
    conc_values <- summary(model)$concordance
    gmodel <- list(
        nevent = model$nevent,
        p_value_log = as.numeric(summary(model)$sctest["pvalue"]),
        AIC = stats::AIC(model),
        concordance = as.numeric(conc_values["C"]),
        concordance.se = as.numeric(conc_values["se(C)"])
    )
    
    ## Format events and AIC with commas
    gmodel$nevent_formatted <- format(gmodel$nevent, big.mark = ",", scientific = FALSE)
    gmodel$AIC_formatted <- format(round(gmodel$AIC, 2), big.mark = ",", scientific = FALSE, nsmall = 2)

    ## Calculate 95% CI for concordance if both concordance and SE are available
    if(!is.null(gmodel$concordance) && !is.na(gmodel$concordance) &&
       !is.null(gmodel$concordance.se) && !is.na(gmodel$concordance.se)) {
        gmodel$concordance.lower <- gmodel$concordance - 1.96 * gmodel$concordance.se
        gmodel$concordance.upper <- gmodel$concordance + 1.96 * gmodel$concordance.se
    } else {
        gmodel$concordance.lower <- NA
        gmodel$concordance.upper <- NA
    }
    
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
            data.table::data.table(var = var, level = "-", Freq = nrow(data), pos = 1, var_order = i)
        }
        else {
            vars = grep(paste0("^", var, "*."), coef$term, value=TRUE)
            data.table::data.table(var = vars, level = "-", Freq = nrow(data),
                                   pos = seq_along(vars), var_order = i)
        }
    })
    
    allTermsDF <- data.table::rbindlist(allTerms)

    data.table::setnames(allTermsDF, c("var", "level", "N", "pos", "var_order"))

## Add events
    formula_terms <- all.vars(model$formula)
    
    if (length(formula_terms) >= 2) {
        ## Assuming standard Surv(time, status) format
        event_var <- formula_terms[2]
        
        ## Calculate events for each term
        allTermsDF[, Events := {
            event_data <- data[[event_var]]
            
            ## Handle factor event indicators (rare but possible)
            if (is.factor(event_data)) {
                ## Convert to binary (assuming second level is the event)
                event_binary <- as.numeric(event_data) == 2
            } else {
                ## Already numeric
                event_binary <- event_data
            }
            
            if (level == "-") {
                ## Continuous variable - total events
                sum(event_binary, na.rm = TRUE)
            } else {
                ## Factor level - events in that group
                sum(event_binary[data[[var]] == level], na.rm = TRUE)
            }
        }, by = .(var, level)]
    } else {
        ## No events column if can't extract event variable
        allTermsDF[, Events := NA_integer_]
    }
    
    ## Create matching index
    allTermsDF[, inds := ifelse(level == "-", var, paste0(var, level))]
    
    ## Process coefficients
    coef[, term := gsub(term, pattern = "`", replacement = "-")]
    coef[, inds := term]
    
    ## Merge data
    toShow <- merge(allTermsDF, coef, by.x = "inds", by.y = "inds", all.x = TRUE, sort = FALSE)
    
    ## Sort by variable order first, then position within variable
    data.table::setorder(toShow, var_order, pos)
    
    ## Add variable-based shading indicator
    toShow[, shade_group := var_order %% 2]
    
    ## Select columns
    toShow <- toShow[, .(var, level, N, Events, p_value, estimate, conf_low, conf_high, pos, var_order, shade_group)]
    
    ## Format the exponential values
    toShowExpClean <- data.table::copy(toShow)
    
    ## Create formatted columns for display
    toShowExpClean[, hr := ifelse(is.na(estimate), 
                                  NA_real_,
                                  exp(estimate))]
    toShowExpClean[, hr_formatted := ifelse(is.na(estimate), 
                                            ref_label,
                                            format(round(exp(estimate), digits), nsmall = digits))]
    toShowExpClean[, conf_low_formatted := ifelse(is.na(conf_low), 
                                                  NA_character_,
                                                  format(round(exp(conf_low), digits), nsmall = digits))]
    toShowExpClean[, conf_high_formatted := ifelse(is.na(conf_high), 
                                                   NA_character_,
                                                   format(round(exp(conf_high), digits), nsmall = digits))]
    
    ## Format p-values
    toShowExpClean[, p_formatted := ifelse(is.na(p_value), 
                                           NA_character_,
                                    ifelse(p_value < 0.001, 
                                           "< 0.001",
                                           format(round(p_value, 3), nsmall = 3)))]
    
    ## Create the combined HR string with expression for italic p
    toShowExpClean[, hr_string_expr := ifelse(
                         is.na(estimate),
                         paste0("'", ref_label, "'"),
                                       ifelse(p_value < 0.001,
                                              paste0("'", hr_formatted, " (", conf_low_formatted, "-", 
                                                     conf_high_formatted, "); '*~italic(p)~'< 0.001'"),
                                              paste0("'", hr_formatted, " (", conf_low_formatted, "-", 
                                                     conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'"))
                     )]
    
    ## Format N and events with thousands separator
    toShowExpClean[, n_formatted := format(N, big.mark = ",", scientific = FALSE)]
    toShowExpClean[, events_formatted := format(Events, big.mark = ",", scientific = FALSE)]
    
    ## Clean up variable names for display
    toShowExpClean[, var_display := as.character(var)]

    ## Apply labels with priority: var_labels > attributes > variable name
    for(v in unique(toShowExpClean$var)) {
        if(v %in% toShowExpClean$var) {
            ## Check custom labels first (highest priority)
            if(!is.null(var_labels) && v %in% names(var_labels)) {
                toShowExpClean[var == v, var_display := var_labels[v]]
            }
            ## Then check for attributes
            else if(!is.null(attr(data[[v]], "label"))) {
                toShowExpClean[var == v, var_display := attr(data[[v]], "label")]
            }
            ## Otherwise keep the variable name as-is
        }
    }

    toShowExpClean[duplicated(var), var_display := ""]
    
    ## Handle missing estimates for plotting
    toShowExpClean[is.na(estimate), estimate := 0]
    
    ## Reorder (flip) - but maintain the variable grouping
    toShowExpClean <- toShowExpClean[order(nrow(toShowExpClean):1)]
    toShowExpClean[, x_pos := .I]
    
    ## Add shade_group based on var_order for alternating by variable
    toShowExpClean[, shade_group := var_order %% 2]
    
    ## Number of rows for later adjustments
    n_rows <- nrow(toShowExpClean)
    
    ## Calculate plot ranges with better handling of extreme cases
    rangeb <- range(toShow$conf_low, toShow$conf_high, na.rm = TRUE)
    breaks <- grDevices::axisTicks(rangeb/2, log = TRUE, nint = 7)
    
    ## Get min and max CI values
    min_ci <- min(toShow$conf_low, na.rm = TRUE)
    max_ci <- max(toShow$conf_high, na.rm = TRUE)
    
    ## Check for one-sided case BEFORE modifying range
    is_one_sided <- (min_ci > 0) || (max_ci < 0)
    
    ## Always ensure HR = 1 is included in breaks and range
    if(!1 %in% breaks) {
        breaks <- sort(unique(c(breaks, 1)))  # Add HR = 1
    }
    
    ## If all values are on one side of 1, extend range to include it properly
    if(min_ci > 0) {  # All HRs > 1
        rangeb[1] <- log(0.9)  # Extend to exactly 0.9
    } else if(max_ci < 0) {  # All HRs < 1
        rangeb[2] <- log(1.1)  # Extend to exactly 1.1
    }
    
    ## Now filter breaks based on the extended range
    breaks <- breaks[breaks >= exp(rangeb[1]) & breaks <= exp(rangeb[2])]
    
    ## Calculate dynamic column widths based on content
    max_var_len <- max(nchar(toShowExpClean$var_display) * 0.75, nchar("Variable"), na.rm = TRUE)
    max_level_len <- max(nchar(toShowExpClean$level) * 0.75, nchar("Group"), na.rm = TRUE)
    max_n_len <- max(nchar(toShowExpClean$n_formatted), nchar("n"), na.rm = TRUE)
    
    ## For HR column, calculate without expressions
    hr_display_len <- ifelse(toShowExpClean$hr_string_expr == paste0("'", ref_label, "'"),
                             nchar(ref_label),
                             nchar(paste0(toShowExpClean$hr_formatted, " (", 
                                          toShowExpClean$conf_low_formatted, "-",
                                          toShowExpClean$conf_high_formatted, "); p = ",
                                          toShowExpClean$p_formatted)))
    max_hr_len <- max(hr_display_len, nchar("HR (95% CI); p-value"), na.rm = TRUE) * 1.1
    
    ## Calculate total character width needed
    total_text_chars <- max_var_len + max_level_len + max_n_len + max_hr_len + 4
    
    ## Calculate optimal tbl_width if not provided
    if(is.null(tbl_width)) {
        
        ## Calculate required width for text in inches
        char_to_inch <- 0.08 * font_size
        text_width_needed <- total_text_chars * char_to_inch
        
        ## Calculate the actual forest plot range based on data
        min_ci_value <- min(toShow$conf_low, na.rm = TRUE)
        max_ci_value <- max(toShow$conf_high, na.rm = TRUE)
        
        ## The actual forest plot width is from min CI to max CI plus some padding
        forest_data_range <- max_ci_value - min_ci_value
        
        ## The forest needs significant space for axis, ticks, labels, and padding
        forest_width_min <- max(5, forest_data_range * 1.2 + 2)
        
        if(!is.null(plot_width)) {
            
            ## If plot width is specified, calculate proportion needed for text
            available_for_text <- plot_width - forest_width_min
            if(available_for_text > text_width_needed) {
                tbl_width <- text_width_needed / plot_width
            } else {
                
                ## Plot width might be too small, maximize text space
                tbl_width <- min(0.5, available_for_text / plot_width)
                warning(paste0("Specified plot width (", plot_width, 
                               ") may be too narrow. Consider width of at least ", 
                               round(text_width_needed + forest_width_min, 1), " inches."))
            }
        } else {
            
            ## No plot width specified, calculate optimal width
            ## Base calculation on total character count and forest needs
            text_proportion_needed <- text_width_needed / (text_width_needed + forest_width_min * 1.2)
            
            ## Adjust based on number of rows
            if(n_rows < 10) {
                tbl_width <- text_proportion_needed * 0.95
            } else if(n_rows > 30) {
                tbl_width <- text_proportion_needed * 0.85
            } else {
                tbl_width <- text_proportion_needed * 0.9
            }
        }
        
        ## Ensure tbl_width is within reasonable bounds
        tbl_width <- max(0.25, min(0.55, tbl_width))  # Increased max to 0.55
        
        ## Apply reduction for one-sided cases
        if(is_one_sided) {
            tbl_width <- tbl_width * 0.85  # Reduce by 15% for one-sided plots
        }
    }
    
    ## Calculate recommended plot dimensions
    ## Ensure minimum height for title and footnotes
    rec_height <- max(5, min(20, 3 + n_rows * 0.25))
    
    if(!is.null(plot_width)) {
        rec_width <- plot_width
        if(!is.null(plot_height)) {
            rec_height <- plot_height
        }
    } else {
        
        ## Calculate width based on actual space needs
        char_to_inch <- 0.08 * font_size
        text_width_needed <- total_text_chars * char_to_inch
        
        ## Calculate forest width based on data range with extra safety margin
        min_ci_value <- min(toShow$conf_low, na.rm = TRUE)
        max_ci_value <- max(toShow$conf_high, na.rm = TRUE)
        forest_data_range <- max_ci_value - min_ci_value
        forest_width_needed <- max(6, forest_data_range * 1.2 + 3)
        
        ## Total width with extra safety margin to prevent overlap
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
    
    ## Calculate how much space we need for the table relative to the forest plot
    ## This creates a virtual space to the left of the smallest HR value
    if(is.null(tbl_width)) {
        table_space_factor <- 1.5  # Default: table takes 1.5x the width of the forest plot range
    } else {
        
        ## Calculate based on tbl_width proportion
        table_space_factor <- tbl_width / (1 - tbl_width)
    }
    
    ## Extend the lower range to accommodate the table
    rangeplot[1] <- rangeb[1] - (range_width * table_space_factor * 1.3)  # Reduced multiplier to bring table closer
    rangeplot[2] <- rangeb[2] + (range_width * 0.05)  # Reduced extension on the right
    
    ## Calculate positions for text columns within the extended range
    ## Place them in the left portion of the plot
    width <- diff(rangeplot)

    ## Handle show_n_events parameter
    if (is.null(show_n_events)) {
        show_n_events <- character(0)
    } else {
        show_n_events[show_n_events == "events"] <- "Events"
    }
    
    ## Start positions from the left edge of the extended range
    y_variable <- rangeplot[1] + col_width_var * width
    y_level <- y_variable + col_width_level * width

    if ("n" %in% show_n_events) {
        y_n <- y_level + col_width_n * width
        next_pos <- y_n
    } else {
        next_pos <- y_level
    }
    
    if ("Events" %in% show_n_events) {
        y_events <- next_pos + col_width_events * width
        y_hr <- y_events + col_width_hr * width
    } else {
        y_hr <- next_pos + col_width_hr * width
    }

    ## Ensure there's a small gap between the HR column and the forest plot
    forest_start <- rangeb[1] - 0.05 * range_width
    
    ## Font size for annotations
    annot_font <- font_size * annot_size  # 11 pt * 0.352778 mm/pt
    
    ## Font size for headers (column headers)
    header_font <- font_size * header_size
    
    ## Custom ticks data
    ticks_df <- data.frame(
        x = -0.5,
        xend = -0.7,
        y = breaks,
        yend = breaks
    )
    
    ## Format the global p-value for display
    global_p_formatted <- if(as.numeric(gmodel$p_value_log) < 0.001) {
                              "< 0.001"
                          } else {
                              format.pval(gmodel$p_value_log, eps = ".001")
                          }
    
    ## Build concordance string based on whether CI is available
    concordance_string <- if(!is.null(gmodel$concordance) && !is.na(gmodel$concordance)) {
                              if(!is.null(gmodel$concordance.lower) && !is.na(gmodel$concordance.lower) && 
                                 !is.null(gmodel$concordance.upper) && !is.na(gmodel$concordance.upper)) {
                                  paste0("Concordance Index: ", round(gmodel$concordance, 2), 
                                         " (95% CI ", round(gmodel$concordance.lower, 2), "-", 
                                         round(gmodel$concordance.upper, 2), ")")
                              } else {
                                  paste0("Concordance Index: ", round(gmodel$concordance, 2))
                              }
                          } else {
                              "Concordance Index: Not available"
                          }
    
    p <- ggplot2::ggplot(toShowExpClean, ggplot2::aes(x_pos, exp(estimate))) +
        
        ## Shading rectangles - extend to cover table area too
        ggplot2::geom_rect(ggplot2::aes(xmin = x_pos - .5, xmax = x_pos + .5,
                                        ymin = exp(rangeplot[1]), ymax = exp(rangeplot[2]),
                                        fill = ordered(shade_group + 1))) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
        ggplot2::scale_size_continuous(range = c(1, 6), guide = "none") +
        ggplot2::scale_fill_manual(values = c("#FFFFFF", "#EEEEEE"), guide = "none") +
        
        ## Forest plot elements
        ggplot2::geom_point(ggplot2::aes(size = N), pch = 22, color = "#000000", fill = color) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = exp(conf_low), ymax = exp(conf_high)), width = 0.15) +
        
        ## Y-axis for forest plot
        ggplot2::annotate(geom = "segment",
                          x = -0.5, xend = -0.5,
                          y = exp(rangeb[1]),
                          yend = exp(rangeb[2]),
                          color = "#000000", linewidth = 1) +
        
        ## Reference line at HR = 1
        ggplot2::annotate(geom = "segment", 
                          x = -0.5, xend = max(toShowExpClean$x_pos) + 0.5, 
                          y = 1, yend = 1, 
                          linetype = "longdash") +
        
        ## Ticks
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
        
        ## Set coordinate system with extended limits
        ggplot2::coord_flip(ylim = exp(rangeplot)) +
        ggplot2::ggtitle(title) +
        ggplot2::scale_y_log10(name = "Hazard Ratio",
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
        
        ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = exp(y_variable),
                          label = toShowExpClean$var_display, fontface = "bold", hjust = 0,
                          size = annot_font) +
        
        ## Group/level column
        
        ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.5, y = exp(y_level),
                          label = "Group", fontface = "bold", hjust = 0,
                          size = header_font) +
        
        ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = exp(y_level),
                          label = toShowExpClean$level, hjust = 0,
                          size = annot_font) +
                
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
        ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.4, y = exp(y_hr),
                          label = "bold('aHR (95% CI); '*bolditalic(p)*'-value')",
                         hjust = 0, size = header_font, parse = TRUE) +
        
        ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = exp(y_hr),
                          label = toShowExpClean$hr_string_expr, hjust = 0,
                          size = annot_font, parse = TRUE) +
        
        ## X-axis label
        ggplot2::annotate(geom = "text", x = -1.5, y = 1,
                          label = "Hazard Ratio", fontface = "bold",
                          hjust = 0.5, vjust = 2, size = annot_font * 1.5) +
        
        ## Model statistics footer
        ggplot2::annotate(geom = "text", x = 0.5, y = exp(y_variable),
                          label = paste0("Total events: ", gmodel$nevent_formatted,
                                         "\nGlobal log-rank p: ", global_p_formatted,
                                         "\n", concordance_string,
                                         "\nAIC: ", gmodel$AIC_formatted),
                          size = annot_font, hjust = 0, vjust = 1.2, fontface = "italic")
    
    ## Add recommended dimensions as an attribute
    attr(p, "recommended_dims") <- list(width = rec_width, height = rec_height)
    
    ## Return the plot
    return(p)
}
