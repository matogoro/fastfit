#' Create a Forest Plot for Linear Models
#'
#' Generates a publication-ready forest plot combining a data table and 
#' graphical representation of coefficients from a linear model.
#' The plot includes variable names, group levels, sample sizes, coefficients
#' with confidence intervals, and p-values.
#'
#' @param model A linear model object of class \code{lm}.
#' @param data A data frame or data.table containing the original data used to 
#'   fit the model. If NULL (default), the function attempts to extract data 
#'   from the model object.
#' @param title Character string specifying the plot title. Default is "Linear Model".
#' @param effect_label Character string for the effect measure label. Default is "Coefficient".
#' @param digits Integer specifying the number of decimal places for coefficients 
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
#' @param col_width_coef Numeric value specifying the proportion of plot width 
#'   allocated to the coefficient column. Default is 0.05.
#' @param show_n Logical. Whether to show the sample size column. Default is TRUE.
#' @param ref_label Character string to display for reference categories. 
#'   Default is "reference".
#' @param var_labels Named character vector for custom variable labels. Names should 
#'   match variable names in the model, values are the display labels.
#' @param color Character string specifying the color for coefficient point estimates. 
#'   Default is "#5A8F5A" (green).
#'
#' @return A ggplot object containing the forest plot. The object has an attribute
#'   "recommended_dims" with suggested width and height in inches.
#'
#' @export
#' @import data.table
#' @import ggplot2
#' @importFrom stats AIC confint nobs pf
#' @importFrom grDevices axisTicks
#'
#' @seealso 
#' \code{\link[stats]{lm}} for fitting linear models,
#' \code{\link[ggplot2]{ggplot}} for the underlying plotting system
#'
lmforest <- function(model, data = NULL,
                     title = "Linear Model",
                     effect_label = "Coefficient",
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
                     col_width_coef = 0.05,
                     show_n = TRUE,
                     ref_label = "reference",
                     var_labels = NULL,
                     color = "#5A8F5A") {
    
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
    
    ## Check model class
    stopifnot(inherits(model, "lm") && !inherits(model, "glm"))
    
    ## Get data
    if(is.null(data)){
        warning("The `data` argument is not provided. Data will be extracted from model fit.")
        data <- model$data
        if (is.null(data))
            data <- model$model
        if (is.null(data))
            stop("The `data` argument should be provided either to lmforest or lm.")
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
                            statistic = coef_summary[, "t value"],
                            p_value = coef_summary[, "Pr(>|t|)"],
                            conf_low = conf_int[, 1],
                            conf_high = conf_int[, 2]
                        )
    
    ## Get model statistics
    model_summary <- summary(model)
    gmodel <- list(
        nobs = nobs(model),
        r_squared = model_summary$r.squared,
        adj_r_squared = model_summary$adj.r.squared,
        f_statistic = model_summary$fstatistic[1],
        f_df1 = model_summary$fstatistic[2],
        f_df2 = model_summary$fstatistic[3],
        AIC = stats::AIC(model)
    )
    
    ## Calculate F-test p-value
    gmodel$f_pvalue <- pf(gmodel$f_statistic, gmodel$f_df1, gmodel$f_df2, lower.tail = FALSE)
    
    ## Format for display
    gmodel$nobs_formatted <- format(gmodel$nobs, big.mark = ",", scientific = FALSE)
    gmodel$AIC_formatted <- format(round(gmodel$AIC, 2), big.mark = ",", scientific = FALSE, nsmall = 2)
    
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
    
    ## Create matching index
    allTermsDF[, inds := ifelse(level == "-", var, paste0(var, level))]
    
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
    toShow <- toShow[, .(var, level, N, p_value, estimate, conf_low, conf_high, pos, var_order, shade_group)]
    
    ## Format the values
    toShowExpClean <- data.table::copy(toShow)
    
    ## Create formatted columns for display
    toShowExpClean[, effect := estimate]
    toShowExpClean[, effect_formatted := ifelse(is.na(estimate), 
                                                ref_label,
                                                format(round(estimate, digits), nsmall = digits))]
    toShowExpClean[, conf_low_formatted := ifelse(is.na(conf_low), 
                                                  NA_character_,
                                                  format(round(conf_low, digits), nsmall = digits))]
    toShowExpClean[, conf_high_formatted := ifelse(is.na(conf_high), 
                                                   NA_character_,
                                                   format(round(conf_high, digits), nsmall = digits))]
    
    ## Format p-values
    toShowExpClean[, p_formatted := ifelse(is.na(p_value), 
                                           NA_character_,
                                    ifelse(p_value < 0.001, 
                                           "< 0.001",
                                           format(round(p_value, 3), nsmall = 3)))]
    
    ## Create the combined effect string with expression for italic p
    toShowExpClean[, effect_string_expr := fcase(
                         is.na(estimate), paste0("'", ref_label, "'"),
                         
                         p_value < 0.001 & (conf_low < 0 | conf_high < 0),
                         paste0("'", effect_formatted, " (", conf_low_formatted, ", ", 
                                conf_high_formatted, "); '*~italic(p)~'< 0.001'"),
                         
                         p_value < 0.001,
                         paste0("'", effect_formatted, " (", conf_low_formatted, "-", 
                                conf_high_formatted, "); '*~italic(p)~'< 0.001'"),
                         
                         conf_low < 0 | conf_high < 0,
                         paste0("'", effect_formatted, " (", conf_low_formatted, ", ", 
                                conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'"),
                         
                         default = paste0("'", effect_formatted, " (", conf_low_formatted, "-", 
                                          conf_high_formatted, "); '*~italic(p)~'= ", p_formatted, "'")
                     )]
    
    ## Format N with thousands separator
    toShowExpClean[, n_formatted := format(N, big.mark = ",", scientific = FALSE)]
    
    ## Clean up variable names for display
    toShowExpClean[, var_display := as.character(var)]
    
    ## Apply labels with priority: var_labels > attributes > variable name
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
    
    ## Calculate plot ranges
    rangeb <- range(toShow$conf_low, toShow$conf_high, na.rm = TRUE)
    breaks <- pretty(rangeb, n = 7)
    breaks <- breaks[breaks >= rangeb[1] & breaks <= rangeb[2]]
    
    ## Check for one-sided case
    is_one_sided <- (min(rangeb) > 0) || (max(rangeb) < 0)
    
    ## Ensure 0 is included for coefficients
    if(!0 %in% breaks) {
        breaks <- sort(unique(c(breaks, 0)))
    }
    
    ## Extend range if needed
    if(min(rangeb) > 0) {
        rangeb[1] <- -0.1 * abs(max(rangeb))
    } else if(max(rangeb) < 0) {
        rangeb[2] <- 0.1 * abs(min(rangeb))
    }
    
    reference_value <- 0
    
    ## Calculate dynamic column widths based on content
    max_var_len <- max(nchar(toShowExpClean$var_display) * 0.75, nchar("Variable"), na.rm = TRUE)
    max_level_len <- max(nchar(toShowExpClean$level) * 0.75, nchar("Group"), na.rm = TRUE)
    max_n_len <- max(nchar(toShowExpClean$n_formatted), nchar("n"), na.rm = TRUE)
    
    ## For effect column
    effect_display_len <- ifelse(toShowExpClean$effect_string_expr == paste0("'", ref_label, "'"),
                                 nchar(ref_label),
                                 nchar(paste0(toShowExpClean$effect_formatted, " (", 
                                              toShowExpClean$conf_low_formatted, ", ",
                                              toShowExpClean$conf_high_formatted, "); p = ",
                                              toShowExpClean$p_formatted)))
    max_effect_len <- max(effect_display_len, nchar("Coefficient (95% CI); p-value"), na.rm = TRUE) * 1.1
    
    ## Calculate total character width needed
    total_text_chars <- max_var_len + max_level_len + (if(show_n) max_n_len else 0) + max_effect_len + 4
    
    ## Calculate optimal tbl_width if not provided
    if(is.null(tbl_width)) {
        char_to_inch <- 0.08 * font_size
        text_width_needed <- total_text_chars * char_to_inch
        
        min_ci_value <- min(toShow$conf_low, na.rm = TRUE)
        max_ci_value <- max(toShow$conf_high, na.rm = TRUE)
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
        
        min_ci_value <- min(toShow$conf_low, na.rm = TRUE) 
        max_ci_value <- max(toShow$conf_high, na.rm = TRUE)
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
    
    ## Start positions from the left edge
    y_variable <- rangeplot[1] + col_width_var * width
    y_level <- y_variable + col_width_level * width
    
    if (show_n) {
        y_n <- y_level + col_width_n * width
        y_coef <- y_n + col_width_coef * width
    } else {
        y_coef <- y_level + col_width_coef * width
    }
    
    ## Font sizes
    annot_font <- font_size * annot_size
    header_font <- font_size * header_size
    
    ## Custom ticks data
    ticks_df <- data.frame(
        x = -0.5,
        xend = -0.7,
        y = breaks,
        yend = breaks
    )
    
    ## Create the plot
    p <- ggplot2::ggplot(toShowExpClean, ggplot2::aes(x_pos, estimate)) +
        
        ## Shading rectangles
        ggplot2::geom_rect(ggplot2::aes(xmin = x_pos - .5, xmax = x_pos + .5,
                                        ymin = rangeplot[1], ymax = rangeplot[2],
                                        fill = ordered(shade_group + 1))) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
        ggplot2::scale_size_continuous(range = c(1, 6), guide = "none") +
        ggplot2::scale_fill_manual(values = c("#FFFFFF", "#EEEEEE"), guide = "none") +
        
        ## Forest plot elements
        ggplot2::geom_point(ggplot2::aes(size = N), pch = 22, color = "#000000", fill = color) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = conf_low, ymax = conf_high), width = 0.15) +
        
        ## Y-axis for forest plot
        ggplot2::annotate(geom = "segment",
                          x = -0.5, xend = -0.5,
                          y = rangeb[1],
                          yend = rangeb[2],
                          color = "#000000", linewidth = 1) +
        
        ## Reference line at 0
        ggplot2::annotate(geom = "segment", 
                          x = -0.5, xend = max(toShowExpClean$x_pos) + 0.5, 
                          y = reference_value, yend = reference_value, 
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
        
        ## Set coordinate system
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
        ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = y_variable,
                          label = toShowExpClean$var_display, fontface = "bold", hjust = 0,
                          size = annot_font) +
        
        ## Group/level column
        ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.5, y = y_level,
                          label = "Group", fontface = "bold", hjust = 0,
                          size = header_font) +
        ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = y_level,
                          label = toShowExpClean$level, hjust = 0,
                          size = annot_font) +
        
        ## N column (conditional)
        {if (show_n) {
             list(
                 ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.5, y = y_n,
                                   label = "n", fontface = "bold.italic", hjust = 0.5,
                                   size = header_font),
                 ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = y_n,
                                   label = toShowExpClean$n_formatted, hjust = 0.5,
                                   size = annot_font)
             )
         }} +
        
        ## Coefficient column
        ggplot2::annotate(geom = "text", x = max(toShowExpClean$x_pos) + 1.4, y = y_coef,
                          label = paste0("bold('Coefficient (95% CI); '*bolditalic(p)*'-value')"),
                          hjust = 0, size = header_font, parse = TRUE) +
        
        ggplot2::annotate(geom = "text", x = toShowExpClean$x_pos, y = y_coef,
                          label = toShowExpClean$effect_string_expr, hjust = 0,
                          size = annot_font, parse = TRUE) +
        
        ## X-axis label
        ggplot2::annotate(geom = "text", x = -1.5, y = reference_value,
                          label = effect_label, fontface = "bold",
                          hjust = 0.5, vjust = 2, size = annot_font * 1.5) +
        
        ## Model statistics at bottom
        ggplot2::annotate(geom = "text", x = 0.5, y = y_variable,
                          label = paste0("Observations: ", gmodel$nobs_formatted,
                                         "\nR²: ", round(gmodel$r_squared, 3),
                                         " (Adjusted: ", round(gmodel$adj_r_squared, 3), ")",
                                         "\nF-statistic: ", round(gmodel$f_statistic, 2),
                                         " (df1=", gmodel$f_df1, ", df2=", gmodel$f_df2,
                                         "); p ", ifelse(gmodel$f_pvalue < 0.001, "< 0.001", 
                                                         paste0("= ", round(gmodel$f_pvalue, 3))),
                                         "\nAIC: ", gmodel$AIC_formatted),
                          size = annot_font * 0.8, hjust = 0, vjust = 1.2, fontface = "italic")
    
    ## Add recommended dimensions as an attribute
    attr(p, "recommended_dims") <- list(width = rec_width, height = rec_height)
    
    ## Return the plot
    return(p)
}
