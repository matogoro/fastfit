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
#' @param show_n Logical. Whether to show the sample size column. Default is TRUE.
#' @param indent_groups Logical. Whether to indent group levels under variables.
#'   Default is FALSE.
#' @param condense_table Logical. Whether to condense binary variables into single rows.
#'   Forces indent_groups = TRUE. Default is FALSE.
#' @param center_padding Numeric value specifying the padding width between the table
#' and the forest plot. Default is 4.
#' @param zebra_stripes Logical. Whether to use alternating background shading 
#'   for different variables to improve readability. Default is TRUE.
#' @param ref_label Character string to display for reference categories. 
#'   Default is "reference".
#' @param var_labels Named character vector for custom variable labels. Names should 
#'   match variable names in the model, values are the display labels.
#' @param units Units used for measurement. Default is inches ("in").
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
                     show_n = TRUE,
                     indent_groups = FALSE,
                     condense_table = FALSE,
                     center_padding = 4,
                     zebra_stripes = TRUE,
                     ref_label = "reference",
                     var_labels = NULL,
                     units = "in",
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

    ## Internally work in inches
    if (!is.null(plot_width) && units != "in") {
        plot_width <- convert_units(plot_width, from = units, to = "inches")
    }
    
    
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
    
    ## Calculate total observations and percentage analyzed
    total_obs <- nrow(data)
    gmodel$pct_analyzed <- (gmodel$nobs / total_obs) * 100
    
    ## Format for display
    gmodel$nobs_formatted <- format(gmodel$nobs, big.mark = ",", scientific = FALSE)
    gmodel$nobs_with_pct <- paste0(gmodel$nobs_formatted, " (", 
                                   sprintf("%.1f%%", gmodel$pct_analyzed), ")")
    gmodel$AIC_formatted <- format(round(gmodel$AIC, 2), big.mark = ",", scientific = FALSE, nsmall = 2)
    
    ## Extract statistics for every variable - preserving order
    all_terms <- lapply(seq_along(terms), function(i){
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
    
    all_terms_df <- data.table::rbindlist(all_terms)
    data.table::setnames(all_terms_df, c("var", "level", "N", "pos", "var_order"))

    ## Apply condensing and indenting if requested
    if (condense_table || indent_groups) {
        if (condense_table) {
            indent_groups <- TRUE
        }
        
        all_terms_df[, inds := ifelse(level == "-", var, paste0(var, level))]
        orig_inds_map <- data.table::copy(all_terms_df[, .(var, level, inds, N)])
        
        processed_rows <- list()
        row_counter <- 1
        unique_vars <- unique(all_terms_df[, var])
        
        for (v in unique_vars) {
            var_rows <- all_terms_df[var == v]
            
            if (nrow(var_rows) == 1) {
                                        # Single-level variable (continuous)
                processed_rows[[row_counter]] <- var_rows
                row_counter <- row_counter + 1
            } else {
                                        # Multi-level variable (categorical)
                is_binary <- nrow(var_rows) == 2
                
                if (condense_table && is_binary) {
                                        # Only condense if BOTH flags are set
                    non_ref_row <- var_rows[2]
                    condensed_row <- data.table::copy(non_ref_row)
                    condensed_row[, var := paste0(v, " (", level, ")")]
                    condensed_row[, level := "-"]
                    processed_rows[[row_counter]] <- condensed_row
                    row_counter <- row_counter + 1
                } else if (indent_groups) {
                                        # Add header for ANY multi-level variable when indenting
                    header_row <- data.table::data.table(
                                                  var = v,
                                                  level = "-",
                                                  N = NA_integer_,
                                                  pos = var_rows$pos[1],
                                                  var_order = var_rows$var_order[1],
                                                  inds = NA_character_
                                              )
                    processed_rows[[row_counter]] <- header_row
                    row_counter <- row_counter + 1
                    
                                        # Add indented levels
                    for (i in 1:nrow(var_rows)) {
                        group_row <- data.table::copy(var_rows[i])
                        group_row[, var := paste0("    ", level)]
                        group_row[, level := "-"]
                        processed_rows[[row_counter]] <- group_row
                        row_counter <- row_counter + 1
                    }
                } else {
                                        # Standard layout
                    for (i in 1:nrow(var_rows)) {
                        processed_rows[[row_counter]] <- var_rows[i]
                        row_counter <- row_counter + 1
                    }
                }
            }
        }
        
        all_terms_df <- data.table::rbindlist(processed_rows, fill = TRUE)
        
        for (i in 1:nrow(all_terms_df)) {
            current_var <- all_terms_df$var[i]
            
            if (is.na(all_terms_df$inds[i]) && !grepl("^    ", current_var) && !grepl("\\(", current_var)) {
                next
            }
            
            if (grepl("\\(", current_var)) {
                orig_var <- gsub(" \\(.*\\)", "", current_var)
                orig_level <- gsub(".*\\((.*)\\)", "\\1", current_var)
                
                matching <- orig_inds_map[var == orig_var & level == orig_level]
                if (nrow(matching) > 0) {
                    all_terms_df[i, `:=`(inds = matching$inds[1],
                                         N = matching$N[1])]
                }
            }
            else if (grepl("^    ", current_var)) {
                clean_level <- gsub("^    ", "", current_var)
                
                parent_var <- NA_character_
                for (j in (i-1):1) {
                    if (!grepl("^    ", all_terms_df$var[j]) && all_terms_df$var[j] != "") {
                        parent_var <- gsub(" \\(.*\\)", "", all_terms_df$var[j])
                        break
                    }
                }
                
                if (!is.na(parent_var)) {
                    matching <- orig_inds_map[var == parent_var & level == clean_level]
                    if (nrow(matching) > 0) {
                        all_terms_df[i, `:=`(inds = matching$inds[1],
                                             N = matching$N[1])]
                    }
                }
            }
        }
    } else {
        all_terms_df[, inds := ifelse(level == "-", var, paste0(var, level))]
    }
    
    ## Process coefficients
    coef[, term := gsub(term, pattern = "`", replacement = "")]
    coef[, inds := term]
    
    ## Merge data
    to_show <- merge(all_terms_df, coef, by.x = "inds", by.y = "inds", all.x = TRUE, sort = FALSE)
    
    ## Sort by variable order first, then position within variable
    data.table::setorder(to_show, var_order, pos)
    
    ## Add variable-based shading indicator (zebra stripes)
    if (zebra_stripes) {
        to_show[, shade_group := var_order %% 2]
        shade_colors <- c("#FFFFFF", "#EEEEEE")
    } else {
        to_show[, shade_group := 0]
        shade_colors <- c("#FFFFFF", "#FFFFFF")
    }
    
    ## Select columns
    to_show <- to_show[, .(var, level, N, p_value, estimate, conf_low, conf_high, pos, var_order, shade_group)]
    
    ## Format the values
    to_show_exp_clean <- data.table::copy(to_show)
    
    ## Create formatted columns for display
    to_show_exp_clean[, effect := estimate]
    to_show_exp_clean[, effect_formatted := ifelse(is.na(N) & is.na(estimate),
                                                   "",
                                            ifelse(is.na(estimate), 
                                                   ref_label,
                                                   format(round(estimate, digits), nsmall = digits)))]
    to_show_exp_clean[, conf_low_formatted := ifelse(is.na(conf_low), 
                                                     NA_character_,
                                                     format(round(conf_low, digits), nsmall = digits))]
    to_show_exp_clean[, conf_high_formatted := ifelse(is.na(conf_high), 
                                                      NA_character_,
                                                      format(round(conf_high, digits), nsmall = digits))]
    
    ## Format p-values
    to_show_exp_clean[, p_formatted := ifelse(is.na(p_value), 
                                              NA_character_,
                                       ifelse(p_value < 0.001, 
                                              "< 0.001",
                                              format(round(p_value, 3), nsmall = 3)))]
    
    ## Create the combined effect string with expression for italic p
    to_show_exp_clean[, effect_string_expr := ifelse(
                            is.na(N) & is.na(estimate),
                            "''",
                            fcase(
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
                            )
                        )]
    
    ## Format N with thousands separator
    to_show_exp_clean[, n_formatted := ifelse(is.na(N), "", format(N, big.mark = ",", scientific = FALSE))]
    
    ## Clean up variable names for display
    to_show_exp_clean[, var_display := as.character(var)]
    
    if (indent_groups || condense_table) {
        to_show_exp_clean[, var_display := var]
        
        for (v in unique(to_show_exp_clean$var)) {
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
                            to_show_exp_clean[var == v, var_display := label]
                        } else {
                            to_show_exp_clean[var == v, var_display := paste0(label, " (", category, ")")]
                        }
                    } else {
                        to_show_exp_clean[var == v, var_display := label]
                    }
                }
            }
        }
        
        to_show_exp_clean[, level := ""]
    } else {
        for(v in unique(to_show_exp_clean$var)) {
            if(v %in% to_show_exp_clean$var) {
                if(!is.null(var_labels) && v %in% names(var_labels)) {
                    to_show_exp_clean[var == v, var_display := var_labels[v]]
                }
                else if(!is.null(attr(data[[v]], "label"))) {
                    to_show_exp_clean[var == v, var_display := attr(data[[v]], "label")]
                }
            }
        }
    }

    if (!indent_groups) {
        to_show_exp_clean[duplicated(var), var_display := ""]
    }
    
    ## Handle missing estimates for plotting
    to_show_exp_clean[is.na(estimate), estimate := 0]
    
    ## Reorder (flip) - but maintain the variable grouping
    to_show_exp_clean <- to_show_exp_clean[order(nrow(to_show_exp_clean):1)]
    to_show_exp_clean[, x_pos := .I]
    
    ## Add shade_group based on var_order for alternating by variable
    to_show_exp_clean[, shade_group := var_order %% 2]
    
    ## Number of rows for later adjustments
    n_rows <- nrow(to_show_exp_clean)
    
    ## Calculate plot ranges
    rangeb <- range(to_show$conf_low, to_show$conf_high, na.rm = TRUE)
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

    ## Calculate plot ranges
    rangeb <- range(to_show$conf_low, to_show$conf_high, na.rm = TRUE)
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

    ## Calculate layout using helper function
    layout <- calculate_forest_layout(
        to_show_exp_clean = to_show_exp_clean,
        show_n = show_n,
        show_events = FALSE,  # lmforest doesn't have events
        indent_groups = indent_groups,
        condense_table = condense_table,
        effect_label = effect_label,
        ref_label = ref_label,
        font_size = font_size,
        tbl_width = ifelse(is.null(tbl_width), 0.6, tbl_width),
        rangeb = rangeb,
        center_padding = center_padding
    )

    ## Set up the extended range for plotting
    rangeplot <- c(layout$rangeplot_start, rangeb[2] + diff(rangeb) * 0.05)

    ## Extract positions
    y_variable <- layout$positions$var
    if (!(indent_groups || condense_table)) {
        y_level <- layout$positions$level
    }
    if (show_n) {
        y_n <- layout$positions$n
    }
    y_coef <- layout$positions$effect  # Note: called y_coef in lmforest, not y_or

    ## Calculate recommended dimensions
    rec_height <- max(5, min(20, 3 + nrow(to_show_exp_clean) * 0.25))

    if(!is.null(plot_width)) {
        rec_width <- plot_width
        if(!is.null(plot_height)) {
            rec_height <- plot_height
        }
    } else {
        ## Use the calculated total width
        rec_width <- layout$total_width + 1.0  # Add margins
        
        ## Apply reasonable bounds based on layout
        if (condense_table || indent_groups) {
            if (!show_n) {
                rec_width <- max(7.5, min(11, rec_width))
            } else {
                rec_width <- max(9, min(14, rec_width))
            }
        } else {
            if (!show_n) {
                rec_width <- max(12, min(15, rec_width))
            } else {
                rec_width <- max(13, min(16, rec_width))
            }
        }
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
    p <- ggplot2::ggplot(to_show_exp_clean, ggplot2::aes(x_pos, estimate)) +
        
        ## Shading rectangles
        ggplot2::geom_rect(ggplot2::aes(xmin = x_pos - .5, xmax = x_pos + .5,
                                        ymin = rangeplot[1], ymax = rangeplot[2],
                                        fill = ordered(shade_group))) +
        ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.10, 0.05))) +
        ggplot2::scale_size_continuous(range = c(1, 6), guide = "none") +
        ggplot2::scale_fill_manual(values = shade_colors, guide = "none") +
        
        ## Forest plot elements
        ggplot2::geom_point(ggplot2::aes(size = N), pch = 22, color = "#000000", fill = color, na.rm = TRUE) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = conf_low, ymax = conf_high), width = 0.15) +
        
        ## Y-axis for forest plot
        ggplot2::annotate(geom = "segment",
                          x = -0.5, xend = -0.5,
                          y = rangeb[1],
                          yend = rangeb[2],
                          color = "#000000", linewidth = 1) +
        
        ## Reference line at 0
        ggplot2::annotate(geom = "segment", 
                          x = -0.5, xend = max(to_show_exp_clean$x_pos) + 0.5, 
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
        ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = y_variable,
                          label = "Variable", fontface = "bold", hjust = 0,
                          size = header_font) +
        
        {if (indent_groups || condense_table) {
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_variable,
                               label = to_show_exp_clean$var_display, 
                               fontface = ifelse(grepl("^    ", to_show_exp_clean$var_display), "plain", "bold"), 
                               hjust = 0,
                               size = annot_font)
         } else {
             ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_variable,
                               label = to_show_exp_clean$var_display, fontface = "bold", hjust = 0,
                               size = annot_font)
         }} +
        
        ## Group/level column
        {if (!(indent_groups || condense_table)) {
             list(
                 ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = y_level,
                                   label = "Group", fontface = "bold", hjust = 0,
                                   size = header_font),
                 ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_level,
                                   label = to_show_exp_clean$level, hjust = 0,
                                   size = annot_font)
             )
         }} +
        
        ## N column (conditional)
        {if (show_n) {
             list(
                 ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.5, y = y_n,
                                   label = "n", fontface = "bold.italic", hjust = 0.5,
                                   size = header_font),
                 ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_n,
                                   label = to_show_exp_clean$n_formatted, hjust = 0.5,
                                   size = annot_font)
             )
         }} +
        
        ## Coefficient column
        ggplot2::annotate(geom = "text", x = max(to_show_exp_clean$x_pos) + 1.4, y = y_coef,
                          label = paste0("bold('Coefficient (95% CI); '*bolditalic(p)*'-value')"),
                          hjust = 0, size = header_font, parse = TRUE) +
        
        ggplot2::annotate(geom = "text", x = to_show_exp_clean$x_pos, y = y_coef,
                          label = to_show_exp_clean$effect_string_expr, hjust = 0,
                          size = annot_font, parse = TRUE) +
        
        ## X-axis label
        ggplot2::annotate(geom = "text", x = -1.5, y = reference_value,
                          label = effect_label, fontface = "bold",
                          hjust = 0.5, vjust = 2, size = annot_font * 1.5) +
        
        ## Model statistics at bottom
        ggplot2::annotate(geom = "text", x = 0.5, y = y_variable,
                          label = paste0("Observations analyzed: ", gmodel$nobs_with_pct,
                                         "\nR²: ", round(gmodel$r_squared, 3),
                                         " (Adjusted: ", round(gmodel$adj_r_squared, 3), ")",
                                         "\nF-statistic: ", round(gmodel$f_statistic, 2),
                                         " (df1 = ", gmodel$f_df1, ", df2 = ", gmodel$f_df2,
                                         "); p ", ifelse(gmodel$f_pvalue < 0.001, "< 0.001", 
                                                         paste0("= ", round(gmodel$f_pvalue, 3))),
                                         "\nAIC: ", gmodel$AIC_formatted),
                          size = annot_font * 0.8, hjust = 0, vjust = 1.2, fontface = "italic")

    ## Convert units back for output if needed
    if (units != "in") {
        rec_width <- convert_units(rec_width, from = "inches", to = units)
        rec_height <- convert_units(rec_height, from = "inches", to = units)
    }

    ## Provide dimension recommendations
    if(is.null(plot_width) || is.null(plot_height)) {
        message(sprintf("Recommended plot dimensions: width = %.1f %s, height = %.1f %s",
                        rec_width, units, rec_height, units))
    }
    
    ## Add recommended dimensions as an attribute
    attr(p, "recommended_dims") <- list(width = rec_width, height = rec_height, units = units)

    ## Return the plot
    return(p)
}
