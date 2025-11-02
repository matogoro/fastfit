#' Convert between units
#' @keywords internal
convert_units <- function(value, from = "in", to = "in", dpi = 96) {
    to_inches <- c(
        "in" = 1,
        "cm" = 1/2.54,    # 1 cm = 1/2.54 inches
        "mm" = 1/25.4,    # 1 mm = 1/25.4 inches
        "px" = 1/dpi,     # 1 px = 1/dpi inches
        "pt" = 1/72       # 1 pt = 1/72 inches
    )
    
    value_in_inches <- value * to_inches[[from]]
    value_in_target <- value_in_inches / to_inches[[to]]
    
    return(value_in_target)
}

#' Calculate table layout for forest plots
#' @keywords internal
calculate_forest_layout <- function(to_show_exp_clean, show_n, show_events, 
                                    indent_groups, condense_table, 
                                    effect_label, ref_label, font_size, 
                                    tbl_width = 0.6, rangeb, center_padding) {
    
    char_to_inch <- 0.08 * font_size
    margin <- 0.3  # inches between columns
    
    ## Build list of active columns with their widths
    columns <- list()
    
    ## Always have variable column
    var_width_chars <- max(nchar(to_show_exp_clean$var_display), nchar("Variable"), na.rm = TRUE)
    columns$var <- var_width_chars * char_to_inch
    
    ## Conditionally add level column
    if (!(indent_groups || condense_table)) {
        level_width_chars <- max(nchar(to_show_exp_clean$level), nchar("Group"), na.rm = TRUE)
        columns$level <- level_width_chars * char_to_inch
    }
    
    ## Conditionally add n column
    if (show_n) {
        n_width_chars <- max(nchar(to_show_exp_clean$n_formatted), nchar("____"), na.rm = TRUE)
        columns$n <- n_width_chars * char_to_inch
    }
    
    ## Conditionally add events column
    if (show_events) {
        events_width_chars <- max(nchar(to_show_exp_clean$events_formatted), nchar("___"), na.rm = TRUE)
        columns$events <- events_width_chars * char_to_inch
    }
    
    ## Always have effect column
    effect_abbrev <- if(effect_label == "Odds Ratio") "aOR" 
                     else if(effect_label == "Risk Ratio") "aRR" 
                     else if(effect_label == "Hazard Ratio") "aHR"
                     else "Coef"
    
    effect_header <- paste0(effect_abbrev, " (95% CI); p-value")
    
    ## Calculate max effect string length (without expression parsing)
    effect_lengths <- nchar(data.table::fifelse(
                                            is.na(to_show_exp_clean$estimate),
                                            ref_label,
                                            paste0(to_show_exp_clean$effect_formatted, " (",
                                                   to_show_exp_clean$conf_low_formatted, "-",
                                                   to_show_exp_clean$conf_high_formatted, "); p = ",
                                                   to_show_exp_clean$p_formatted)
                                        )) + center_padding
    
    effect_width_chars <- max(effect_lengths, nchar(effect_header), na.rm = TRUE)
    columns$effect <- effect_width_chars * char_to_inch
    
    ## Calculate total table width
    table_width <- sum(unlist(columns)) + length(columns) * margin * 2
    
    ## Calculate forest width based on tbl_width proportion
    forest_width <- table_width * ((1 - tbl_width) / tbl_width)
    
    ## Convert table width to log scale units
    range_width <- diff(rangeb)
    table_width_log_units <- (table_width / forest_width) * range_width
    
    ## Calculate positions in log scale
    ## Start from the extended left edge
    rangeplot_start <- rangeb[1] - table_width_log_units
    
    ## Convert inch measurements to log scale units
    inch_to_log <- range_width / forest_width
    
    positions <- list()
    current_pos <- rangeplot_start
    
    for (name in names(columns)) {
        current_pos <- current_pos + margin * inch_to_log
        positions[[name]] <- current_pos
        current_pos <- current_pos + columns[[name]] * inch_to_log + margin * inch_to_log
    }
    
    return(list(
        table_width = table_width,
        forest_width = forest_width,
        positions = positions,
        rangeplot_start = rangeplot_start,
        total_width = table_width + forest_width,
        effect_abbrev = effect_abbrev
    ))
}
