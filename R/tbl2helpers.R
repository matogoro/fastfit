#' Add padding to exported table headers
#' @keywords internal
add_header_padding <- function (col_names) {
    padded_names <- paste0("\\rule{0pt}{3ex}", col_names, "\\rule[-1.5ex]{0pt}{0pt}")
    return(as.character(padded_names))
}

#' Add padding to exported table variables
#' @keywords internal
add_variable_padding <- function (df) {
    var_col <- if ("Variable" %in% names(df)) 
                   "Variable"
               else if ("variable" %in% names(df)) 
                   "variable"
               else NULL
    if (is.null(var_col)) 
        return(df)
    var_starts <- which(df[[var_col]] != "" & !grepl("^\\\\hspace", 
                                                     df[[var_col]]) & !grepl("^\\s+", df[[var_col]]))
    if (length(var_starts) == 0) 
        return(df)
    new_df <- data.frame()
    padding_row <- df[1, ]
    padding_row[1, ] <- ""
    new_df <- rbind(new_df, padding_row)
    for (i in seq_along(var_starts)) {
        if (i < length(var_starts)) {
            var_rows <- var_starts[i]:(var_starts[i + 1] - 1)
        }
        else {
            var_rows <- var_starts[i]:nrow(df)
        }
        new_df <- rbind(new_df, df[var_rows, ])
        padding_row <- df[1, ]
        padding_row[1, ] <- ""
        new_df <- rbind(new_df, padding_row)
    }
    return(new_df)
}

#' Check LaTeX installation
#' @keywords internal
check_latex <- function () {
    pdflatex_check <- Sys.which("pdflatex")
    if (pdflatex_check != "") 
        return(TRUE)
    xelatex_check <- Sys.which("xelatex")
    if (xelatex_check != "") 
        return(TRUE)
    return(FALSE)
}

#' Determine alignment for exported tables
#' @keywords internal
determine_alignment <- function(df) {
    align <- "r"  # Start with row numbers column (if any)
    for (col in names(df)) {
        if (col %in% c("Variable", "Group")) {
            align <- paste0(align, "l")
        } else {
            align <- paste0(align, "c")
        }
    }
    return(align)
}

#' Apply formatting to column headers in exported tables
#' @keywords internal
format_column_headers <- function (col_names, add_header_space = TRUE) {
    if (!is.character(col_names)) {
        col_names <- as.character(col_names)
    }
    col_names <- gsub("%", "\\\\%", col_names)
    col_names <- gsub("^n$", "\\\\textit{n}", col_names, ignore.case = TRUE)
    col_names <- gsub("p-value", "\\\\textit{p}-value", col_names, 
                      ignore.case = TRUE)
    col_names <- gsub("Uni p", "Uni \\\\textit{p}", col_names)
    col_names <- gsub("Multi p", "Multi \\\\textit{p}", col_names)
    if (add_header_space && !is.null(col_names)) {
        col_names <- add_header_padding(col_names)
    }
    return(as.character(col_names))
}

#' Apply formatting to column headers in exported tables (HTML-specific)
#' @keywords internal
format_column_headers_html <- function (col_names) {
    col_names <- gsub("^n$", "<i>n</i>", col_names, ignore.case = TRUE)
    col_names <- gsub("p-value", "<i>p</i>-value", col_names, 
                      ignore.case = TRUE)
    return(col_names)
}

#' Apply formatting to column headers in exported tables (HTML)
#' @keywords internal
format_indented_groups <- function (df, indent_string = "    ") {
    if (!("Variable" %in% names(df) && "Group" %in% names(df))) {
        return(df)
    }
    latex_indent <- "\\hspace{1em}"
    var_rows <- which(df$Variable != "")
    p_cols <- grep("p-value|Uni p|Multi p", names(df), value = TRUE)
    effect_cols <- c("OR (95% CI)", "HR (95% CI)", "RR (95% CI)","Estimate (95% CI)",
                     "Univariable OR (95% CI)", "Multivariable aOR (95% CI)", 
                     "Univariable HR (95% CI)", "Multivariable aHR (95% CI)",
                     "Univariable RR (95% CI)", "Multivariable aRR (95% CI)")
    is_regression_table <- any(effect_cols %in% names(df))
    is_fastfit_table <- any(c("Uni p", "Multi p") %in% names(df))
    new_df <- data.frame()
    for (i in seq_along(var_rows)) {
        current <- var_rows[i]
        next_var <- if (i < length(var_rows)) 
                        var_rows[i + 1]
                    else nrow(df) + 1
        var_name <- df$Variable[current]
        has_group <- df$Group[current] != "" && !is.na(df$Group[current]) && 
            df$Group[current] != "-"
        if (has_group) {
            if (is_regression_table) {
                var_row <- df[current, ]
                var_row$Variable <- var_name
                var_row$Group <- ""
                data_cols <- setdiff(names(df), c("Variable", 
                                                  "Group"))
                for (col in data_cols) {
                    var_row[[col]] <- ""
                }
                new_df <- rbind(new_df, var_row)
                for (j in current:(next_var - 1)) {
                    group_row <- df[j, ]
                    group_row$Variable <- paste0(latex_indent, 
                                                 df$Group[j])
                    if (is_fastfit_table) {
                        for (ec in effect_cols) {
                            if (ec %in% names(group_row)) {
                                val <- group_row[[ec]]
                                if (!is.na(val) && (val == "-" || grepl("Reference", 
                                                                        val))) {
                                    if (grepl("Univariable", ec) && "Uni p" %in% 
                                        names(group_row)) {
                                        group_row[["Uni p"]] <- ""
                                    }
                                    else if (grepl("Multivariable", ec) && 
                                             "Multi p" %in% names(group_row)) {
                                        group_row[["Multi p"]] <- ""
                                    }
                                }
                            }
                        }
                    }
                    else {
                        is_reference <- FALSE
                        for (ec in effect_cols) {
                            if (ec %in% names(group_row)) {
                                val <- group_row[[ec]]
                                if (!is.na(val) && (val == "-" || grepl("Reference", 
                                                                        val))) {
                                    is_reference <- TRUE
                                    break
                                }
                            }
                        }
                        if (is_reference) {
                            for (p_col in p_cols) {
                                if (p_col %in% names(group_row)) {
                                    group_row[[p_col]] <- ""
                                }
                            }
                        }
                    }
                    new_df <- rbind(new_df, group_row)
                }
            }
            else {
                var_row <- df[current, ]
                var_row$Variable <- var_name
                data_cols <- setdiff(names(df), c("Variable", 
                                                  "Group", p_cols))
                for (col in data_cols) {
                    var_row[[col]] <- ""
                }
                new_df <- rbind(new_df, var_row)
                for (j in current:(next_var - 1)) {
                    group_row <- df[j, ]
                    group_row$Variable <- paste0(latex_indent, 
                                                 df$Group[j])
                    for (p_col in p_cols) {
                        if (p_col %in% names(group_row)) {
                            group_row[[p_col]] <- ""
                        }
                    }
                    new_df <- rbind(new_df, group_row)
                }
            }
        }
        else {
            new_df <- rbind(new_df, df[current, ])
        }
    }
    new_df$Group <- NULL
    rownames(new_df) <- NULL
    return(new_df)
}

#' Format p-values for exported tables
#' @keywords internal
format_pvalues_export_tex <- function (df, sig_threshold = 0.05) {
    for (col in names(df)) {
        if (col == "p-value" || col == "Uni p" || col == "Multi p" || 
            grepl("p.value|pvalue", col, ignore.case = TRUE)) {
            for (i in seq_len(nrow(df))) {
                cell_value <- as.character(df[[col]][i])
                if (is.na(cell_value) || cell_value == "" || 
                    cell_value == "NA" || grepl("\\\\textbf", cell_value)) {
                    next
                }
                is_significant <- FALSE
                if (grepl("^<\\s*0\\.001", cell_value)) {
                    is_significant <- TRUE
                }
                else if (grepl("^[0-9]\\.[0-9]", cell_value) || 
                         grepl("^0\\.[0-9]", cell_value)) {
                    p_numeric <- suppressWarnings(as.numeric(gsub("[^0-9.]", 
                                                                  "", cell_value)))
                    if (!is.na(p_numeric) && p_numeric < sig_threshold) {
                        is_significant <- TRUE
                    }
                }
                if (is_significant) {
                    df[[col]][i] <- paste0("\\textbf{", cell_value, 
                                           "}")
                }
            }
        }
    }
    return(df)
}

#' Format p-values for exported tables (HTML)
#' @keywords internal
format_pvalues_export_html <- function (df, sig_threshold = 0.05) {
    for (col in names(df)) {
        if (col == "p-value" || col == "Uni p" || col == "Multi p" || 
            grepl("p.value|pvalue", col, ignore.case = TRUE)) {
            for (i in seq_len(nrow(df))) {
                cell_value <- as.character(df[[col]][i])
                if (is.na(cell_value) || cell_value == "" || 
                    cell_value == "NA") {
                    next
                }
                is_significant <- FALSE
                if (grepl("^<\\s*0\\.001", cell_value)) {
                    is_significant <- TRUE
                }
                else if (grepl("^[0-9]\\.[0-9]", cell_value) || 
                         grepl("^0\\.[0-9]", cell_value)) {
                    p_numeric <- suppressWarnings(as.numeric(gsub("[^0-9.]", 
                                                                  "", cell_value)))
                    if (!is.na(p_numeric) && p_numeric < sig_threshold) {
                        is_significant <- TRUE
                    }
                }
                if (is_significant) {
                    df[[col]][i] <- paste0("<b>", cell_value, "</b>")
                }
            }
        }
    }
    return(df)
}

#' Get paper size for PDF/LaTeX export
#' @keywords internal
get_paper_settings <- function (paper, margins = NULL) {
    paper <- match.arg(paper, c("letter", "a4", "auto"))
    settings <- switch(paper,
                       letter = list(latex_paper = "letterpaper", 
                                     width = 8.5,
                                     height = 11,
                                     default_margins = c(1, 1, 1, 1)),
                       a4 = list(latex_paper = "a4paper",
                                 width = 8.27, 
                                 height = 11.69,
                                 default_margins = c(1, 1, 1, 1)),
                       auto = list(latex_paper = "letterpaper", 
                                   width = NULL, height = NULL, default_margins = c(0.5, 0.5, 0.5, 0.5)))
    if (!is.null(margins)) {
        if (length(margins) == 1) {
            margins <- rep(margins, 4)
        }
        else if (length(margins) != 4) {
            stop("margins must be length 1 or 4")
        }
        settings$margins <- margins
    }
    else {
        settings$margins <- settings$default_margins
    }
    return(settings)
}

#' Sanitize certain symbols for LaTeX
#' @keywords internal
sanitize_for_latex <- function (x) {
    if (is.null(x) || length(x) == 0) 
        return(x)
    has_latex <- grepl("\\\\text(bf|it)\\{", x)
    needs_sanitizing <- !has_latex & !is.na(x)
    if (any(needs_sanitizing)) {
        x[needs_sanitizing] <- gsub("%", "\\\\%", x[needs_sanitizing])
        x[needs_sanitizing] <- gsub("&", "\\\\&", x[needs_sanitizing])
        x[needs_sanitizing] <- gsub("#", "\\\\#", x[needs_sanitizing])
        x[needs_sanitizing] <- gsub("_", "\\\\_", x[needs_sanitizing])
        x[needs_sanitizing] <- gsub("\\$", "\\\\$", x[needs_sanitizing])
    }
    return(x)
}

#' Format column headers with n counts (TeX)
#' @keywords internal
format_column_headers_with_n_tex <- function(col_names, n_row_data) {
    new_names <- character(length(col_names))
    
    for (i in seq_along(col_names)) {
        col <- col_names[i]
        
        ## Check if this column should have n count
        has_n <- col %in% names(n_row_data) && 
            !is.na(n_row_data[[col]]) && 
            n_row_data[[col]] != "" && 
            n_row_data[[col]] != "0"
        
        if (col %in% c("Variable", "Group")) {
            ## Standard formatting for these columns
            new_names[i] <- format_column_headers(col)
        } else if (col %in% c("p-value", "p value")) {
            ## Format p-value with italics
            new_names[i] <- format_column_headers(col)
        } else if (has_n) {
            ## Column with n count - use tabular for text mode
            n_value <- n_row_data[[col]]
            ## Add vertical padding with struts
            new_names[i] <- paste0("\\begin{tabular}{@{}c@{}}\\rule{0pt}{2.5ex}", 
                                   col, 
                                   "\\\\[-0ex] (\\textit{n} = ", 
                                   format(as.numeric(gsub(",", "", n_value)), big.mark = ","), 
                                   ")\\rule[-1ex]{0pt}{0pt}\\end{tabular}")
        } else {
            ## Regular column without n count
            new_names[i] <- format_column_headers(col)
        }
    }
    
    return(new_names)
}

#' Format column headers with n counts (HTML)
#' @keywords internal
format_column_headers_with_n_html <- function(col_names, n_row_data) {
    new_names <- col_names
    for (i in seq_along(col_names)) {
        col <- col_names[i]
        
        ## Skip Variable, Group, and p-value columns
        if (col %in% c("Variable", "Group", "p-value", "p value")) {
            new_names[i] <- format_column_headers_html(col)
            next
        }
        
        ## Get n value from the N row
        if (col %in% names(n_row_data)) {
            n_value <- n_row_data[[col]]
            if (!is.na(n_value) && n_value != "" && n_value != "0") {
                ## Format header with only n italicized
                clean_col <- format_column_headers_html(col)
                new_names[i] <- paste0(clean_col, "<br>(<i>n</i> = ", n_value, ")")
            } else {
                new_names[i] <- format_column_headers_html(col)
            }
        } else {
            new_names[i] <- format_column_headers_html(col)
        }
    }
    return(new_names)
}
