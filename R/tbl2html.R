#' Export Table to HTML Format
#'
#' Converts a data frame or data.table to HTML format with optional CSS styling
#' for web display or inclusion in HTML documents.
#'
#' @param table A data.frame, data.table, or matrix to export.
#' @param file Character string specifying the output HTML filename. Must have
#'   .html or .htm extension.
#' @param caption Character string. Optional caption to display below table.
#' @param format_headers Logical. If TRUE, formats column headers for better
#'   display (converts underscores to spaces, applies proper casing). Default
#'   is TRUE.
#' @param variable_padding Logical. If TRUE, adds padding around variable groups
#'   for improved readability. Default is FALSE.
#' @param bold_significant Logical. If TRUE, wraps significant p-values in
#'   HTML bold tags. Default is TRUE.
#' @param sig_threshold Numeric. P-value threshold for significance highlighting.
#'   Default is 0.05.
#' @param indent_groups Logical. If TRUE, indents grouped rows using non-breaking
#'   spaces for hierarchical display. Default is FALSE.
#' @param condense_table Logical. If TRUE, decreases table vertical height by
#'   reducing continuous, survival, and binary categorical variables to single
#'   rows. Makes indent_groups = TRUE by default. Default is FALSE.
#' @param include_css Logical. If TRUE, includes basic CSS styling in the output
#'   file. Set to FALSE if including in existing HTML with its own styles.
#'   Default is TRUE.
#' @param ... Additional arguments passed to xtable::xtable.
#'
#' @return Invisibly returns NULL. Creates an HTML file at the specified location.
#'
#' @details
#' The function generates an HTML table that can be viewed directly in web
#' browsers or embedded in HTML documents. When \code{include_css = TRUE}, adds
#' minimal styling for borders, padding, and header formatting.
#' 
#' The indent_groups option uses HTML non-breaking spaces (&nbsp;) to create
#' visual hierarchy, useful for categorical variables in regression output.
#' 
#' P-values can be automatically highlighted using bold formatting when they
#' fall below the significance threshold, making important results easy to identify.
#'
#' This function is tailored to outputs from \code{\link{desctbl()}},
#' \code{\link{fastfit()}}, \code{\link{uscreen()}}, \code{\link{fit()}},
#' and \code{\link{compfit()}}, although it can theoretically be applied to
#' any data frame or data.table object.
#'
#' @examples
#' if (FALSE) {
#' # Basic HTML export
#' tbl2html(results, "results.html")
#' 
#' # Without CSS for embedding
#' tbl2html(table_data, "fragment.html",
#'          include_css = FALSE)
#' 
#' # Hierarchical display with significance highlighting
#' tbl2html(regression_output, "regression.html",
#'          indent_groups = TRUE,
#'          bold_significant = TRUE,
#'          sig_threshold = 0.01)
#' 
#' # Clean headers without special formatting
#' tbl2html(summary_table, "summary.html",
#'          format_headers = FALSE)
#' }
#' @seealso
#' \code{\link{tbl2pdf}} for PDF output, 
#' \code{\link{tbl2tex}} for LaTeX output
#' 
#' @export
tbl2html <- function(table,
                     file,
                     caption = NULL,
                     format_headers = TRUE,
                     variable_padding = FALSE, 
                     bold_significant = TRUE,
                     sig_threshold = 0.05,
                     indent_groups = FALSE,
                     condense_table = FALSE,
                     zebra_stripes = FALSE,
                     stripe_color = "#EEEEEE",
                     dark_header = FALSE,
                     include_css = TRUE,
                     ...) {
    
    if (!requireNamespace("xtable", quietly = TRUE)) {
        stop("Package 'xtable' required. Install with: install.packages('xtable')")
    }
    
    if (!grepl("\\.(html|htm)$", tolower(file))) {
        stop("File must have .html or .htm extension")
    }
    
    df <- as.data.frame(table)

    has_n_row <- FALSE
    n_row_data <- NULL
    if (nrow(df) > 0 && "Variable" %in% names(df) && df$Variable[1] == "N") {
        has_n_row <- TRUE
        n_row_data <- df[1, ]
        df <- df[-1, ]
    }

                                        # Apply transformations first
    if (variable_padding && !condense_table && ("Variable" %in% names(df) || "variable" %in% names(df))) {
        df <- add_variable_padding(df)
    }

    if (condense_table) {
        indent_groups <- TRUE
        df <- condense_table_rows(df, indent_groups = indent_groups)
        df <- format_indented_groups(df, indent_string = "&nbsp;&nbsp;&nbsp;&nbsp;")
    } else if (indent_groups) {
        df <- format_indented_groups(df, indent_string = "&nbsp;&nbsp;&nbsp;&nbsp;")
    }

    if (bold_significant) {
        df <- format_pvalues_export_html(df, sig_threshold)
    }
    
                                        # NOW calculate var_groups on the FINAL dataframe structure
    var_groups <- NULL
    if (zebra_stripes && "Variable" %in% names(df)) {
                                        # For indented tables, non-indented rows mark variable starts
        if (indent_groups || condense_table) {
            var_starts <- which(!grepl("^&nbsp;", df$Variable) & 
                                df$Variable != "" & 
                                !is.na(df$Variable))
        } else {
                                        # For regular tables, any non-empty Variable marks a start
            var_starts <- which(df$Variable != "" & !is.na(df$Variable))
        }
        
        if (length(var_starts) > 0) {
            var_groups <- integer(nrow(df))
            for (i in seq_along(var_starts)) {
                start_idx <- var_starts[i]
                end_idx <- if (i < length(var_starts)) {
                               var_starts[i + 1] - 1
                           } else {
                               nrow(df)
                           }
                var_groups[start_idx:end_idx] <- i
            }
            
                                        # Apply zebra stripes by tracking groups
            df$.zebra_group <- var_groups
        }
    }
    
    if (format_headers) {
        if (has_n_row) {
            names(df) <- format_column_headers_with_n_html(names(df), n_row_data)
        } else {
            names(df) <- format_column_headers_html(names(df))
        }
    }
    
                                        # Remove zebra tracking column from display
    display_df <- df
    if (zebra_stripes && ".zebra_group" %in% names(display_df)) {
        zebra_groups <- display_df$.zebra_group
        display_df$.zebra_group <- NULL
    } else {
        zebra_groups <- NULL
    }
    
                                        # Remove debug print
                                        # print(zebra_groups)
    
    xt <- xtable::xtable(display_df, caption = caption, ...)
    
                                        # Build CSS (unchanged)
    if (include_css) {
        css <- paste0("<style>\n",
                      "table { \n",
                      "border-collapse: collapse; \n",
                      "font-family: Arial, sans-serif;\n",
                      "margin: 20px;\n",
                      "}\n",
                      "th, td { \n",
                      "padding: 8px 12px; \n",
                      "text-align: left; \n",
                      "border: 1px solid #ddd;\n",
                      "}\n")
        
        if (indent_groups) {
            css <- paste0(css,
                          "th:not(:first-child), \n",
                          "td:not(:first-child) { \n",
                          "text-align: center; \n",
                          "}\n")
        } else {
            css <- paste0(css,
                          "th:not(:nth-child(1)):not(:nth-child(2)), \n",
                          "td:not(:nth-child(1)):not(:nth-child(2)) { \n",
                          "text-align: center; \n",
                          "}\n")
        }
        
                                        # Header styling - dark or light
        if (dark_header) {
            css <- paste0(css,
                          "th { \n",
                          "background-color: #000000; \n",
                          "color: #FFFFFF; \n",
                          "font-weight: bold;\n",
                          "}\n")
        } else {
            css <- paste0(css,
                          "th { \n",
                          "background-color: #f2f2f2; \n",
                          "font-weight: bold;\n",
                          "}\n")
        }
        
                                        # Zebra stripe styling
        if (zebra_stripes) {
            css <- paste0(css,
                          "tr.zebra-stripe { \n",
                          "background-color: ", stripe_color, ";\n",
                          "}\n")
        }
        
        css <- paste0(css,
                      "caption {\n",
                      "text-align: left;\n",
                      "margin-top: 10px;\n",
                      "margin-bottom: 10px;\n",
                      "font-weight: bold;\n",
                      "font-size: 1.1em;\n",
                      "}\n",
                      "</style>\n")
    }
    
                                        # Generate HTML table
    if (zebra_stripes && !is.null(zebra_groups)) {
                                        # Capture the HTML output
        html_output <- capture.output(
            print(xt,
                  type = "html",
                  include.rownames = FALSE,
                  sanitize.text.function = identity,
                  sanitize.rownames.function = identity,
                  sanitize.colnames.function = identity,
                  ...)
        )
        
                                        # Add classes to TR elements based on zebra groups
        tr_count <- 0
        for (i in seq_along(html_output)) {
            if (grepl("^  <tr>", html_output[i])) {
                tr_count <- tr_count + 1
                if (tr_count <= length(zebra_groups)) {
                    group_num <- zebra_groups[tr_count]
                    if (!is.na(group_num) && group_num %% 2 == 1) {
                        html_output[i] <- gsub("  <tr>", '  <tr class="zebra-stripe">', html_output[i])
                    }
                }
            }
        }
        
                                        # Write everything to file
        if (include_css) {
            writeLines(c(css, html_output), file)
        } else {
            writeLines(html_output, file)
        }
    } else {
                                        # Standard output without zebra stripes
        if (include_css) {
            cat(css, file = file)
            print(xt,
                  type = "html",
                  file = file,
                  append = TRUE,
                  include.rownames = FALSE,
                  sanitize.text.function = identity,
                  sanitize.rownames.function = identity,
                  sanitize.colnames.function = identity,
                  ...)
        } else {
            print(xt,
                  type = "html",
                  file = file,
                  append = FALSE,
                  include.rownames = FALSE,
                  sanitize.text.function = identity,
                  sanitize.rownames.function = identity,
                  sanitize.colnames.function = identity,
                  ...)
        }
    }
    
    message(sprintf("Table exported to %s", file))
    
    invisible(NULL)
}
