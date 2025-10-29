#' Export Table to LaTeX Format
#'
#' Converts a data frame or data.table to LaTeX format suitable for inclusion
#' in LaTeX documents. Supports various table packages and formatting options.
#'
#' @param table A data.frame, data.table, or matrix to export.
#' @param file Character string specifying the output .tex filename. Must have
#'   .tex extension.
#' @param caption Character string. Table caption for \\caption{} command.
#' @param format_headers Logical. If TRUE, formats column headers by converting
#'   underscores to spaces and applying title case. Default is TRUE.
#' @param var_padding Logical. If TRUE, adds vertical spacing around variable
#'   groups for better readability. Default is TRUE.
#' @param cell_padding Character string or numeric. Add vertical padding to text
#'   within cells. Choices are "none", "normal" (default), "relaxed", "loose", or
#'   a number representing a custom multiplier.
#' @param bold_significant Logical. If TRUE, wraps significant p-values in
#'   \\textbf{} commands. Default is TRUE.
#' @param sig_threshold Numeric. P-value threshold for bolding (e.g., 0.05).
#'   Values below this are considered significant. Default is 0.05.
#' @param align Character string or vector specifying column alignment. If NULL,
#'   automatically determines based on column content. Options: "l" (left),
#'   "c" (center), "r" (right), or "p{width}" for paragraph columns.
#' @param indent_groups Logical. If TRUE, uses \\hspace to indent grouped rows
#'   creating a hierarchical display. Default is FALSE.
#' @param condense_table Logical. If TRUE, decreases table vertical height by
#'   reducing continuous, survival, and binary categorical variables to single
#'   rows. Makes indent_groups = TRUE by default. Default is FALSE.
#' @param booktabs Logical. If TRUE, uses booktabs package commands for
#'   professional-quality rules. Default is FALSE.
#' @param zebra_stripes Logical. If TRUE, adds alternating row colors for variable
#'   groups. Requires \\usepackage[table]{xcolor} in preamble. Default is FALSE.
#' @param stripe_color Character string. LaTeX color specification for zebra stripes
#'   (e.g., "gray!20", "blue!10"). Default is "gray!20".
#' @param dark_header Logical. If TRUE, creates white text on black background for
#'   header row. Requires \\usepackage[table]{xcolor} in preamble. Default is FALSE.
#' @param label Character string. LaTeX label for cross-references using \\ref{}.
#' @param ... Additional arguments passed to xtable::xtable.
#'
#' @return Invisibly returns NULL. Creates a .tex file at the specified location.
#'
#' @details
#' The function generates a LaTeX tabular environment that can be directly
#' included in LaTeX documents using \\input{} or \\include{}. Special characters
#' are automatically escaped for LaTeX compatibility.
#'
#' A comprehensive list of LaTeX packages needed for full functionality is as
#' follows (installed by default in most LaTeX distributions, with some
#' exceptions):
#' \itemize{
#'   \item \code{fontenc}
#'   \item \code{inputenc}
#'   \item \code{array}
#'   \item \code{booktabs}
#'   \item \code{longtable}
#'   \item \code{graphicx}
#'   \item \code{geometry}
#'   \item \code{pdflscape}
#'   \item \code{lscape}
#'   \item \code{helvet}
#'   \item \code{standalone}
#'   \item \code{varwidth}
#'   \item \code{float}
#'   \item \code{caption}
#'   \item \code{xcolor} (required for zebra_stripes or dark_header)
#'   \item \code{colortbl} (required for zebra_stripes or dark_header)
#' }
#' When \code{booktabs = TRUE}, uses \code{\\toprule}, \code{\\midrule}, and
#'     \code{\\bottomrule} for cleaner table appearance. This requires
#'     \code{\\usepackage{booktabs}} in the LaTeX preamble.
#' 
#' When \code{dark_header = TRUE} or \code{zebra_stripes = TRUE}, requires
#'     \code{\\usepackage[table]{xcolor}} in the LaTeX preamble.
#' 
#' The indent_groups option is particularly useful for regression tables with
#'     categorical variables, creating a clear visual hierarchy.
#' 
#' This function is tailored to outputs from \code{\link{desctbl()}},
#' \code{\link{fastfit()}}, \code{\link{uscreen()}}, \code{\link{fit()}},
#' and \code{\link{compfit()}}, although it can theoretically be applied to
#' any data frame or data.table object.
#' 
#' @examples
#' if (FALSE) {
#' # Basic LaTeX export
#' tbl2tex(results, "results.tex")
#' 
#' # Publication-quality table with booktabs
#' tbl2tex(final_results, "publication.tex",
#'         booktabs = TRUE,
#'         caption = "Multivariable regression analysis",
#'         label = "tab:regression")
#' 
#' # Multi-line caption
#' tbl2tex(table, "results.tex",
#'         caption = "Table 1 - Regression results\\\\
#'                   aHR = adjusted hazard ratio\\\\
#'                   \\textsuperscript{1}Note goes here")
#' 
#' # Hierarchical table with indentation
#' tbl2tex(model_output, "model.tex",
#'         indent_groups = TRUE,
#'         bold_significant = TRUE)
#' 
#' # Custom alignment
#' tbl2tex(wide_table, "aligned.tex",
#'         align = c("l", "l", "r", "r", "c"))
#'         
#' # Table with dark header
#' tbl2tex(results, "dark_header.tex",
#'         dark_header = TRUE,
#'         bold_significant = TRUE)
#'         
#' # Table with zebra stripes
#' tbl2tex(results, "striped.tex",
#'         zebra_stripes = TRUE,
#'         stripe_color = "blue!10")
#' }
#' @seealso
#' \code{\link{tbl2pdf}} for PDF output,
#' \code{\link{tbl2html}} for HTML output
#'
#' @export
tbl2tex <- function (table,
                     file,
                     format_headers = TRUE,
                     var_padding = TRUE, 
                     cell_padding = "normal",
                     bold_significant = TRUE,
                     sig_threshold = 0.05,
                     align = NULL, 
                     indent_groups = FALSE,
                     condense_table = FALSE,
                     booktabs = FALSE,
                     zebra_stripes = FALSE,
                     stripe_color = "gray!20",
                     dark_header = FALSE,
                     caption = NULL,
                     label = NULL, 
                     ...) {
    
    if (!requireNamespace("xtable", quietly = TRUE)) {
        stop("Package 'xtable' required. Install with: install.packages('xtable')")
    }
    
    if (!grepl("\\.tex$", tolower(file))) {
        stop("File must have .tex extension")
    }
    
    df <- as.data.frame(table)

                                        # Detect variable groups BEFORE any processing (for zebra stripes)
    var_groups <- NULL
    if (zebra_stripes) {
        if ("Variable" %in% names(df)) {
            var_starts <- which(df$Variable != "" & !is.na(df$Variable))
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
            }
        }
    }

    has_n_row <- FALSE
    n_row_data <- NULL
    if (nrow(df) > 0 && "Variable" %in% names(df) && df$Variable[1] == "N") {
        has_n_row <- TRUE
        n_row_data <- df[1, ]
        df <- df[-1, ]  # Remove N row from data
                                        # Adjust var_groups if N row was removed
        if (!is.null(var_groups) && length(var_groups) > 1) {
            var_groups <- var_groups[-1]
        }
    }

    if (condense_table) {
        indent_groups <- TRUE
        df <- condense_table_rows(df, indent_groups = indent_groups)
                                        # Re-detect groups after condensing
        if (zebra_stripes && "Variable" %in% names(df)) {
            var_starts <- which(df$Variable != "" & !is.na(df$Variable))
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
            }
        }
        df <- format_indented_groups(df, indent_string = "\\hspace{1em}")
    } else if (indent_groups) {
        df <- format_indented_groups(df, indent_string = "\\hspace{1em}")
    }

    if (bold_significant) {
        df <- format_pvalues_export_tex(df, sig_threshold)
    }
    
                                        # Track padding rows if adding padding
    if (var_padding && ("Variable" %in% names(df) || "variable" %in% names(df))) {
                                        # Before adding padding, adjust var_groups if necessary
        if (!is.null(var_groups)) {
                                        # Find where padding rows will be inserted
            padding_positions <- which(diff(var_groups) != 0)
                                        # Adjust var_groups for padding
            new_var_groups <- integer(nrow(df) + length(padding_positions))
            current_pos <- 1
            for (i in 1:nrow(df)) {
                new_var_groups[current_pos] <- var_groups[i]
                current_pos <- current_pos + 1
                if (i %in% padding_positions) {
                                        # Padding row gets group 0 (no shading)
                    new_var_groups[current_pos] <- 0
                    current_pos <- current_pos + 1
                }
            }
            var_groups <- new_var_groups
        }
        df <- add_variable_padding(df)
    }
    
                                        # Set up add.to.row for zebra stripes
    add.to.row <- NULL
    if (zebra_stripes && "Variable" %in% names(df)) {
                                        # Check if table has been indented
        is_indented <- indent_groups || condense_table
        
        if (is_indented) {
                                        # For indented tables - look for rows where Variable is not empty and not indented
            var_starts <- which(!grepl("\\\\hspace", df$Variable) & 
                                df$Variable != "" & 
                                !is.na(df$Variable))
            
            commands <- character()
            positions <- numeric()
            
            for (i in seq_along(var_starts)) {
                if (i %% 2 != 0) {  # Odd variable groups get shading
                    start_idx <- var_starts[i]
                    end_idx <- if (i < length(var_starts)) {
                                   if (var_padding) {
                                       var_starts[i + 1] - 2
                                   } else {
                                       var_starts[i + 1] - 1
                                   }
                               } else {
                                   nrow(df)
                               }
                    
                    for (row in start_idx:end_idx) {
                        if (!is.na(df$Variable[row])) {
                            commands <- c(commands, paste0("\\rowcolor{", stripe_color, "}"))
                            positions <- c(positions, row - 1)
                        }
                    }
                }
            }
            
            if (length(commands) > 0) {
                add.to.row <- list(pos = as.list(positions), command = commands)
            }
        } else if ("Group" %in% names(df)) {
                                        # For non-indented tables with Group column
            var_starts <- which(df$Variable != "" & !is.na(df$Variable))
            
            commands <- character()
            positions <- numeric()
            
            for (i in seq_along(var_starts)) {
                if (i %% 2 != 0) {
                    start_idx <- var_starts[i]
                    end_idx <- if (i < length(var_starts)) {
                                   if (var_padding) {
                                       var_starts[i + 1] - 2
                                   } else {
                                       var_starts[i + 1] - 1
                                   }
                               } else {
                                   nrow(df)
                               }
                    
                    for (row in start_idx:end_idx) {
                        if (!is.na(df$Variable[row])) {
                            commands <- c(commands, paste0("\\rowcolor{", stripe_color, "}"))
                            positions <- c(positions, row - 1)
                        }
                    }
                }
            }
            
            if (length(commands) > 0) {
                add.to.row <- list(pos = as.list(positions), command = commands)
            }
        } else if (!is.null(var_groups)) {
                                        # Use var_groups for simpler tables
            commands <- character()
            positions <- numeric()
            
            for (i in 1:nrow(df)) {
                if (var_groups[i] %% 2 != 0 && var_groups[i] > 0) {
                    commands <- c(commands, paste0("\\rowcolor{", stripe_color, "}"))
                    positions <- c(positions, i - 1)
                }
            }
            
            if (length(commands) > 0) {
                add.to.row <- list(pos = as.list(positions), command = commands)
            }
        }
    }
    
                                        # Add cell padding command if requested
                                        # Determine array stretch value for cell padding
    array_stretch <- switch(as.character(cell_padding),
                            "none" = NULL,
                            "normal" = "1.3",
                            "relaxed" = "1.5",
                            "loose" = "1.8",
                            {
                                        # If numeric value provided
                                if (is.numeric(cell_padding)) {
                                    as.character(cell_padding)
                                } else {
                                    "1.3"  # Default to normal
                                }
                            })
    
                                        # Add arraystretch command at the beginning of the table
    if (!is.null(array_stretch)) {
        arraystretch_command <- sprintf("\\renewcommand{\\arraystretch}{%s}", array_stretch)
        
        if (!is.null(add.to.row) && length(add.to.row$pos) > 0) {
                                        # Check if position -1 already exists
            if (-1 %in% unlist(add.to.row$pos)) {
                                        # Find the index and prepend to that command
                idx <- which(unlist(add.to.row$pos) == -1)
                add.to.row$command[idx] <- paste0(arraystretch_command, " ", add.to.row$command[idx])
            } else {
                                        # Add as first command
                add.to.row$pos <- c(list(-1), add.to.row$pos)
                add.to.row$command <- c(arraystretch_command, add.to.row$command)
            }
        } else {
                                        # Create new add.to.row with arraystretch
            add.to.row <- list(
                pos = list(-1),
                command = arraystretch_command
            )
        }
    }
    
    if (is.null(align)) {
        align <- determine_alignment(df)
    }
    
                                        # Store original column names before formatting
    original_col_names <- names(df)
    
    if (format_headers) {
        if (has_n_row) {
            names(df) <- format_column_headers_with_n_tex(names(df), n_row_data)
        } else {
            names(df) <- format_column_headers(names(df))
        }
    }
    
                                        # PROPERLY IMPLEMENTED DARK HEADER
    if (dark_header) {
                                        # Wrap each column name with the white color command
        col_names <- names(df)
        for (i in seq_along(col_names)) {
            col_names[i] <- paste0("\\color{white}", col_names[i])
        }
        names(df) <- col_names
        
                                        # Add the rowcolor command (without the color command)
        header_command <- "\\rowcolor{black}"
        
        if (!is.null(add.to.row) && length(add.to.row$pos) > 0) {
                                        # Check if we need to combine with arraystretch at position -1
            if (-1 %in% unlist(add.to.row$pos)) {
                                        # Find the index and append to that command
                idx <- which(unlist(add.to.row$pos) == -1)
                add.to.row$command[idx] <- paste0(add.to.row$command[idx], " ", header_command)
            } else {
                                        # Prepend header command
                add.to.row$pos <- c(list(-1), add.to.row$pos)
                add.to.row$command <- c(header_command, add.to.row$command)
            }
        } else {
                                        # Create new add.to.row
            add.to.row <- list(
                pos = list(-1),
                command = header_command
            )
        }
    }
    
    xt <- xtable::xtable(df,
                         align = align,
                         caption = caption, 
                         label = label, ...)
    
                                        # Prepare file output with necessary package declarations if needed
    file_conn <- file(file, "w")
    
                                        # Write comment about required packages if using special features
    if (dark_header || zebra_stripes) {
        writeLines("% This table requires \\usepackage[table]{xcolor} in your LaTeX preamble", file_conn)
    }
    if (!is.null(array_stretch)) {
        writeLines(sprintf("%% Cell padding applied with \\arraystretch{%s}", array_stretch), file_conn)
    }
    
    close(file_conn)
    
                                        # Print with or without add.to.row commands
    if (!is.null(add.to.row) && length(add.to.row$pos) > 0) {
        print(xt,
              file = file,
              append = TRUE,  # Append after our comments
              include.rownames = FALSE,
              booktabs = booktabs, 
              floating = FALSE,
              add.to.row = add.to.row,
              sanitize.text.function = sanitize_for_latex, 
              sanitize.rownames.function = sanitize_for_latex,
              sanitize.colnames.function = function(x) x,  # Don't sanitize since we added color commands
              hline.after = c(-1, 0, nrow(xt)))
    } else {
        print(xt,
              file = file,
              append = TRUE,  # Append after our comments
              include.rownames = FALSE,
              booktabs = booktabs, 
              floating = FALSE,
              sanitize.text.function = sanitize_for_latex, 
              sanitize.rownames.function = sanitize_for_latex,
              sanitize.colnames.function = function(x) x,  # Don't sanitize since we may have formatted headers
              hline.after = c(-1, 0, nrow(xt)))
    }
    
                                        # Add a note about required packages if using special features
    if (dark_header || zebra_stripes) {
        message("Note: This table requires \\usepackage[table]{xcolor} in your LaTeX preamble")
    }
    if (!is.null(array_stretch) && cell_padding != "none") {
        message(sprintf("Note: Cell padding applied with \\arraystretch{%s}", array_stretch))
    }
    
    message(sprintf("Table exported to %s", file))
    
    invisible(NULL)
}
