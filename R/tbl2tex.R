#' Export Table to LaTeX Format
#'
#' Converts a data frame, data.table, or matrix to LaTeX source code suitable for 
#' inclusion in LaTeX documents. Generates publication-quality table markup with 
#' extensive formatting options including booktabs styling, color schemes, and 
#' hierarchical displays. Output can be directly \code{\\input{}} or \code{\\include{}} 
#' into LaTeX manuscripts.
#'
#' @param table A data.frame, data.table, or matrix to export. Can be output from 
#'   \code{\link{desctbl}}, \code{\link{fit}}, \code{\link{uscreen}}, 
#'   \code{\link{fastfit}}, \code{\link{compfit}}, or any tabular data.
#'   
#' @param file Character string specifying the output .tex filename. Must have 
#'   \code{.tex} extension. Example: \code{"results.tex"}, \code{"table1.tex"}.
#'   
#' @param caption Character string. Table caption for \code{\\caption\{\}} command. 
#'   Supports multi-line captions using \code{\\\\}. Default is \code{NULL}.
#'   
#' @param format_headers Logical. If \code{TRUE}, formats column headers by 
#'   converting underscores to spaces, italicizing statistical notation (\emph{n}, 
#'   \emph{p}), and applying title case. Default is \code{TRUE}.
#'   
#' @param var_padding Logical. If \code{TRUE}, adds vertical spacing around 
#'   variable groups using \code{\\addlinespace} for improved readability. 
#'   Default is \code{TRUE}.
#'   
#' @param cell_padding Character string or numeric. Vertical padding within cells:
#'   \itemize{
#'     \item \code{"none"} - No extra padding
#'     \item \code{"normal"} - Standard padding [default]
#'     \item \code{"relaxed"} - Increased padding
#'     \item \code{"loose"} - Maximum padding
#'     \item Numeric - Custom \code{\\arraystretch} value
#'   }
#'   
#' @param bold_significant Logical. If \code{TRUE}, wraps significant p-values 
#'   in \code{\\textbf\{\}} commands for bold display. Default is \code{TRUE}.
#'   
#' @param sig_threshold Numeric. P-value threshold for bolding. Default is 0.05.
#'   
#' @param align Character string or vector specifying column alignment:
#'   \itemize{
#'     \item \code{"l"} - Left
#'     \item \code{"c"} - Center  
#'     \item \code{"r"} - Right
#'     \item \code{"p\{width\}"} - Paragraph column with specified width
#'   }
#'   If \code{NULL}, automatically determines based on content. Can specify 
#'   per-column as vector. Default is \code{NULL}.
#'   
#' @param indent_groups Logical. If \code{TRUE}, uses \code{\\hspace\{\}} to 
#'   indent grouped rows, creating hierarchical display. Useful for factor 
#'   variables in regression tables. Default is \code{FALSE}.
#'   
#' @param condense_table Logical. If \code{TRUE}, condenses table by showing 
#'   only essential rows (single row for continuous, non-reference for binary). 
#'   Automatically sets \code{indent_groups = TRUE}. Default is \code{FALSE}.
#'   
#' @param booktabs Logical. If \code{TRUE}, uses booktabs package commands 
#'   (\code{\\toprule}, \code{\\midrule}, \code{\\bottomrule}) for professional 
#'   table rules. Requires \code{\\usepackage\{booktabs\}} in LaTeX preamble. 
#'   Default is \code{FALSE}.
#'   
#' @param zebra_stripes Logical. If \code{TRUE}, adds alternating row colors 
#'   for variable groups using \code{\\rowcolor\{\}}. Requires 
#'   \code{\\usepackage[table]\{xcolor\}} in preamble. Default is \code{FALSE}.
#'   
#' @param stripe_color Character string. LaTeX color specification for zebra 
#'   stripes (e.g., \code{"gray!20"}, \code{"blue!10"}). Only used when 
#'   \code{zebra_stripes = TRUE}. Default is \code{"gray!20"}.
#'   
#' @param dark_header Logical. If \code{TRUE}, creates white text on black 
#'   background for header row using \code{\\rowcolor\{black\}} and 
#'   \code{\\color\{white\}}. Requires \code{\\usepackage[table]\{xcolor\}}. 
#'   Default is \code{FALSE}.
#'   
#' @param label Character string. LaTeX label for cross-references using 
#'   \code{\\label\{\}} and \code{\\ref\{\}}. Example: \code{"tab:regression"}. 
#'   Default is \code{NULL}.
#'   
#' @param ... Additional arguments passed to \code{\link[xtable]{xtable}}.
#'
#' @return Invisibly returns \code{NULL}. Creates a .tex file at the specified 
#'   location containing a LaTeX tabular environment.
#'
#' @details
#' \strong{Output Format:}
#' 
#' The function generates a standalone LaTeX tabular environment that can be:
#' \enumerate{
#'   \item Included in documents: \code{\\input\{results.tex\}}
#'   \item Embedded in table/figure environments
#'   \item Used in manuscript classes (article, report, etc.)
#' }
#' 
#' The output includes:
#' \itemize{
#'   \item Complete tabular environment with proper alignment
#'   \item Horizontal rules (\code{\\hline} or booktabs rules)
#'   \item Column headers with optional formatting
#'   \item Data rows with automatic escaping of special characters
#'   \item Optional caption and label commands
#' }
#' 
#' \strong{Required LaTeX Packages:}
#' 
#' Add these to your LaTeX document preamble:
#' 
#' \emph{Always required:}
#' \preformatted{
#' \\usepackage[T1]{fontenc}
#' \\usepackage[utf8]{inputenc}
#' \\usepackage{array}
#' \\usepackage{graphicx}  % If using resizebox
#' }
#' 
#' \emph{Optional (based on parameters):}
#' \preformatted{
#' \\usepackage{booktabs}  % For booktabs = TRUE
#' \\usepackage[table]{xcolor}  % For zebra_stripes or dark_header
#' }
#' 
#' \strong{Booktabs Style:}
#' 
#' When \code{booktabs = TRUE}, the table uses publication-quality rules:
#' \itemize{
#'   \item \code{\\toprule} - Heavy rule at top
#'   \item \code{\\midrule} - Medium rule below headers
#'   \item \code{\\bottomrule} - Heavy rule at bottom
#'   \item No vertical rules (booktabs philosophy)
#'   \item Better spacing around rules
#' }
#' 
#' This is the preferred style for most academic journals.
#' 
#' \strong{Color Features:}
#' 
#' \emph{Zebra Stripes:}
#' Creates alternating background colors for visual grouping:
#' \preformatted{
#' zebra_stripes = TRUE
#' stripe_color = "gray!20"  # 20\% gray
#' stripe_color = "blue!10"  # 10\% blue  
#' }
#' 
#' \emph{Dark Header:}
#' Creates high-contrast header row:
#' \preformatted{
#' dark_header = TRUE  # Black background, white text
#' }
#' 
#' Both require \code{\\usepackage[table]{xcolor}} in your document.
#' 
#' \strong{Integration with LaTeX Documents:}
#' 
#' \emph{Basic inclusion:}
#' \preformatted{
#' \\begin{table}[htbp]
#'   \\centering
#'   \\caption{Regression Results}
#'   \\label{tab:regression}
#'   \\input{results.tex}
#' \\end{table}
#' }
#' 
#' \emph{With resizing:}
#' \preformatted{
#' \\begin{table}[htbp]
#'   \\centering
#'   \\caption{Results}
#'   \\resizebox{\\textwidth}{!}{\\input{results.tex}}
#' \\end{table}
#' }
#' 
#' \emph{Landscape orientation:}
#' \preformatted{
#' \\usepackage{pdflscape}
#' \\begin{landscape}
#'   \\begin{table}[htbp]
#'     \\centering
#'     \\input{wide_results.tex}
#'   \\end{table}
#' \\end{landscape}
#' }
#' 
#' \strong{Caption Formatting:}
#' 
#' Captions in the \code{caption} parameter are written as LaTeX comments in 
#' the output file for reference. For actual LaTeX captions, wrap the table 
#' in a table environment (see examples above).
#' 
#' \strong{Special Characters:}
#' 
#' The function automatically escapes LaTeX special characters in your data:
#' \itemize{
#'   \item Ampersand becomes \\&
#'   \item Percent becomes \\%
#'   \item Dollar becomes \\$
#'   \item Hash becomes \\#
#'   \item Underscore becomes \\_
#'   \item Left brace becomes \\{
#'   \item Right brace becomes \\}
#'   \item Tilde becomes \\textasciitilde{}
#'   \item Caret becomes \\textasciicircum{}
#' }
#' 
#' Variable names and labels should not include these characters unless 
#' intentionally using LaTeX commands.
#' 
#' \strong{Hierarchical Display:}
#' 
#' The \code{indent_groups} option is particularly useful for regression tables 
#' with categorical variables:
#' \itemize{
#'   \item Variable names appear unindented and bold
#'   \item Category levels appear indented beneath
#'   \item Reference categories clearly identified
#'   \item Creates professional hierarchical structure
#' }
#' 
#' \strong{Table Width Management:}
#' 
#' For tables wider than \code{\\textwidth}:
#' \enumerate{
#'   \item Use landscape orientation in LaTeX document
#'   \item Use \code{\\resizebox} to scale table
#'   \item Reduce font size in LaTeX: \code{\\small \\input{table.tex}}
#'   \item Use \code{condense_table = TRUE} to reduce columns
#'   \item Consider breaking across multiple tables
#' }
#' 
#' \strong{Workflow:}
#' 
#' Typical workflow for journal submission:
#' \enumerate{
#'   \item Generate table: \code{tbl2tex(results, "table1.tex")}
#'   \item Create LaTeX document with proper preamble
#'   \item Include table: \code{\\input{table1.tex}}
#'   \item Compile with pdflatex or other LaTeX engine
#'   \item Adjust formatting parameters as needed
#'   \item Regenerate and recompile
#' }
#'
#' @seealso
#' \code{\link{tbl2pdf}} for direct PDF output,
#' \code{\link{tbl2html}} for HTML tables,
#' \code{\link{tbl2docx}} for Word documents,
#' \code{\link{tbl2pptx}} for PowerPoint,
#' \code{\link{fit}} for regression tables,
#' \code{\link{desctbl}} for descriptive tables
#'
#' @examples
#' # Load data and create regression table
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' results <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     var_labels = clintrial_labels
#' )
#' 
#' # Example 1: Basic LaTeX export
#' tbl2tex(results, "basic.tex")
#' 
#' # Example 2: With booktabs for publication
#' tbl2tex(results, "publication.tex",
#'         booktabs = TRUE,
#'         caption = "Multivariable logistic regression results",
#'         label = "tab:regression")
#' 
#' # Example 3: Multi-line caption with abbreviations
#' tbl2tex(results, "detailed.tex",
#'         booktabs = TRUE,
#'         caption = "Table 1: Risk Factors for Mortality\\\\
#'                   aOR = adjusted odds ratio; CI = confidence interval\\\\
#'                   Model adjusted for age, sex, treatment, and disease stage",
#'         label = "tab:mortality")
#' 
#' # Example 4: Hierarchical display with indentation
#' tbl2tex(results, "indented.tex",
#'         indent_groups = TRUE,
#'         booktabs = TRUE)
#' 
#' # Example 5: Condensed table (reduced height)
#' tbl2tex(results, "condensed.tex",
#'         condense_table = TRUE,
#'         booktabs = TRUE)
#' 
#' # Example 6: With zebra stripes
#' tbl2tex(results, "striped.tex",
#'         zebra_stripes = TRUE,
#'         stripe_color = "gray!15",
#'         booktabs = TRUE)
#' # Remember to add \usepackage[table]{xcolor} to your LaTeX document
#' 
#' # Example 7: Dark header style
#' tbl2tex(results, "dark_header.tex",
#'         dark_header = TRUE,
#'         booktabs = TRUE)
#' # Requires \usepackage[table]{xcolor}
#' 
#' # Example 8: Custom cell padding
#' tbl2tex(results, "relaxed.tex",
#'         cell_padding = "relaxed",
#'         booktabs = TRUE)
#' 
#' # Example 9: Custom column alignment
#' tbl2tex(results, "aligned.tex",
#'         align = c("l", "l", "r", "r", "c"),
#'         booktabs = TRUE)
#' 
#' # Example 10: No header formatting (keep original names)
#' tbl2tex(results, "raw_headers.tex",
#'         format_headers = FALSE)
#' 
#' # Example 11: Disable significance bolding
#' tbl2tex(results, "no_bold.tex",
#'         bold_significant = FALSE,
#'         booktabs = TRUE)
#' 
#' # Example 12: Stricter significance threshold
#' tbl2tex(results, "strict_sig.tex",
#'         bold_significant = TRUE,
#'         sig_threshold = 0.01,  # Bold only if p < 0.01
#'         booktabs = TRUE)
#' 
#' # Example 13: Complete publication-ready table
#' tbl2tex(results, "final_table1.tex",
#'         booktabs = TRUE,
#'         caption = "Table 1: Multivariable Analysis of Mortality Risk Factors",
#'         label = "tab:main_results",
#'         indent_groups = TRUE,
#'         zebra_stripes = FALSE,  # Many journals prefer no stripes
#'         bold_significant = TRUE,
#'         cell_padding = "normal")
#' 
#' # Example 14: Descriptive statistics table
#' desc_table <- desctbl(
#'     data = clintrial,
#'     strata = "treatment",
#'     vars = c("age", "sex", "bmi"),
#'     var_labels = clintrial_labels
#' )
#' 
#' tbl2tex(desc_table, "table1_descriptive.tex",
#'         booktabs = TRUE,
#'         caption = "Table 1: Baseline Characteristics",
#'         label = "tab:baseline")
#' 
#' # Example 15: Model comparison table
#' models <- list(
#'     base = c("age", "sex"),
#'     full = c("age", "sex", "treatment", "stage")
#' )
#' 
#' comparison <- compfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     model_list = models
#' )
#' 
#' tbl2tex(comparison, "model_comparison.tex",
#'         booktabs = TRUE,
#'         caption = "Model Comparison Statistics",
#'         label = "tab:models")
#' 
#' # Example 16: For integration in LaTeX document
#' # After running tbl2tex(), use in LaTeX:
#' #
#' # \begin{table}[htbp]
#' #   \centering
#' #   \caption{Your Caption}
#' #   \label{tab:yourlabel}
#' #   \input{final_table1.tex}
#' # \end{table}
#' 
#' # Example 17: With resizing in LaTeX
#' # \begin{table}[htbp]
#' #   \centering
#' #   \caption{Wide Table}
#' #   \resizebox{\textwidth}{!}{\input{wide_results.tex}}
#' # \end{table}
#' 
#' # Example 18: Landscape table in LaTeX
#' # \usepackage{pdflscape}
#' # ...
#' # \begin{landscape}
#' #   \begin{table}[htbp]
#' #     \centering
#' #     \input{landscape_table.tex}
#' #   \end{table}
#' # \end{landscape}
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
