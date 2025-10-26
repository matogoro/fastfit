#' Export Table to PDF Format
#'
#' Converts a data frame or data.table to a formatted PDF document using LaTeX
#' as an intermediate format. Supports automatic page sizing, landscape orientation,
#' and various formatting options for publication-ready output.
#'
#' @param table A data.frame, data.table, or matrix to export.
#' @param file Character string specifying the output PDF filename. Must have
#'   .pdf extension.
#' @param orientation Character string: "portrait" (default) or "landscape".
#' @param paper Character string specifying paper size: "letter" (default), "a4",
#'   "legal", or "auto" for content-based sizing.
#' @param margins Numeric vector of length 4 specifying margins in inches as
#'   c(top, right, bottom, left). Default is c(1, 1, 1, 1).
#' @param fit_to_page Logical. If TRUE, scales table to fit page width. Default
#'   is TRUE.
#' @param font_size Numeric. Base font size in points. Default is 8.
#' @param caption Character string. Optional caption to display below table.
#' @param format_headers Logical. If TRUE, formats column headers (e.g., converts
#'   underscores to spaces). Default is TRUE.
#' @param variable_padding Logical. If TRUE, adds spacing around variable groups.
#'   Default is TRUE.
#' @param cell_padding Character string or numeric. Add vertical padding to text
#'   within cells. Choices are "none", "normal" (default), "relaxed", "loose", or
#'   a number representing a custom multiplier.
#' @param bold_significant Logical. If TRUE, bolds p-values below significance
#'   threshold. Default is TRUE.
#' @param sig_threshold Numeric. P-value threshold for boldface. Default is 0.05.
#' @param align Character string or vector specifying column alignment. If NULL,
#'   automatically determines based on content. Use "l" (left), "c" (center),
#'   or "r" (right).
#' @param indent_groups Logical. If TRUE, indents grouped rows hierarchically.
#'   Default is FALSE.
#' @param condense_table Logical. If TRUE, decreases table vertical height by
#'   reducing continuous, survival, and binary categorical variables to single
#'   rows. Makes indent_groups = TRUE by default. Default is FALSE.
#' @param zebra_stripes Logical. If TRUE, adds alternating row shading. Default is FALSE.
#' @param stripe_color Character string. Uses LaTeX color names. Default is "gray!20".
#' @param show_logs Logical. If TRUE, keeps log files after PDF creation.
#'   Default is FALSE.
#' @param ... Additional arguments passed to xtable::xtable.
#'
#' @return Invisibly returns NULL. Creates a PDF file at the specified location.
#'
#' @details
#' This function uses LaTeX as an intermediate format to generate PDFs.  As such,
#' a working installation of a LaTeX distribution is required (e.g., TinyTeX,
#' TeX Live, MiKTeX, MacTeX, etc.).  Prior to execution, this function checks
#' for LaTeX availability and provides installation guidance if missing.
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
#'   \item \code{xcolor}
#'   \item \code{colortbl}
#' }
#' When \code{paper = "auto"}, the function attempts to size the PDF to content using
#' either the standalone document class or pdfcrop utility if available.
#' 
#' While this function attempts to gracefully handle special characters that need
#' escaping in LaTeX, user-inputted text needs to be properly escaped for the output
#' text to render appropriately (\emph{see} the multi-line caption example included
#' below).
#'
#' This function is tailored to outputs from \code{\link{desctbl()}},
#' \code{\link{fastfit()}}, \code{\link{uscreen()}}, \code{\link{fit()}},
#' and \code{\link{compfit()}}, although it can theoretically be applied to
#' any data frame or data.table object.
#' 
#' @examples
#' if (FALSE) {
#' # Basic export
#' tbl2pdf(results_table, "results.pdf")
#' 
#' # Landscape orientation with caption
#' tbl2pdf(wide_table, "wide_results.pdf",
#'         orientation = "landscape",
#'         caption = "Table 1: Regression Results")
#' 
#' # Multi-line caption
#' tbl2pdf(table, "results.pdf",
#'         caption = "Table 1 - Regression results\\\\
#'                   aHR = adjusted hazard ratio\\\\
#'                   \\textsuperscript{1}Note goes here")
#' 
#' # Auto-sized output with group indentation
#' tbl2pdf(grouped_table, "grouped.pdf",
#'         paper = "auto",
#'         indent_groups = TRUE)
#' 
#' # Custom margins and font size
#' tbl2pdf(table, "custom.pdf",
#'         margins = c(0.5, 0.5, 0.5, 0.5),
#'         font_size = 10)
#' }
#' @seealso
#' \code{\link{tbl2tex}} for LaTeX output,
#' \code{\link{tbl2html}} for HTML output.
#' 
#' @export
tbl2pdf <- function (table,
                     file,
                     orientation = "portrait",
                     paper = "letter", 
                     margins = NULL,
                     fit_to_page = TRUE,
                     font_size = 8,
                     caption = NULL, 
                     format_headers = TRUE,
                     variable_padding = TRUE,
                     cell_padding = "normal",
                     bold_significant = TRUE, 
                     sig_threshold = 0.05,
                     align = NULL,
                     indent_groups = FALSE,
                     condense_table = FALSE,
                     zebra_stripes = FALSE,
                     stripe_color = "gray!20",
                     dark_header = FALSE,
                     show_logs = TRUE,
                     ...) {
    
    if (!requireNamespace("xtable", quietly = TRUE)) {
        stop("Package 'xtable' required. Install with: install.packages('xtable')")
    }
    
    if (!grepl("\\.pdf$", tolower(file))) {
        stop("File must have .pdf extension")
    }
    
    if (!check_latex()) {
        if (requireNamespace("tinytex", quietly = TRUE) && tinytex::is_tinytex()) {
        } else {
            stop("PDF compilation requires LaTeX. Install TinyTeX with: tinytex::install_tinytex()")
        }
    }
    
    orientation <- match.arg(orientation, c("portrait", "landscape"))
    paper_settings <- get_paper_settings(paper, margins)
    df <- as.data.frame(table)

                                        # Detect variable groups BEFORE any processing
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
        df <- df[-1, ]
        if (!is.null(var_groups) && length(var_groups) > 1) {
            var_groups <- var_groups[-1]
        }
    }

    if (bold_significant) {
        df <- format_pvalues_export_tex(df, sig_threshold)
    }
    
    original_nrow <- nrow(df)

    if (condense_table) {
        indent_groups <- TRUE
        df <- condense_table_rows(df, indent_groups = indent_groups)
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
    
    padding_rows <- NULL
    if (variable_padding && ("Variable" %in% names(df) || "variable" %in% names(df))) {
        if (!is.null(var_groups)) {
            padding_positions <- which(diff(var_groups) != 0)
            new_var_groups <- integer(nrow(df) + length(padding_positions))
            current_pos <- 1
            for (i in 1:nrow(df)) {
                new_var_groups[current_pos] <- var_groups[i]
                current_pos <- current_pos + 1
                if (i %in% padding_positions) {
                    new_var_groups[current_pos] <- 0
                    current_pos <- current_pos + 1
                }
            }
            var_groups <- new_var_groups
        }
        df <- add_variable_padding(df)
    }
    
    add.to.row <- NULL
    if (zebra_stripes && "Variable" %in% names(df)) {
        is_indented <- indent_groups || condense_table
        
        if (is_indented) {
            var_starts <- which(!grepl("\\\\hspace", df$Variable) & 
                                df$Variable != "" & 
                                !is.na(df$Variable))
            
            commands <- character()
            positions <- numeric()
            
            for (i in seq_along(var_starts)) {
                if (i %% 2 != 0) {
                    start_idx <- var_starts[i]
                    end_idx <- if (i < length(var_starts)) {
                                   if (variable_padding) {
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
        } else {
            if ("Group" %in% names(df)) {
                var_starts <- which(df$Variable != "" & !is.na(df$Variable))
                
                commands <- character()
                positions <- numeric()
                
                for (i in seq_along(var_starts)) {
                    if (i %% 2 != 0) {
                        start_idx <- var_starts[i]
                        end_idx <- if (i < length(var_starts)) {
                                       if (variable_padding) {
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
    }

    if (is.null(align)) {
        align <- determine_alignment(df)
    }
    
    original_col_names <- names(df)
    
    if (format_headers) {
        if (has_n_row) {
            names(df) <- format_column_headers_with_n_tex(names(df), n_row_data)
        } else {
            names(df) <- format_column_headers(names(df))
        }
    }

    if (dark_header) {
        col_names <- names(df)
        for (i in seq_along(col_names)) {
            col_names[i] <- paste0("\\color{white}", col_names[i])
        }
        names(df) <- col_names
        
        header_command <- "\\rowcolor{black}"
        
        if (!is.null(add.to.row) && length(add.to.row$pos) > 0) {
            add.to.row$pos <- c(list(-1), add.to.row$pos)
            add.to.row$command <- c(header_command, add.to.row$command)
        } else {
            add.to.row <- list(
                pos = list(-1),
                command = header_command
            )
        }
    }
    
                                        # Create xtable object
    xt <- xtable::xtable(df, align = align, ...)
    
    file_base <- tools::file_path_sans_ext(file)
    tex_file <- paste0(file_base, ".tex")
    use_standalone <- FALSE
    
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
    
    if (is.null(paper_settings$width)) {
        standalone_check <- system2("kpsewhich", args = "standalone.cls", 
                                    stdout = TRUE, stderr = FALSE)
        use_standalone <- length(standalone_check) > 0 && standalone_check != ""
        
        if (use_standalone) {
            message("Using standalone class for auto-sized output")
            xcolor_line <- if (zebra_stripes || dark_header) "\\usepackage[table]{xcolor}\n" else ""
                                        # Add arraystretch command if padding requested
            array_stretch_line <- if (!is.null(array_stretch)) 
                                      sprintf("\\renewcommand{\\arraystretch}{%s}\n", array_stretch) else ""
            
            cat(sprintf("\\documentclass[%dpt,border=10pt,varwidth=\\maxdimen]{standalone}\n
                         \\usepackage[T1]{fontenc}\n
                         \\usepackage[utf8]{inputenc}\n
                         %s\\usepackage{helvet}\n
                         \\renewcommand{\\familydefault}{\\sfdefault}\n
                         \\usepackage{array}\n
                         \\usepackage{graphicx}\n
                         \\usepackage{varwidth}\n
                         %s\\begin{document}\n
                         \\begin{varwidth}{\\linewidth}\n", 
                        font_size, xcolor_line, array_stretch_line), file = tex_file)
            
        } else {
            
            has_pdfcrop <- Sys.which("pdfcrop") != ""
            
            if (has_pdfcrop) {
                message("Standalone class not found. Will use pdfcrop for auto-sizing")
            } else {
                warning("Auto-sizing requested but neither standalone class nor pdfcrop available.\n", 
                        "Install standalone with: tlmgr install standalone\n", 
                        "Using minimal margins instead.")
            }
            
            xcolor_line <- if (zebra_stripes || dark_header) "\\usepackage[table]{xcolor}\n" else ""
            array_stretch_line <- if (!is.null(array_stretch)) 
                                      sprintf("\\renewcommand{\\arraystretch}{%s}\n", array_stretch) else ""
            
            cat(sprintf("\\documentclass[%dpt]{article}\n
                         \\usepackage[T1]{fontenc}\n
                         \\usepackage[utf8]{inputenc}\n
                         \\usepackage[paperwidth=15in,paperheight=15in,margin=0.5in]{geometry}\n
                         %s\\usepackage{helvet}\n
                         \\renewcommand{\\familydefault}{\\sfdefault}\n
                         \\usepackage{array}\n
                         \\usepackage{graphicx}\n
                         \\pagestyle{empty}\n
                         %s\\begin{document}\n
                         \\noindent\n", 
                        font_size, xcolor_line, array_stretch_line), file = tex_file)
        }
        
    } else {
        
        margin_str <- sprintf("top=%.1fin,right=%.1fin,bottom=%.1fin,left=%.1fin", 
                              paper_settings$margins[1], paper_settings$margins[2], 
                              paper_settings$margins[3], paper_settings$margins[4])
        xcolor_line <- if (zebra_stripes || dark_header) "\\usepackage[table]{xcolor}\n" else ""
        array_stretch_line <- if (!is.null(array_stretch)) 
                                  sprintf("\\renewcommand{\\arraystretch}{%s}\n", array_stretch) else ""
        
        cat(sprintf("\\documentclass[%dpt]{article}\n
                     \\usepackage[T1]{fontenc}\n
                     \\usepackage[utf8]{inputenc}\n
                     \\usepackage[%s,%s,%s]{geometry}\n
                     %s\\usepackage{helvet}\n
                     \\renewcommand{\\familydefault}{\\sfdefault}\n
                     \\usepackage{array}\n
                     \\usepackage{graphicx}\n
                     \\pagestyle{empty}\n
                     %s\\begin{document}\n", 
                    font_size, paper_settings$latex_paper, orientation, 
                    margin_str, xcolor_line, array_stretch_line), file = tex_file)
    }
    
    if (fit_to_page && !is.null(paper_settings$width)) {
        cat("\n\\noindent\\resizebox{\\textwidth}{!}{%\n", file = tex_file, 
            append = TRUE)
    }

    if ((zebra_stripes || dark_header) && length(add.to.row$pos) > 0) {
        print(xt,
              file = tex_file,
              append = TRUE,
              include.rownames = FALSE, 
              booktabs = FALSE,
              floating = FALSE,
              tabular.environment = "tabular", 
              hline.after = c(-1, 0, nrow(xt)),
              add.to.row = add.to.row,
              sanitize.text.function = sanitize_for_latex, 
              sanitize.rownames.function = sanitize_for_latex,
              sanitize.colnames.function = function(x) x)
    } else {
        print(xt,
              file = tex_file,
              append = TRUE,
              include.rownames = FALSE, 
              booktabs = FALSE,
              floating = FALSE,
              tabular.environment = "tabular", 
              hline.after = c(-1, 0, nrow(xt)),
              sanitize.text.function = sanitize_for_latex, 
              sanitize.rownames.function = sanitize_for_latex,
              sanitize.colnames.function = function(x) x)
    }
    
    if (fit_to_page && !is.null(paper_settings$width)) {
        cat("}\n", file = tex_file, append = TRUE)
    }
    
    if (!is.null(caption)) {
        cat(sprintf("\n\n\\vspace{1em}\\noindent %s", caption), 
            file = tex_file, append = TRUE)
    }
    
    if (use_standalone) {
        cat("\n\\end{varwidth}\n", file = tex_file, append = TRUE)
    }
    
    cat("\n\\end{document}", file = tex_file, append = TRUE)

    message("Compiling PDF...")
    
    result <- system2("pdflatex",
                      args = c("-interaction=nonstopmode", tex_file),
                      stdout = FALSE,
                      stderr = FALSE)
    
    pdf_file <- paste0(file_base, ".pdf")
    
    if (is.null(paper_settings$width) && !use_standalone && file.exists(pdf_file) && Sys.which("pdfcrop") != "") {
        message("Cropping PDF to content...")
        temp_pdf <- paste0(file_base, "_temp.pdf")
        file.rename(pdf_file, temp_pdf)
        system2("pdfcrop", args = c("--margins", "10", temp_pdf, pdf_file), stdout = FALSE, stderr = FALSE)
        file.remove(temp_pdf)
    }
    
    if (!file.exists(pdf_file)) {
        warning("PDF compilation failed. Ensure show_logs = TRUE and check ", file_base, 
                ".log for errors")
    }
    
    aux_files <- paste0(file_base, c(".aux", ".log"))

    if (!show_logs) {
        for (f in aux_files) {
            if (file.exists(f))
                file.remove(f)
        }
    }
    
    message(sprintf("Table exported to %s", file))
    
    invisible(NULL)
}
