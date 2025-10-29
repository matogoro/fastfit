#' Export Table to PDF Format
#'
#' Converts a data frame, data.table, or matrix to a professionally formatted PDF 
#' document using LaTeX as an intermediate format. Provides extensive control over 
#' page layout, typography, and formatting for publication-ready output. Particularly 
#' well-suited for tables from regression analyses, descriptive statistics, and 
#' model comparisons.
#'
#' @param table A data.frame, data.table, or matrix to export. Can be output from 
#'   \code{\link{desctbl}}, \code{\link{fit}}, \code{\link{uscreen}}, 
#'   \code{\link{fastfit}}, \code{\link{compfit}}, or any tabular data structure.
#'   
#' @param file Character string specifying the output PDF filename. Must have 
#'   \code{.pdf} extension. Example: \code{"results.pdf"}, \code{"Table1.pdf"}.
#'   
#' @param orientation Character string specifying page orientation: 
#'   \itemize{
#'     \item \code{"portrait"} - Vertical orientation [default]
#'     \item \code{"landscape"} - Horizontal orientation (recommended for wide tables)
#'   }
#'   
#' @param paper Character string specifying paper size:
#'   \itemize{
#'     \item \code{"letter"} - US Letter (8.5" × 11") [default]
#'     \item \code{"a4"} - A4 (210mm × 297mm)
#'     \item \code{"legal"} - US Legal (8.5" × 14")
#'     \item \code{"auto"} - Auto-size to content (no margins, crops to fit)
#'   }
#'   
#' @param margins Numeric vector of length 4 specifying margins in inches as 
#'   \code{c(top, right, bottom, left)}. Default is \code{c(1, 1, 1, 1)}. 
#'   Ignored when \code{paper = "auto"}.
#'   
#' @param fit_to_page Logical. If \code{TRUE}, scales table to fit within the 
#'   text width (respects margins). Useful for wide tables that would otherwise 
#'   overflow. Default is \code{TRUE}.
#'   
#' @param font_size Numeric. Base font size in points. Default is 8. Smaller 
#'   values accommodate more content; larger values improve readability. 
#'   Typical range: 6-12 points.
#'   
#' @param caption Character string. Optional caption displayed below the table. 
#'   Supports LaTeX formatting for multi-line captions, superscripts, italics, etc. 
#'   See Details for formatting guidance. Default is \code{NULL}.
#'   
#' @param format_headers Logical. If \code{TRUE}, applies automatic formatting 
#'   to column headers: converts underscores to spaces, italicizes statistical 
#'   notation (\emph{n}, \emph{p}), and improves readability. Default is \code{TRUE}.
#'   
#' @param var_padding Logical. If \code{TRUE}, adds vertical spacing between 
#'   different variables in the table, creating visual grouping. Particularly 
#'   useful for regression tables with multiple predictors. Default is \code{TRUE}.
#'   
#' @param cell_padding Character string or numeric specifying vertical padding 
#'   within table cells:
#'   \itemize{
#'     \item \code{"none"} - No extra padding (most compact)
#'     \item \code{"normal"} - Standard padding [default]
#'     \item \code{"relaxed"} - Increased padding
#'     \item \code{"loose"} - Maximum padding
#'     \item Numeric value - Custom multiplier (e.g., \code{1.5})
#'   }
#'   Adjusts \code{\\arraystretch} in LaTeX.
#'   
#' @param bold_significant Logical. If \code{TRUE}, applies bold formatting to 
#'   p-values below the significance threshold, making significant results stand 
#'   out visually. Default is \code{TRUE}.
#'   
#' @param sig_threshold Numeric. P-value threshold for bold formatting. Only 
#'   used when \code{bold_significant = TRUE}. Default is 0.05. Common alternatives: 
#'   0.01, 0.10.
#'   
#' @param align Character string or vector specifying column alignment. Options:
#'   \itemize{
#'     \item \code{"l"} - Left aligned
#'     \item \code{"c"} - Center aligned
#'     \item \code{"r"} - Right aligned
#'   }
#'   If \code{NULL} (default), automatically determines alignment based on content 
#'   (text left, numbers right). Can specify per-column: \code{c("l", "l", "r", "r")}.
#'   
#' @param indent_groups Logical. If \code{TRUE}, indents factor levels/groups 
#'   under their parent variable using horizontal space, creating a hierarchical 
#'   display. Useful for factor variables in regression tables. Default is \code{FALSE}.
#'   
#' @param condense_table Logical. If \code{TRUE}, condenses the table by:
#'   \itemize{
#'     \item Showing only one row per continuous variable (estimate + CI)
#'     \item Showing only non-reference categories for binary factors
#'     \item Automatically setting \code{indent_groups = TRUE}
#'   }
#'   Significantly reduces table height. Default is \code{FALSE}.
#'   
#' @param zebra_stripes Logical. If \code{TRUE}, applies alternating gray 
#'   background shading to different variables (not individual rows) for improved 
#'   visual grouping and readability. Default is \code{FALSE}.
#'   
#' @param stripe_color Character string. LaTeX color specification for zebra 
#'   stripes. Default is \code{"gray!20"} (20\% gray). Can use other colors like 
#'   \code{"blue!10"}, \code{"red!15"}. Requires \code{zebra_stripes = TRUE}.
#'   
#' @param dark_header Logical. If \code{TRUE}, creates a dark (black) background 
#'   with white text for the header row. Provides strong visual contrast. 
#'   Default is \code{FALSE}.
#'   
#' @param show_logs Logical. If \code{TRUE}, retains LaTeX log and auxiliary 
#'   files after PDF compilation for troubleshooting. If \code{FALSE}, deletes 
#'   these files. Default is \code{TRUE}.
#'   
#' @param ... Additional arguments passed to \code{\link[xtable]{xtable}} for 
#'   advanced LaTeX table customization.
#'
#' @return Invisibly returns \code{NULL}. Creates a PDF file at the specified 
#'   location. If compilation fails, check the \code{.log} file (if 
#'   \code{show_logs = TRUE}) for error details.
#'
#' @details
#' \strong{LaTeX Requirements:}
#' 
#' This function requires a working LaTeX installation. The function checks for 
#' LaTeX availability and provides installation guidance if missing.
#' 
#' \emph{Recommended LaTeX distributions:}
#' \itemize{
#'   \item \strong{TinyTeX} (lightweight, R-integrated): Install via 
#'     \code{tinytex::install_tinytex()}
#'   \item \strong{TeX Live} (comprehensive, cross-platform)
#'   \item \strong{MiKTeX} (Windows)
#'   \item \strong{MacTeX} (macOS)
#' }
#' 
#' \emph{Required LaTeX packages} (auto-installed with most distributions):
#' \itemize{
#'   \item fontenc, inputenc - Character encoding
#'   \item array, booktabs, longtable - Table formatting
#'   \item graphicx - Scaling tables
#'   \item geometry - Page layout
#'   \item pdflscape, lscape - Landscape orientation
#'   \item helvet - Sans-serif fonts
#'   \item standalone, varwidth - Auto-sizing (for \code{paper = "auto"})
#'   \item float, caption - Floats and captions
#'   \item xcolor, colortbl - Colors (for \code{zebra_stripes} or \code{dark_header})
#' }
#' 
#' \strong{Caption Formatting:}
#' 
#' Captions support LaTeX commands for rich formatting:
#' ```r
#' # Multi-line caption with line breaks
#' caption = "Table 1: Multivariable Analysis\\\\
#'           OR = odds ratio; CI = confidence interval"
#' 
#' # With superscripts and italics
#' caption = "Table 1: Results\\\\
#'           \\textsuperscript{a}Adjusted for age\\\\
#'           \\textit{p}-values from Wald tests"
#' 
#' # With special characters (must escape)
#' caption = "Results for \\$1000 increase in income (\\%)"
#' ```
#' 
#' LaTeX special characters that need escaping: \code{& % $ # _ { } ~ ^ \\}
#' 
#' \strong{Auto-Sizing (paper = "auto"):}
#' 
#' When \code{paper = "auto"}, the function attempts to create a minimal PDF 
#' sized exactly to the table content:
#' \enumerate{
#'   \item Tries to use the \code{standalone} LaTeX class (cleanest output)
#'   \item Falls back to \code{pdfcrop} utility if standalone unavailable
#'   \item Falls back to minimal margins if neither available
#' }
#' 
#' Auto-sized PDFs are ideal for:
#' \itemize{
#'   \item Including in other documents
#'   \item Flexible embedding in presentations
#'   \item Maximum space efficiency
#' }
#' 
#' \strong{Table Width Management:}
#' 
#' For wide tables that don't fit on the page:
#' \enumerate{
#'   \item Use \code{orientation = "landscape"}
#'   \item Use \code{fit_to_page = TRUE} (default) to auto-scale
#'   \item Reduce \code{font_size} (e.g., 7 or 6)
#'   \item Use \code{condense_table = TRUE} to reduce columns/rows
#'   \item Consider \code{paper = "auto"} for maximum flexibility
#' }
#' 
#' \strong{Optimized for FastFit Tables:}
#' 
#' This function is specifically designed for tables produced by:
#' \itemize{
#'   \item \code{\link{desctbl}} - Descriptive statistics tables
#'   \item \code{\link{fit}} - Single model regression tables
#'   \item \code{\link{uscreen}} - Univariable screening tables
#'   \item \code{\link{fastfit}} - Combined univariable/multivariable tables
#'   \item \code{\link{compfit}} - Model comparison tables
#' }
#' 
#' The function automatically detects and formats these table types, including:
#' \itemize{
#'   \item Sample size rows (N = X)
#'   \item Variable grouping
#'   \item P-value highlighting
#'   \item Proper alignment of estimates and CIs
#' }
#' 
#' \strong{Troubleshooting:}
#' 
#' If PDF compilation fails:
#' \enumerate{
#'   \item Check that LaTeX is installed: Run \code{Sys.which("pdflatex")}
#'   \item Set \code{show_logs = TRUE} and examine the .log file
#'   \item Common issues:
#'     \itemize{
#'       \item Missing LaTeX packages: Install via package manager
#'       \item Special characters in text: Escape properly
#'       \item Very wide tables: Use landscape or reduce font size
#'       \item Caption formatting: Check LaTeX syntax
#'     }
#' }
#'
#' @seealso
#' \code{\link{tbl2tex}} for LaTeX source files,
#' \code{\link{tbl2html}} for HTML output,
#' \code{\link{tbl2docx}} for Microsoft Word,
#' \code{\link{tbl2pptx}} for PowerPoint,
#' \code{\link{desctbl}} for descriptive tables,
#' \code{\link{fit}} for regression tables
#'
#' @examples
#' # Load example data and create a regression table
#' data(clintrial)
#' data(clintrial_labels)
#' library(survival)
#' 
#' # Create a regression table
#' results <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     var_labels = clintrial_labels
#' )
#' 
#' # Example 1: Basic PDF export
#' tbl2pdf(results, "basic_results.pdf")
#' 
#' # Example 2: Landscape orientation for wide tables
#' tbl2pdf(results, "wide_results.pdf",
#'         orientation = "landscape")
#' 
#' # Example 3: With caption
#' tbl2pdf(results, "captioned.pdf",
#'         caption = "Table 1: Multivariable logistic regression results")
#' 
#' # Example 4: Multi-line caption with formatting
#' tbl2pdf(results, "formatted_caption.pdf",
#'         caption = "Table 1: Risk Factors for Mortality\\\\
#'                   aOR = adjusted odds ratio; CI = confidence interval\\\\
#'                   \\textsuperscript{a}Adjusted for age, sex, and treatment")
#' 
#' # Example 5: Auto-sized PDF (no fixed page dimensions)
#' tbl2pdf(results, "autosize.pdf",
#'         paper = "auto")
#' 
#' # Example 6: A4 paper with custom margins
#' tbl2pdf(results, "a4_custom.pdf",
#'         paper = "a4",
#'         margins = c(0.75, 0.75, 0.75, 0.75))
#' 
#' # Example 7: Larger font for readability
#' tbl2pdf(results, "large_font.pdf",
#'         font_size = 11)
#' 
#' # Example 8: Indented hierarchical display
#' tbl2pdf(results, "indented.pdf",
#'         indent_groups = TRUE)
#' 
#' # Example 9: Condensed table (reduced height)
#' tbl2pdf(results, "condensed.pdf",
#'         condense_table = TRUE)
#' 
#' # Example 10: With zebra stripes
#' tbl2pdf(results, "striped.pdf",
#'         zebra_stripes = TRUE,
#'         stripe_color = "gray!15")
#' 
#' # Example 11: Dark header style
#' tbl2pdf(results, "dark_header.pdf",
#'         dark_header = TRUE)
#' 
#' # Example 12: Combination of formatting options
#' tbl2pdf(results, "publication_ready.pdf",
#'         orientation = "portrait",
#'         paper = "letter",
#'         font_size = 9,
#'         caption = "Table 2: Multivariable Analysis\\\\
#'                   Model adjusted for age, sex, and clinical factors",
#'         indent_groups = TRUE,
#'         zebra_stripes = TRUE,
#'         bold_significant = TRUE,
#'         sig_threshold = 0.05)
#' 
#' # Example 13: Adjust cell padding
#' tbl2pdf(results, "relaxed_padding.pdf",
#'         cell_padding = "relaxed")  # More spacious
#' 
#' # Example 14: No scaling (natural table width)
#' tbl2pdf(results, "no_scale.pdf",
#'         fit_to_page = FALSE,
#'         font_size = 10)
#' 
#' # Example 15: Hide significance bolding
#' tbl2pdf(results, "no_bold.pdf",
#'         bold_significant = FALSE)
#' 
#' # Example 16: Custom column alignment
#' tbl2pdf(results, "custom_align.pdf",
#'         align = c("l", "l", "r", "r", "c"))
#' 
#' # Example 17: Descriptive statistics table
#' desc_table <- desctbl(
#'     data = clintrial,
#'     strata = "treatment",
#'     vars = c("age", "sex", "bmi", "stage"),
#'     var_labels = clintrial_labels
#' )
#' 
#' tbl2pdf(desc_table, "descriptive.pdf",
#'         caption = "Table 1: Baseline Characteristics by Treatment Group",
#'         landscape = TRUE)
#' 
#' # Example 18: Model comparison table
#' models <- list(
#'     base = c("age", "sex"),
#'     full = c("age", "sex", "bmi", "treatment")
#' )
#' 
#' comparison <- compfit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     model_list = models
#' )
#' 
#' tbl2pdf(comparison, "model_comparison.pdf",
#'         caption = "Table 3: Model Comparison Statistics")
#' 
#' # Example 19: Very wide table with aggressive fitting
#' wide_model <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "race", "bmi", "smoking", 
#'                   "hypertension", "diabetes", "treatment", "stage")
#' )
#' 
#' tbl2pdf(wide_model, "very_wide.pdf",
#'         orientation = "landscape",
#'         paper = "legal",  # Even wider than letter
#'         font_size = 7,
#'         fit_to_page = TRUE,
#'         condense_table = TRUE)
#' 
#' # Example 20: Troubleshooting - keep logs
#' tbl2pdf(results, "debug.pdf",
#'         show_logs = TRUE)
#' # If it fails, check debug.log for error messages
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
                     var_padding = TRUE,
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
    if (var_padding && ("Variable" %in% names(df) || "variable" %in% names(df))) {
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
        } else {
            if ("Group" %in% names(df)) {
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
