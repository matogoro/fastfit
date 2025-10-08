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
#' @param add_padding Logical. If TRUE, adds spacing around variable groups.
#'   Default is TRUE.
#' @param bold_significant Logical. If TRUE, bolds p-values below significance
#'   threshold. Default is TRUE.
#' @param sig_threshold Numeric. P-value threshold for boldface. Default is 0.05.
#' @param align Character string or vector specifying column alignment. If NULL,
#'   automatically determines based on content. Use "l" (left), "c" (center),
#'   or "r" (right).
#' @param indent_groups Logical. If TRUE, indents grouped rows hierarchically.
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
                     add_padding = TRUE,
                     bold_significant = TRUE, 
                     sig_threshold = 0.05,
                     align = NULL,
                     indent_groups = FALSE, 
                     ...) {
    
    if (!requireNamespace("xtable", quietly = TRUE)) {
        stop("Package 'xtable' required. Install with: install.packages('xtable')")
    }
    
    if (!grepl("\\.pdf$", tolower(file))) {
        stop("File must have .pdf extension")
    }
    
    if (!check_latex()) {
        if (requireNamespace("tinytex", quietly = TRUE) && tinytex::is_tinytex()) {
        }
        else {
            stop("PDF compilation requires LaTeX. Install TinyTeX with: tinytex::install_tinytex()")
        }
    }
    
    orientation <- match.arg(orientation, c("portrait", "landscape"))
    paper_settings <- get_paper_settings(paper, margins)
    df <- as.data.frame(table)

    has_n_row <- FALSE
    n_row_data <- NULL
    if (nrow(df) > 0 && "Variable" %in% names(df) && df$Variable[1] == "N") {
        has_n_row <- TRUE
        n_row_data <- df[1, ]
        df <- df[-1, ]
    }
    
    if (indent_groups) {
        df <- format_indented_groups(df)
    }
    
    if (bold_significant) {
        df <- format_pvalues_export_tex(df, sig_threshold)
    }
    
    if (add_padding && ("Variable" %in% names(df) || "variable" %in% 
                        names(df))) {
        df <- add_variable_padding(df)
    }
    
    if (is.null(align)) {
        align <- determine_alignment(df)
    }
    
    if (format_headers) {
        if (has_n_row) {
            names(df) <- format_column_headers_with_n_tex(names(df), n_row_data)
        } else {
            names(df) <- format_column_headers(names(df))
        }
    }
    
    xt <- xtable::xtable(df, align = align, ...)
    file_base <- tools::file_path_sans_ext(file)
    tex_file <- paste0(file_base, ".tex")
    use_standalone <- FALSE
    
    if (is.null(paper_settings$width)) {
        standalone_check <- system2("kpsewhich", args = "standalone.cls", 
                                    stdout = TRUE, stderr = FALSE)
        use_standalone <- length(standalone_check) > 0 && standalone_check != ""
        
        if (use_standalone) {
            message("Using standalone class for auto-sized output")
            cat(sprintf("\\documentclass[%dpt,border=10pt,varwidth=\\maxdimen]{standalone}\n
                         \\usepackage[T1]{fontenc}\n
                         \\usepackage[utf8]{inputenc}\n
                         \\usepackage{helvet}\n
                         \\renewcommand{\\familydefault}{\\sfdefault}\n
                         \\usepackage{array}\n
                         \\usepackage{graphicx}\n
                         \\usepackage{varwidth}\n
                         \\begin{document}\n
                         \\begin{varwidth}{\\linewidth}\n", 
                        font_size), file = tex_file)
            
        } else {
            
            has_pdfcrop <- Sys.which("pdfcrop") != ""
            
            if (has_pdfcrop) {
                message("Standalone class not found. Will use pdfcrop for auto-sizing")
            } else {
                warning("Auto-sizing requested but neither standalone class nor pdfcrop available.\n", 
                        "Install standalone with: tlmgr install standalone\n", 
                        "Using minimal margins instead.")
            }
            
            cat(sprintf("\\documentclass[%dpt]{article}\n
                         \\usepackage[T1]{fontenc}\n
                         \\usepackage[utf8]{inputenc}\n
                         \\usepackage[paperwidth=15in,paperheight=15in,margin=0.5in]{geometry}\n
                         \\usepackage{helvet}\n
                         \\renewcommand{\\familydefault}{\\sfdefault}\n
                         \\usepackage{array}\n
                         \\usepackage{graphicx}\n
                         \\pagestyle{empty}\n
                         \\begin{document}\n
                         \\noindent\n", 
                        font_size), file = tex_file)
        }
        
    } else {
        
        margin_str <- sprintf("top=%.1fin,right=%.1fin,bottom=%.1fin,left=%.1fin", 
                              paper_settings$margins[1], paper_settings$margins[2], 
                              paper_settings$margins[3], paper_settings$margins[4])
        cat(sprintf("\\documentclass[%dpt]{article}\n
                     \\usepackage[T1]{fontenc}\n
                     \\usepackage[utf8]{inputenc}\n
                     \\usepackage[%s,%s,%s]{geometry}\n
                     \\usepackage{helvet}\n
                     \\renewcommand{\\familydefault}{\\sfdefault}\n
                     \\usepackage{array}\n
                     \\usepackage{graphicx}\n
                     \\pagestyle{empty}\n
                     \\begin{document}\n", 
                    font_size, paper_settings$latex_paper, orientation, 
                    margin_str), file = tex_file)
    }
    
    if (fit_to_page && !is.null(paper_settings$width)) {
        cat("\n\\noindent\\resizebox{\\textwidth}{!}{%\n", file = tex_file, 
            append = TRUE)
    }
    
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
          sanitize.colnames.function = identity)
    
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
        warning("PDF compilation failed. Check ", file_base, 
                ".log for errors")
    }
    
    aux_files <- paste0(file_base, c(".aux", ".log"))
    
    for (f in aux_files) {
        if (file.exists(f)) 
            file.remove(f)
    }
    
    message(sprintf("Table exported to %s", file))
    
    invisible(NULL)
}
