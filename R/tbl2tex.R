#' Export Table to LaTeX Format
#'
#' Converts a data frame or data.table to LaTeX format suitable for inclusion
#' in LaTeX documents. Supports various table packages and formatting options.
#'
#' Tailored to outputs from desctbl(), fastfit(), uscreen(), fit(), and compfit().
#'
#' @param table A data.frame, data.table, or matrix to export.
#' @param file Character string specifying the output .tex filename. Must have
#'   .tex extension.
#' @param caption Character string. Table caption for \\caption{} command.
#' @param format_headers Logical. If TRUE, formats column headers by converting
#'   underscores to spaces and applying title case. Default is TRUE.
#' @param add_padding Logical. If TRUE, adds vertical spacing around variable
#'   groups for better readability. Default is TRUE.
#' @param bold_significant Logical. If TRUE, wraps significant p-values in
#'   \\textbf{} commands. Default is TRUE.
#' @param sig_threshold Numeric. P-value threshold for bolding (e.g., 0.05).
#'   Values below this are considered significant. Default is 0.05.
#' @param align Character string or vector specifying column alignment. If NULL,
#'   automatically determines based on column content. Options: "l" (left),
#'   "c" (center), "r" (right), or "p{width}" for paragraph columns.
#' @param indent_groups Logical. If TRUE, uses \\hspace to indent grouped rows
#'   creating a hierarchical display. Default is FALSE.
#' @param booktabs Logical. If TRUE, uses booktabs package commands for
#'   professional-quality rules. Default is FALSE.
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
#' When booktabs = TRUE, uses \\toprule, \\midrule, and \\bottomrule for
#' cleaner table appearance. This requires \\usepackage{booktabs} in the
#' LaTeX preamble.
#' 
#' The indent_groups option is particularly useful for regression tables with
#' categorical variables, creating a clear visual hierarchy.
#'
#' @examples
#' \dontrun{
#' # Basic LaTeX export
#' tbl2tex(results, "results.tex")
#' 
#' # Publication-quality table with booktabs
#' tbl2tex(final_results, "publication.tex",
#'         booktabs = TRUE,
#'         caption = "Multivariable regression analysis",
#'         label = "tab:regression")
#' 
#' # Hierarchical table with indentation
#' tbl2tex(model_output, "model.tex",
#'         indent_groups = TRUE,
#'         bold_significant = TRUE)
#' 
#' # Custom alignment
#' tbl2tex(wide_table, "aligned.tex",
#'         align = c("l", "l", "r", "r", "c"))
#' }
#'
#' @seealso
#' \code{\link{tbl2pdf}} for PDF output,
#' \code{\link{tbl2html}} for HTML output,
#'
#' @export
tbl2tex <- function (table,
                     file,
                     format_headers = TRUE,
                     add_padding = TRUE, 
                     bold_significant = TRUE,
                     sig_threshold = 0.05,
                     align = NULL, 
                     indent_groups = FALSE,
                     booktabs = FALSE,
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
    
    if (indent_groups) {
        df <- format_indented_groups(df)
    }
    
    if (bold_significant) {
        df <- format_pvalues_export_tex(df, sig_threshold)
    }
    
    if (add_padding && ("Variable" %in% names(df) || "variable" %in% names(df))) {
        df <- add_variable_padding(df)
    }
    
    if (is.null(align)) {
        align <- determine_alignment(df)
    }
    
    if (format_headers) {
        names(df) <- format_column_headers(names(df))
    }
    
    xt <- xtable::xtable(df,
                         align = align,
                         caption = caption, 
                         label = label, ...)
    
    print(xt,
          file = file,
          include.rownames = FALSE,
          booktabs = booktabs, 
          floating = FALSE,
          sanitize.text.function = sanitize_for_latex, 
          sanitize.rownames.function = sanitize_for_latex,
          sanitize.colnames.function = identity, 
          hline.after = c(-1, 0, nrow(xt)))
    
    message(sprintf("Table exported to %s", file))
    
    invisible(NULL)
}
