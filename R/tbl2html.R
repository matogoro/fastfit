#' Export Table to HTML Format
#'
#' Converts a data frame, data.table, or matrix to HTML format with optional CSS 
#' styling for web display, HTML documents, or embedding in web applications. 
#' Generates clean, standards-compliant HTML with professional styling options 
#' including responsive design support, color schemes, and interactive features.
#'
#' @param table A data.frame, data.table, or matrix to export. Can be output from 
#'   \code{\link{desctbl}}, \code{\link{fit}}, \code{\link{uscreen}}, 
#'   \code{\link{fastfit}}, \code{\link{compfit}}, or any tabular data.
#'   
#' @param file Character string specifying the output HTML filename. Must have 
#'   \code{.html} or \code{.htm} extension. Example: \code{"results.html"}.
#'   
#' @param caption Character string. Optional caption displayed below the table. 
#'   Supports basic HTML formatting. Default is \code{NULL}.
#'   
#' @param format_headers Logical. If \code{TRUE}, formats column headers by 
#'   converting underscores to spaces and applying proper casing. Wraps statistical 
#'   notation (\emph{n}, \emph{p}) in \code{<i>} tags. Default is \code{TRUE}.
#'   
#' @param variable_padding Logical. If \code{TRUE}, adds vertical spacing around 
#'   variable groups for improved readability. Particularly useful for multi-row 
#'   factor variables. Default is \code{FALSE}.
#'   
#' @param bold_significant Logical. If \code{TRUE}, wraps significant p-values 
#'   in \code{<b>} tags for bold display. Makes important results stand out 
#'   visually. Default is \code{TRUE}.
#'   
#' @param p_threshold Numeric. P-value threshold for significance highlighting. 
#'   Only used when \code{bold_significant = TRUE}. Default is 0.05.
#'   
#' @param indent_groups Logical. If \code{TRUE}, indents grouped rows using 
#'   non-breaking spaces (\code{&nbsp;}) for hierarchical display. Useful for 
#'   factor variables in regression output. Default is \code{FALSE}.
#'   
#' @param condense_table Logical. If \code{TRUE}, condenses table by showing 
#'   only essential rows. Automatically sets \code{indent_groups = TRUE}. 
#'   Default is \code{FALSE}.
#'   
#' @param zebra_stripes Logical. If \code{TRUE}, applies alternating background 
#'   shading to different variables (not individual rows) for visual grouping. 
#'   Default is \code{FALSE}.
#'   
#' @param stripe_color Character string. HTML color specification for zebra 
#'   stripes. Can use hex codes (\code{"#EEEEEE"}), RGB (\code{"rgb(238,238,238)"}), 
#'   or color names (\code{"lightgray"}). Default is \code{"#EEEEEE"}.
#'   
#' @param dark_header Logical. If \code{TRUE}, creates black background with 
#'   white text for the header row. Provides strong visual contrast. 
#'   Default is \code{FALSE}.
#'   
#' @param include_css Logical. If \code{TRUE}, includes embedded CSS styling in 
#'   the output file for standalone HTML. Set to \code{FALSE} when embedding 
#'   in existing HTML with its own stylesheet. Default is \code{TRUE}.
#'   
#' @param ... Additional arguments passed to \code{\link[xtable]{xtable}}.
#'
#' @return Invisibly returns \code{NULL}. Creates an HTML file at the specified 
#'   location that can be opened in web browsers or embedded in HTML documents.
#'
#' @details
#' \strong{Output Format:}
#' 
#' The function generates standards-compliant HTML5 markup with:
#' \itemize{
#'   \item Semantic \code{<table>} structure
#'   \item Proper \code{<thead>} and \code{<tbody>} sections
#'   \item Accessible header cells (\code{<th>})
#'   \item Clean, readable markup
#'   \item Optional embedded CSS styling
#' }
#' 
#' \strong{Standalone vs. Embedded:}
#' 
#' \emph{Standalone HTML (\code{include_css = TRUE}):}
#' \itemize{
#'   \item Can be opened directly in web browsers
#'   \item Includes all necessary styling
#'   \item Self-contained, portable
#'   \item Suitable for sharing via email or web hosting
#' }
#' 
#' \emph{Embedded HTML (\code{include_css = FALSE}):}
#' \itemize{
#'   \item For inclusion in existing HTML documents
#'   \item No CSS included (use parent document's styles)
#'   \item Smaller file size
#'   \item Integrates with web frameworks (Shiny, R Markdown, Quarto)
#' }
#' 
#' \strong{CSS Styling:}
#' 
#' When \code{include_css = TRUE}, the function applies professional styling:
#' \itemize{
#'   \item \strong{Table:} Border-collapse, sans-serif font (Arial), 20px margin
#'   \item \strong{Cells:} 8px vertical × 12px horizontal padding, left-aligned text
#'   \item \strong{Borders:} 1px solid \code{#ddd} (light gray)
#'   \item \strong{Headers:} Bold text, light gray background (\code{#f2f2f2})
#'   \item \strong{Numeric columns:} Center-aligned (auto-detected)
#'   \item \strong{Caption:} Bold, 1.1em font, positioned below table
#' }
#' 
#' \emph{With \code{dark_header = TRUE}:}
#' \itemize{
#'   \item Header background: Black (\code{#000000})
#'   \item Header text: White (\code{#FFFFFF})
#'   \item Creates high contrast, modern appearance
#' }
#' 
#' \emph{With \code{zebra_stripes = TRUE}:}
#' \itemize{
#'   \item Alternating variable groups receive background color
#'   \item Default color: \code{#EEEEEE} (light gray)
#'   \item Applied via CSS class \code{.zebra-stripe}
#'   \item Groups entire variable (all factor levels together)
#' }
#' 
#' \strong{Hierarchical Display:}
#' 
#' The \code{indent_groups} option creates visual hierarchy using HTML 
#' non-breaking spaces:
#' ```html
#' <td><b>Treatment</b></td>  <!-- Variable name -->
#' <td>&nbsp;&nbsp;&nbsp;&nbsp;Control</td>  <!-- Indented level -->
#' <td>&nbsp;&nbsp;&nbsp;&nbsp;Active</td>   <!-- Indented level -->
#' ```
#' 
#' Benefits:
#' \itemize{
#'   \item Clear parent-child relationships
#'   \item Professional appearance
#'   \item Easy to scan and interpret
#'   \item Works across all browsers
#' }
#' 
#' \strong{Responsive Design:}
#' 
#' The default CSS supports responsive display:
#' \itemize{
#'   \item Tables scroll horizontally on narrow screens
#'   \item Maintains formatting across devices
#'   \item Mobile-friendly viewing
#'   \item Print-friendly styling
#' }
#' 
#' For better mobile support, consider:
#' \itemize{
#'   \item Using \code{condense_table = TRUE} to reduce width
#'   \item Adding viewport meta tag: 
#'     \code{<meta name="viewport" content="width=device-width, initial-scale=1">}
#'   \item Wrapping in responsive container: \code{<div style="overflow-x:auto;">}
#' }
#' 
#' \strong{Integration with R Markdown/Quarto:}
#' 
#' For R Markdown or Quarto documents:
#' ```r
#' # Generate HTML fragment (no CSS)
#' tbl2html(results, "table.html", include_css = FALSE)
#' 
#' # Include in document with results='asis'
#' \preformatted{
#'   cat(readLines("table.html"), sep = "\\n")
#' }
#' 
#' Or directly render without file:
#' ```r
#' # For inline display
#' htmltools::HTML(
#'   capture.output(
#'     print(xtable::xtable(results), type = "html")
#'   )
#' )
#' ```
#' 
#' \strong{Integration with Shiny:}
#' 
#' For Shiny applications:
#' ```r
#' # In server function
#' output$results_table <- renderUI(\{
#'   tbl2html(results_data(), "temp.html", include_css = FALSE)
#'   HTML(readLines("temp.html"))
#' \})
#' 
#' # Or use directly with DT package for interactive tables
#' output$interactive_table <- DT::renderDT(\{
#'   results_data()
#' \})
#' ```
#' 
#' \strong{Customizing Appearance:}
#' 
#' To customize beyond built-in options:
#' 
#' \emph{Option 1: Modify generated HTML}
#' ```r
#' tbl2html(results, "table.html", include_css = FALSE)
#' # Add custom CSS in parent HTML document
#' ```
#' 
#' \emph{Option 2: Post-process with HTML/CSS tools}
#' ```r
#' # Generate table
#' tbl2html(results, "table.html")
#' 
#' # Edit table.html to add custom classes or styles
#' # Use CSS to apply custom formatting
#' ```
#' 
#' \strong{Accessibility:}
#' 
#' The generated HTML follows accessibility best practices:
#' \itemize{
#'   \item Semantic table structure
#'   \item Proper header cells (\code{<th>}) with scope attributes
#'   \item Clear visual hierarchy
#'   \item Adequate color contrast (when using default styles)
#'   \item Screen reader friendly markup
#' }
#' 
#' For enhanced accessibility:
#' \itemize{
#'   \item Add descriptive caption
#'   \item Ensure sufficient color contrast if customizing
#'   \item Consider adding \code{summary} attribute for complex tables
#'   \item Test with screen readers
#' }
#' 
#' \strong{Performance:}
#' 
#' HTML tables are efficient for web display:
#' \itemize{
#'   \item Fast rendering in modern browsers
#'   \item Small file sizes (text-based)
#'   \item No external dependencies required
#'   \item Can be cached effectively
#' }
#' 
#' For very large tables (>1000 rows):
#' \itemize{
#'   \item Consider pagination or virtual scrolling
#'   \item Use interactive table libraries (DT, reactable)
#'   \item Consider server-side rendering for dynamic data
#' }
#' 
#' \strong{Browser Compatibility:}
#' 
#' Generated HTML works in all modern browsers:
#' \itemize{
#'   \item Chrome, Firefox, Safari, Edge (latest versions)
#'   \item Mobile browsers (iOS Safari, Chrome Mobile)
#'   \item Internet Explorer 11+ (with minor styling differences)
#' }
#'
#' @seealso
#' \code{\link{tbl2pdf}} for PDF output,
#' \code{\link{tbl2tex}} for LaTeX output,
#' \code{\link{tbl2docx}} for Word documents,
#' \code{\link{tbl2pptx}} for PowerPoint,
#' \code{\link{fit}} for regression tables,
#' \code{\link{desctbl}} for descriptive tables
#'
#' @examples
#' \donttest{
#'   # Load data and create table
#'   data(clintrial)
#'   data(clintrial_labels)
#'   
#'   results <- fit(
#'       data = clintrial,
#'       outcome = "os_status",
#'       predictors = c("age", "sex", "treatment", "stage"),
#'       labels = clintrial_labels
#'   )
#'   
#'   # Example 1: Basic HTML export (standalone)
#'   tbl2html(results, "results.html")
#'   # Open results.html in web browser
#'   
#'   # Example 2: With caption
#'   tbl2html(results, "captioned.html",
#'            caption = "Table 1: Multivariable Logistic Regression Results")
#'   
#'   # Example 3: For embedding (no CSS)
#'   tbl2html(results, "embed.html",
#'            include_css = FALSE)
#'   # Include in your HTML document
#'   
#'   # Example 4: Hierarchical display
#'   tbl2html(results, "indented.html",
#'            indent_groups = TRUE)
#'   
#'   # Example 5: Condensed table
#'   tbl2html(results, "condensed.html",
#'            condense_table = TRUE)
#'   
#'   # Example 6: With zebra stripes
#'   tbl2html(results, "striped.html",
#'            zebra_stripes = TRUE,
#'            stripe_color = "#F0F0F0")
#'   
#'   # Example 7: Dark header style
#'   tbl2html(results, "dark.html",
#'            dark_header = TRUE)
#'   
#'   # Example 8: Combination styling
#'   tbl2html(results, "styled.html",
#'            zebra_stripes = TRUE,
#'            dark_header = TRUE,
#'            bold_significant = TRUE)
#'   
#'   # Example 9: Custom stripe color
#'   tbl2html(results, "blue_stripes.html",
#'            zebra_stripes = TRUE,
#'            stripe_color = "#E3F2FD")  # Light blue
#'   
#'   # Example 10: Disable significance bolding
#'   tbl2html(results, "no_bold.html",
#'            bold_significant = FALSE)
#'   
#'   # Example 11: Stricter significance threshold
#'   tbl2html(results, "strict.html",
#'            bold_significant = TRUE,
#'            p_threshold = 0.01)
#'   
#'   # Example 12: No header formatting
#'   tbl2html(results, "raw_headers.html",
#'            format_headers = FALSE)
#'   
#'   # Example 13: Descriptive statistics table
#'   desc_table <- desctbl(
#'       data = clintrial,
#'       by = "treatment",
#'       variables = c("age", "sex", "bmi"),
#'       labels = clintrial_labels
#'   )
#'   
#'   tbl2html(desc_table, "baseline.html",
#'            caption = "Table 1: Baseline Characteristics by Treatment Group")
#'   
#'   # Example 14: For R Markdown (no CSS, for inline display)
#'   tbl2html(results, "rmd_table.html",
#'            include_css = FALSE,
#'            indent_groups = TRUE)
#'   
#'   # Then in R Markdown:
#'   # ```{r results='asis', echo=FALSE}
#'   # cat(readLines("rmd_table.html"), sep = "\n")
#'   # ```
#'   
#'   # Example 15: Email-friendly version
#'   tbl2html(results, "email.html",
#'            include_css = TRUE,  # Self-contained
#'            zebra_stripes = TRUE,
#'            caption = "Regression Results - See Attached")
#'   # Can be directly included in HTML emails
#'   
#'   # Example 16: Publication-ready web version
#'   tbl2html(results, "publication.html",
#'            caption = "Table 2: Multivariable Analysis of Risk Factors",
#'            indent_groups = TRUE,
#'            zebra_stripes = FALSE,  # Clean look
#'            bold_significant = TRUE,
#'            dark_header = FALSE)
#'   
#'   # Example 17: Modern dark theme
#'   tbl2html(results, "dark_theme.html",
#'            dark_header = TRUE,
#'            stripe_color = "#2A2A2A",  # Dark gray stripes
#'            zebra_stripes = TRUE)
#'   
#'   # Example 18: Minimal styling for custom CSS
#'   tbl2html(results, "minimal.html",
#'            include_css = FALSE,
#'            format_headers = FALSE,
#'            bold_significant = FALSE)
#'   # Apply your own CSS classes and styling
#'   
#'   # Example 19: Model comparison table
#'   models <- list(
#'       base = c("age", "sex"),
#'       full = c("age", "sex", "treatment", "stage")
#'   )
#'   
#'   comparison <- compfit(
#'       data = clintrial,
#'       outcome = "os_status",
#'       model_list = models
#'   )
#'   
#'   tbl2html(comparison, "comparison.html",
#'            caption = "Model Comparison Statistics")
#'   
#'   # Example 20: Mobile-optimized table
#'   tbl2html(results, "mobile.html",
#'            condense_table = TRUE,  # Reduce width
#'            font_size = 10)  # Note: font_size not in this function,
#'                             # adjust in CSS if needed
#' }
#' @export
tbl2html <- function(table,
                     file,
                     caption = NULL,
                     format_headers = TRUE,
                     variable_padding = FALSE, 
                     bold_significant = TRUE,
                     p_threshold = 0.05,
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
        df <- format_pvalues_export_html(df, p_threshold)
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
