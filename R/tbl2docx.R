#' Export Table to Microsoft Word Format (DOCX)
#'
#' Converts a data frame, data.table, or matrix to a fully editable Microsoft Word 
#' document (.docx) using the flextable and officer packages. Creates professional, 
#' publication-ready tables with extensive formatting options including typography, 
#' alignment, colors, and page layout. Tables can be further edited in Microsoft 
#' Word after creation.
#'
#' @param table A data.frame, data.table, or matrix to export. Can be output from 
#'   \code{\link{desctbl}}, \code{\link{fit}}, \code{\link{uscreen}}, 
#'   \code{\link{fastfit}}, \code{\link{compfit}}, or any tabular data.
#'   
#' @param file Character string specifying the output DOCX filename. Must have 
#'   \code{.docx} extension. Example: \code{"results.docx"}, \code{"Table1.docx"}.
#'   
#' @param caption Character string. Optional caption displayed above the table 
#'   in the Word document. Default is \code{NULL}.
#'   
#' @param font_size Numeric. Base font size in points for table content. 
#'   Default is 8. Typical range: 8-12 points. Headers use slightly larger size.
#'   
#' @param font_family Character string. Font family name for the table. Must be 
#'   a font installed on the system. Default is \code{"Arial"}. Common options: 
#'   \code{"Times New Roman"}, \code{"Calibri"}, \code{"Helvetica"}.
#'   
#' @param format_headers Logical. If \code{TRUE}, formats column headers by 
#'   italicizing statistical notation (\emph{n}, \emph{p}), converting underscores 
#'   to spaces, and improving readability. Default is \code{TRUE}.
#'   
#' @param bold_significant Logical. If \code{TRUE}, applies bold formatting to 
#'   p-values below the significance threshold. Makes significant results stand 
#'   out. Default is \code{TRUE}.
#'   
#' @param sig_threshold Numeric. P-value threshold for bold formatting. Only 
#'   used when \code{bold_significant = TRUE}. Default is 0.05.
#'   
#' @param indent_groups Logical. If \code{TRUE}, indents factor levels under 
#'   their parent variable using horizontal spacing, creating hierarchical display. 
#'   Useful for categorical variables in regression tables. Default is \code{FALSE}.
#'   
#' @param condense_table Logical. If \code{TRUE}, condenses table by showing 
#'   only essential rows (single row for continuous, non-reference for binary). 
#'   Automatically sets \code{indent_groups = TRUE}. Significantly reduces table 
#'   height. Default is \code{FALSE}.
#'   
#' @param zebra_stripes Logical. If \code{TRUE}, applies alternating row shading 
#'   to different variables (not individual rows) for visual grouping. 
#'   Default is \code{FALSE}.
#'   
#' @param dark_header Logical. If \code{TRUE}, creates a dark background with 
#'   light text for the header row, providing strong visual contrast. 
#'   Default is \code{FALSE}.
#'   
#' @param paper Character string specifying paper size:
#'   \itemize{
#'     \item \code{"letter"} - US Letter (8.5" × 11") [default]
#'     \item \code{"a4"} - A4 (210mm × 297mm)
#'     \item \code{"legal"} - US Legal (8.5" × 14")
#'   }
#'   
#' @param orientation Character string specifying page orientation:
#'   \itemize{
#'     \item \code{"portrait"} - Vertical [default]
#'     \item \code{"landscape"} - Horizontal (for wide tables)
#'   }
#'   
#' @param width Numeric. Table width in inches. If \code{NULL} (default), 
#'   automatically fits to content and page width. Specify to control exactly.
#'   
#' @param align Character vector specifying column alignment for each column. 
#'   Options: \code{"left"}, \code{"center"}, or \code{"right"}. If \code{NULL} 
#'   (default), automatically determines based on content (text left, numbers right).
#'   Example: \code{c("left", "left", "center", "right", "right")}.
#'   
#' @param return_ft Logical. If \code{TRUE}, returns the flextable object 
#'   directly for further customization. If \code{FALSE} (default), returns 
#'   invisibly with flextable as attribute. See Details for usage. 
#'   Default is \code{FALSE}.
#'   
#' @param ... Additional arguments passed to \code{\link[officer]{read_docx}} 
#'   for document initialization.
#'
#' @return Behavior depends on \code{return_ft}:
#'   \describe{
#'     \item{\code{return_ft = FALSE}}{Invisibly returns a list with components:
#'       \itemize{
#'         \item \code{file} - Path to created file
#'         \item \code{caption} - Caption text (if provided)
#'       }
#'       The flextable object is accessible via \code{attr(result, "flextable")}
#'     }
#'     \item{\code{return_ft = TRUE}}{Directly returns the flextable object for 
#'       immediate further customization}
#'   }
#'   
#'   In both cases, creates a .docx file at the specified location.
#'
#' @details
#' \strong{Package Requirements:}
#' 
#' This function requires:
#' \itemize{
#'   \item \strong{flextable} - For creating formatted tables
#'   \item \strong{officer} - For Word document manipulation
#' }
#' 
#' Install if needed:
#' ```r
#' install.packages(c("flextable", "officer"))
#' ```
#' 
#' \strong{Output Features:}
#' 
#' The generated Word document contains:
#' \itemize{
#'   \item Fully editable table (native Word table, not image)
#'   \item Professional typography and spacing
#'   \item Proper page setup (size, orientation, margins)
#'   \item Caption (if provided) as separate paragraph above table
#'   \item All formatting preserved but editable
#'   \item Compatible with Word 2007 and later
#' }
#' 
#' \strong{Editability in Word:}
#' 
#' Once created, you can edit the table in Microsoft Word:
#' \itemize{
#'   \item Modify cell contents
#'   \item Adjust column widths
#'   \item Change fonts and colors
#'   \item Add/remove rows or columns
#'   \item Apply Word table styles
#'   \item Copy/paste into other documents
#'   \item Convert to plain text or other formats
#' }
#' 
#' This makes the function ideal for creating initial drafts that require 
#' manual refinement.
#' 
#' \strong{Further Customization with Flextable:}
#' 
#' For programmatic customization beyond the built-in options, access the 
#' flextable object:
#' 
#' \emph{Method 1: Via attribute (default)}
#' ```r
#' result <- tbl2docx(table, "output.docx")
#' ft <- attr(result, "flextable")
#' 
#' # Customize flextable
#' ft <- flextable::bold(ft, i = 1, j = 1, part = "body")
#' ft <- flextable::color(ft, i = 2, j = 3, color = "red")
#' 
#' # Re-save if needed
#' doc <- officer::read_docx()
#' doc <- flextable::body_add_flextable(doc, ft)
#' print(doc, target = "customized.docx")
#' ```
#' 
#' \emph{Method 2: Direct return}
#' ```r
#' ft <- tbl2docx(table, "output.docx", return_ft = TRUE)
#' 
#' # Customize immediately
#' ft <- flextable::bg(ft, bg = "yellow", part = "header")
#' ft <- flextable::autofit(ft)
#' 
#' # Save to new document
#' doc <- officer::read_docx()
#' doc <- flextable::body_add_flextable(doc, ft)
#' print(doc, target = "custom.docx")
#' ```
#' 
#' \strong{Page Layout:}
#' 
#' The function automatically sets up the Word document with:
#' \itemize{
#'   \item Specified paper size and orientation
#'   \item Standard margins (1 inch by default)
#'   \item Continuous section (no page breaks before table)
#'   \item Left-aligned table placement
#' }
#' 
#' For landscape orientation:
#' \itemize{
#'   \item Automatically swaps page width and height
#'   \item Applies landscape property to section
#'   \item Useful for wide tables with many columns
#' }
#' 
#' \strong{Table Width Management:}
#' 
#' Width behavior:
#' \itemize{
#'   \item \code{width = NULL} - Auto-fits to content and page width
#'   \item \code{width = 6} - Exactly 6 inches wide
#'   \item Width distributed evenly across columns by default
#'   \item Can adjust individual column widths in Word after creation
#' }
#' 
#' For very wide tables:
#' \enumerate{
#'   \item Use \code{orientation = "landscape"}
#'   \item Use \code{paper = "legal"} for extra width
#'   \item Reduce \code{font_size}
#'   \item Use \code{condense_table = TRUE}
#'   \item Consider breaking across multiple tables
#' }
#' 
#' \strong{Typography:}
#' 
#' The function applies professional typography:
#' \itemize{
#'   \item Column headers: Bold, slightly larger font
#'   \item Body text: Regular weight, specified font size
#'   \item Numbers: Right-aligned for easy comparison
#'   \item Text: Left-aligned for readability
#'   \item Consistent spacing: Adequate padding in cells
#' }
#' 
#' Font family must be installed on the system where Word opens the document. 
#' Common cross-platform choices:
#' \itemize{
#'   \item Arial - Sans-serif, highly readable
#'   \item Times New Roman - Serif, traditional
#'   \item Calibri - Microsoft default, modern
#'   \item Helvetica - Sans-serif, professional
#' }
#' 
#' \strong{Zebra Striping:}
#' 
#' When \code{zebra_stripes = TRUE}:
#' \itemize{
#'   \item Alternating variables receive light gray background
#'   \item All rows of same variable share same shading
#'   \item Improves visual grouping
#'   \item Particularly useful for tables with many factor variables
#'   \item Color can be changed in Word after creation
#' }
#' 
#' \strong{Dark Header:}
#' 
#' When \code{dark_header = TRUE}:
#' \itemize{
#'   \item Header row: Dark gray/black background
#'   \item Header text: White for high contrast
#'   \item Modern, professional appearance
#'   \item Draws attention to column names
#' }
#' 
#' \strong{Integration with R Markdown/Quarto:}
#' 
#' For R Markdown/Quarto Word output:
#' ```r
#' # Create flextable for inline display
#' ft <- tbl2docx(results, "temp.docx", return_ft = TRUE)
#' 
#' # Display in R Markdown chunk
#' ft  # Renders in Word output
#' ```
#' 
#' Or use flextable directly in chunks:
#' ````markdown
#' ```{r}
#' flextable::flextable(results)
#' ```
#' ````
#' 
#' \strong{Optimized for FastFit Tables:}
#' 
#' Specifically designed for tables from:
#' \itemize{
#'   \item \code{\link{desctbl}} - Descriptive statistics
#'   \item \code{\link{fit}} - Regression results
#'   \item \code{\link{uscreen}} - Univariable screening
#'   \item \code{\link{fastfit}} - Combined analyses
#'   \item \code{\link{compfit}} - Model comparisons
#' }
#' 
#' Automatic handling of:
#' \itemize{
#'   \item Sample size rows (N = X)
#'   \item Variable grouping
#'   \item P-value formatting and bolding
#'   \item Confidence intervals
#'   \item Multi-level categorical variables
#' }
#' 
#' \strong{Troubleshooting:}
#' 
#' Common issues and solutions:
#' 
#' \emph{Table too wide:}
#' \itemize{
#'   \item Use \code{orientation = "landscape"}
#'   \item Reduce \code{font_size}
#'   \item Use \code{condense_table = TRUE}
#'   \item Manually adjust column widths in Word
#' }
#' 
#' \emph{Font not displaying correctly:}
#' \itemize{
#'   \item Ensure font is installed on system
#'   \item Use common fonts (Arial, Times, Calibri)
#'   \item Check font name spelling
#' }
#' 
#' \emph{Need more customization:}
#' \itemize{
#'   \item Use \code{return_ft = TRUE}
#'   \item Modify flextable object with flextable package functions
#'   \item Edit directly in Word after creation
#' }
#'
#' @seealso
#' \code{\link{tbl2pptx}} for PowerPoint slides,
#' \code{\link{tbl2pdf}} for PDF output,
#' \code{\link{tbl2html}} for HTML tables,
#' \code{\link{tbl2tex}} for LaTeX output,
#' \code{\link[flextable]{flextable}} for the underlying table object,
#' \code{\link[officer]{read_docx}} for Word document manipulation
#'
#' @examples
#' \dontrun{
#' # Load data
#' data(clintrial)
#' data(clintrial_labels)
#' 
#' # Create regression table
#' results <- fit(
#'     data = clintrial,
#'     outcome = "os_status",
#'     predictors = c("age", "sex", "treatment", "stage"),
#'     var_labels = clintrial_labels
#' )
#' 
#' # Example 1: Basic Word export
#' tbl2docx(results, "results.docx")
#' 
#' # Example 2: With caption
#' tbl2docx(results, "captioned.docx",
#'          caption = "Table 1: Multivariable Logistic Regression Results")
#' 
#' # Example 3: Landscape orientation for wide tables
#' tbl2docx(results, "wide.docx",
#'          orientation = "landscape")
#' 
#' # Example 4: Custom font and size
#' tbl2docx(results, "custom_font.docx",
#'          font_family = "Times New Roman",
#'          font_size = 11)
#' 
#' # Example 5: Hierarchical display
#' tbl2docx(results, "indented.docx",
#'          indent_groups = TRUE)
#' 
#' # Example 6: Condensed table
#' tbl2docx(results, "condensed.docx",
#'          condense_table = TRUE)
#' 
#' # Example 7: With zebra stripes
#' tbl2docx(results, "striped.docx",
#'          zebra_stripes = TRUE)
#' 
#' # Example 8: Dark header style
#' tbl2docx(results, "dark.docx",
#'          dark_header = TRUE)
#' 
#' # Example 9: A4 paper for international journals
#' tbl2docx(results, "a4.docx",
#'          paper = "a4")
#' 
#' # Example 10: Get flextable for customization
#' result <- tbl2docx(results, "base.docx")
#' ft <- attr(result, "flextable")
#' 
#' # Customize the flextable
#' ft <- flextable::bold(ft, i = 1, part = "body")
#' ft <- flextable::color(ft, j = "p-value", color = "blue")
#' 
#' # Example 11: Direct flextable return
#' ft <- tbl2docx(results, "direct.docx", return_ft = TRUE)
#' ft <- flextable::bg(ft, bg = "yellow", part = "header")
#' 
#' # Example 12: Publication-ready table
#' tbl2docx(results, "publication.docx",
#'          caption = "Table 2: Adjusted Odds Ratios for Mortality",
#'          font_family = "Times New Roman",
#'          font_size = 10,
#'          indent_groups = TRUE,
#'          zebra_stripes = FALSE,
#'          bold_significant = TRUE)
#' 
#' # Example 13: Custom column alignment
#' tbl2docx(results, "aligned.docx",
#'          align = c("left", "left", "center", "right", "right"))
#' 
#' # Example 14: Disable significance bolding
#' tbl2docx(results, "no_bold.docx",
#'          bold_significant = FALSE)
#' 
#' # Example 15: Stricter significance threshold
#' tbl2docx(results, "strict.docx",
#'          bold_significant = TRUE,
#'          sig_threshold = 0.01)
#' }
#'
#' @export
tbl2docx <- function(table,
                     file,
                     caption = NULL,
                     font_size = 8,
                     font_family = "Arial",
                     format_headers = TRUE,
                     bold_significant = TRUE,
                     sig_threshold = 0.05,
                     indent_groups = FALSE,
                     condense_table = FALSE,
                     zebra_stripes = FALSE,
                     dark_header = FALSE,
                     paper = "letter",
                     orientation = "portrait",
                     width = NULL,
                     align = NULL,
                     return_ft = FALSE,
                     ...) {
    
    if (!requireNamespace("flextable", quietly = TRUE)) {
        stop("Package 'flextable' required. Install with: install.packages('flextable')")
    }
    if (!requireNamespace("officer", quietly = TRUE)) {
        stop("Package 'officer' required. Install with: install.packages('officer')")
    }
    
    if (!grepl("\\.docx$", tolower(file))) {
        stop("File must have .docx extension")
    }
    
                                        # Validate paper and orientation
    paper <- match.arg(paper, c("letter", "a4", "legal"))
    orientation <- match.arg(orientation, c("portrait", "landscape"))
    
                                        # Process table using shared function
    result <- process_table_for_flextable(
        table = table,
        caption = caption,
        font_size = font_size,
        font_family = font_family,
        format_headers = format_headers,
        bold_significant = bold_significant,
        sig_threshold = sig_threshold,
        indent_groups = indent_groups,
        condense_table = condense_table,
        zebra_stripes = zebra_stripes,
        dark_header = dark_header,
        paper = paper,
        orientation = orientation,
        width = width,
        align = align
    )
    
                                        # Extract flextable for easy access
    ft <- result$ft
    
                                        # Create Word document
    doc <- officer::read_docx(...)
    
                                        # Set page size and orientation using prop_section
    page_width <- if(paper == "a4") 8.27 else 8.5
    page_height <- if(paper == "a4") 11.69 else if(paper == "legal") 14 else 11
    
    if (orientation == "landscape") {
        temp <- page_width
        page_width <- page_height
        page_height <- temp
    }
    
    sect_properties <- officer::prop_section(
                                    page_size = officer::page_size(width = page_width, 
                                                                   height = page_height, 
                                                                   orient = orientation),
                                    type = "continuous",
                                    page_margins = officer::page_mar()
                                )
    
    doc <- officer::body_set_default_section(doc, value = sect_properties)
    
                                        # Add caption if provided
    if (!is.null(result$caption)) {
        doc <- officer::body_add_par(doc, result$caption, 
                                     style = "Normal", pos = "after")
        doc <- officer::body_add_par(doc, "", style = "Normal")
    }
    
                                        # Add table
    doc <- flextable::body_add_flextable(doc, ft, align = "left")
    
                                        # Save document
    print(doc, target = file)
    
    message(sprintf("Table exported to %s", file))
    
                                        # Return based on user preference
    if (return_ft) {
                                        # Return flextable directly for easy access
        return(ft)
    } else {
                                        # FIX: Create a proper list to hold attributes
                                        # Return invisibly with flextable as attribute for backward compatibility
        result_obj <- list(
            file = file,
            caption = result$caption
        )
        attr(result_obj, "flextable") <- ft
        class(result_obj) <- c("tbl2docx_result", "list")
        return(invisible(result_obj))
    }
}

#' Print method for tbl2docx results
#' @export
print.tbl2docx_result <- function(x, ...) {
    cat("Table exported to:", x$file, "\n")
    if (!is.null(x$caption)) {
        cat("Caption:", x$caption, "\n")
    }
    cat("Flextable object available via: attr(result, 'flextable')\n")
    invisible(x)
}
