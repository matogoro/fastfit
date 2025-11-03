#' Export Table to Microsoft PowerPoint Format (PPTX)
#'
#' Converts a data frame, data.table, or matrix to a Microsoft PowerPoint slide 
#' (.pptx) with a formatted table using the flextable and officer packages. Creates 
#' presentation-ready slides with extensive control over table formatting, positioning, 
#' and layout. Tables can be further edited in PowerPoint after creation. Ideal for 
#' creating data-driven presentations, conference talks, and stakeholder reports.
#'
#' @param table A data.frame, data.table, or matrix to export. Can be output from 
#'   \code{\link{desctbl}}, \code{\link{fit}}, \code{\link{uscreen}}, 
#'   \code{\link{fastfit}}, \code{\link{compfit}}, or any tabular data.
#'   
#' @param file Character string specifying the output PPTX filename. Must have 
#'   \code{.pptx} extension. Example: \code{"results.pptx"}, \code{"Slide1.pptx"}.
#'   
#' @param caption Character string. Optional title displayed in the slide's title 
#'   placeholder or as text box above the table. Default is \code{NULL}.
#'   
#' @param font_size Numeric. Base font size in points for table content. 
#'   Default is 10. Typical range for presentations: 10-14 points. Larger than 
#'   print documents for visibility at distance.
#'   
#' @param font_family Character string. Font family name for the table. Must be 
#'   installed on the system. Default is \code{"Arial"}. Common presentation 
#'   fonts: \code{"Calibri"}, \code{"Helvetica"}, \code{"Arial"}.
#'   
#' @param format_headers Logical. If \code{TRUE}, formats column headers by 
#'   italicizing statistical notation (\emph{n}, \emph{p}) and improving readability. 
#'   Default is \code{TRUE}.
#'   
#' @param bold_significant Logical. If \code{TRUE}, applies bold formatting to 
#'   p-values below the significance threshold. Makes important results stand 
#'   out in presentations. Default is \code{TRUE}.
#'   
#' @param sig_threshold Numeric. P-value threshold for bold formatting. 
#'   Default is 0.05.
#'   
#' @param indent_groups Logical. If \code{TRUE}, indents factor levels under 
#'   their parent variable, creating hierarchical display. Useful for categorical 
#'   variables. Default is \code{FALSE}.
#'   
#' @param condense_table Logical. If \code{TRUE}, condenses table to essential 
#'   rows only. Automatically sets \code{indent_groups = TRUE}. Crucial for 
#'   fitting content on slides. Default is \code{FALSE}.
#'   
#' @param zebra_stripes Logical. If \code{TRUE}, applies alternating row shading 
#'   to different variables for visual grouping. Improves readability during 
#'   presentations. Default is \code{FALSE}.
#'   
#' @param dark_header Logical. If \code{TRUE}, creates dark background with 
#'   light text for header row. Provides strong contrast visible from distance. 
#'   Default is \code{FALSE}.
#'   
#' @param width Numeric. Table width in inches. If \code{NULL} (default), 
#'   auto-fits to slide width (approximately 9 inches for standard 10-inch slide 
#'   with margins). Specify for exact control.
#'   
#' @param align Character vector specifying column alignment. Options: 
#'   \code{"left"}, \code{"center"}, or \code{"right"}. If \code{NULL} (default), 
#'   automatically determines based on content.
#'   
#' @param template Character string. Path to custom PPTX template file. If 
#'   \code{NULL} (default), uses officer's default blank template. Use to match 
#'   corporate branding or conference themes.
#'   
#' @param layout Character string. Name of slide layout to use from template. 
#'   Default is \code{"Title and Content"}. Common layouts:
#'   \itemize{
#'     \item \code{"Title and Content"} - Standard layout [default]
#'     \item \code{"Blank"} - Empty slide for maximum control
#'     \item \code{"Title Only"} - Title area only
#'     \item \code{"Two Content"} - Title with two content areas
#'   }
#'   
#' @param master Character string. Name of slide master to use. Default is 
#'   \code{"Office Theme"}. Varies by template. Check template for available masters.
#'   
#' @param left Numeric. Horizontal position from left edge of slide in inches. 
#'   Default is 0.5. Standard slide is 10 inches wide.
#'   
#' @param top Numeric. Vertical position from top edge of slide in inches. 
#'   Default is 1.5 (leaves room for title). Standard slide is 7.5 inches tall. 
#'   Adjust based on table size and layout.
#'   
#' @param return_ft Logical. If \code{TRUE}, returns the flextable object 
#'   directly for further customization. If \code{FALSE} (default), returns 
#'   invisibly with metadata. Default is \code{FALSE}.
#'   
#' @param ... Additional arguments passed to \code{\link[officer]{read_pptx}}.
#'
#' @return Behavior depends on \code{return_ft}:
#'   \describe{
#'     \item{\code{return_ft = FALSE}}{Invisibly returns a list with:
#'       \itemize{
#'         \item \code{file} - Path to created file
#'         \item \code{caption} - Caption/title text
#'         \item \code{layout} - Layout name used
#'         \item \code{master} - Master name used
#'         \item \code{template} - Template path (if provided)
#'         \item \code{position} - List with \code{left} and \code{top} coordinates
#'       }
#'       Flextable accessible via \code{attr(result, "flextable")}
#'     }
#'     \item{\code{return_ft = TRUE}}{Directly returns the flextable object}
#'   }
#'   
#'   Always creates a .pptx file at the specified location.
#'
#' @details
#' \strong{Package Requirements:}
#' 
#' Requires:
#' \itemize{
#'   \item \strong{flextable} - Table creation and formatting
#'   \item \strong{officer} - PowerPoint manipulation
#' }
#' 
#' Install: \code{install.packages(c("flextable", "officer"))}
#' 
#' \strong{Presentation Considerations:}
#' 
#' Tables in presentations differ from print:
#' \itemize{
#'   \item \strong{Larger fonts} - Visible from back of room (10-14pt minimum)
#'   \item \strong{Less content} - Limit rows to fit slide (use \code{condense_table})
#'   \item \strong{High contrast} - Consider \code{dark_header = TRUE}
#'   \item \strong{Simple layout} - Avoid excessive detail
#'   \item \strong{Key findings} - Show only essential results
#' }
#' 
#' \strong{Slide Dimensions:}
#' 
#' Standard PowerPoint slide:
#' \itemize{
#'   \item Width: 10 inches (25.4 cm)
#'   \item Height: 7.5 inches (19.05 cm)
#'   \item Aspect ratio: 4:3 (standard) or 16:9 (widescreen)
#' }
#' 
#' Safe content area (with margins):
#' \itemize{
#'   \item Width: ~9 inches
#'   \item Height: ~6 inches (accounting for title)
#' }
#' 
#' \strong{Positioning:}
#' 
#' The \code{left} and \code{top} parameters control table placement:
#' \itemize{
#'   \item (0, 0) = Top-left corner of slide
#'   \item Default (0.5, 1.5) = Standard position with title room
#'   \item Center: \code{left = (10 - table_width) / 2}
#' }
#' 
#' When caption is provided:
#' \itemize{
#'   \item Attempts to use title placeholder (if layout supports)
#'   \item Falls back to text box above table
#'   \item Automatically adjusts table position downward
#' }
#' 
#' \strong{Slide Layouts:}
#' 
#' Different layouts serve different purposes:
#' 
#' \emph{Title and Content (default):}
#' \itemize{
#'   \item Has title and content placeholders
#'   \item Caption goes in title area
#'   \item Table in content area
#'   \item Most common for data slides
#' }
#' 
#' \emph{Blank:}
#' \itemize{
#'   \item No predefined areas
#'   \item Maximum flexibility
#'   \item Use absolute positioning (\code{left}, \code{top})
#'   \item Good for custom layouts
#' }
#' 
#' \emph{Title Only:}
#' \itemize{
#'   \item Title area only
#'   \item Large space for table
#'   \item Good for data-heavy slides
#' }
#' 
#' \strong{Custom Templates:}
#' 
#' Use organizational or conference templates:
#' ```r
#' tbl2pptx(table, "branded.pptx",
#'          template = "company_template.pptx",
#'          layout = "Content Layout",  # Name from template
#'          master = "Company Theme")   # Name from template
#' ```
#' 
#' To find layout and master names in template:
#' ```r
#' pres <- officer::read_pptx("template.pptx")
#' officer::layout_summary(pres)
#' ```
#' 
#' \strong{Table Size Management:}
#' 
#' Fitting tables on slides:
#' \enumerate{
#'   \item \strong{Condense content}: \code{condense_table = TRUE}
#'   \item \strong{Reduce font}: \code{font_size = 9} or smaller
#'   \item \strong{Limit rows}: Filter data before export
#'   \item \strong{Limit columns}: Select key variables only
#'   \item \strong{Split tables}: Multiple slides for large datasets
#'   \item \strong{Adjust width}: \code{width = 9} for maximum
#' }
#' 
#' \strong{Visual Design Tips:}
#' 
#' For effective presentation tables:
#' \itemize{
#'   \item Use \code{zebra_stripes = TRUE} to group related rows
#'   \item Use \code{dark_header = TRUE} for emphasis
#'   \item Bold significant results (\code{bold_significant = TRUE})
#'   \item Limit to 5-8 rows maximum for readability
#'   \item Use \code{indent_groups = TRUE} for hierarchy
#'   \item Ensure sufficient contrast for projectors
#'   \item Test visibility from distance
#' }
#' 
#' \strong{Multiple Slides:}
#' 
#' Creating presentations with multiple tables:
#' ```r
#' # Each call creates new presentation - combine after
#' tbl2pptx(table1, "slide1.pptx", caption = "Results Part 1")
#' tbl2pptx(table2, "slide2.pptx", caption = "Results Part 2")
#' 
#' # Then manually combine in PowerPoint, or:
#' # Use officer to create multi-slide presentation
#' pres <- officer::read_pptx()
#' 
#' # Add first table
#' ft1 <- tbl2pptx(table1, "temp1.pptx", return_ft = TRUE)
#' pres <- officer::add_slide(pres)
#' pres <- officer::ph_with(pres, ft1, 
#'                          location = officer::ph_location(left = 0.5, top = 1.5))
#' 
#' # Add second table
#' ft2 <- tbl2pptx(table2, "temp2.pptx", return_ft = TRUE)
#' pres <- officer::add_slide(pres)
#' pres <- officer::ph_with(pres, ft2,
#'                          location = officer::ph_location(left = 0.5, top = 1.5))
#' 
#' print(pres, target = "combined.pptx")
#' ```
#' 
#' \strong{Further Customization:}
#' 
#' Access flextable for advanced formatting:
#' ```r
#' ft <- tbl2pptx(table, "base.pptx", return_ft = TRUE)
#' 
#' # Customize
#' ft <- flextable::color(ft, j = "p-value", color = "red")
#' ft <- flextable::bg(ft, i = 1, bg = "yellow")
#' ft <- flextable::bold(ft, i = ~ estimate > 0, j = "estimate")
#' 
#' # Save to new slide
#' pres <- officer::read_pptx()
#' pres <- officer::add_slide(pres)
#' pres <- officer::ph_with(pres, ft, 
#'                          location = officer::ph_location(left = 0.5, top = 1.5))
#' print(pres, target = "custom.pptx")
#' ```
#' 
#' \strong{Editability in PowerPoint:}
#' 
#' Created tables are fully editable:
#' \itemize{
#'   \item Modify cell contents
#'   \item Adjust column widths
#'   \item Change colors and fonts
#'   \item Add/remove rows or columns
#'   \item Reposition on slide
#'   \item Copy to other slides
#'   \item Apply PowerPoint table styles
#' }
#' 
#' \strong{Best Practices:}
#' 
#' For professional presentations:
#' \enumerate{
#'   \item Keep tables simple (5-8 rows maximum)
#'   \item Use large fonts (≥10pt)
#'   \item Highlight key findings with bold or color
#'   \item Add clear, concise titles
#'   \item Test visibility on projector
#'   \item Consider audience distance
#'   \item Use consistent styling across slides
#'   \item Provide slide notes with details
#' }
#'
#' @seealso
#' \code{\link{tbl2docx}} for Word documents,
#' \code{\link{tbl2pdf}} for PDF output,
#' \code{\link{tbl2html}} for HTML tables,
#' \code{\link{tbl2tex}} for LaTeX output,
#' \code{\link[flextable]{flextable}} for table customization,
#' \code{\link[officer]{read_pptx}} for PowerPoint manipulation
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
#'     predictors = c("age", "sex", "treatment"),
#'     var_labels = clintrial_labels
#' )
#' 
#' # Example 1: Basic PowerPoint slide
#' tbl2pptx(results, "results.pptx")
#' 
#' # Example 2: With title
#' tbl2pptx(results, "titled.pptx",
#'          caption = "Multivariable Regression Results")
#' 
#' # Example 3: Larger font for visibility
#' tbl2pptx(results, "large_font.pptx",
#'          font_size = 12,
#'          caption = "Main Findings")
#' 
#' # Example 4: Condensed for slide space
#' tbl2pptx(results, "condensed.pptx",
#'          condense_table = TRUE,
#'          caption = "Key Results")
#' 
#' # Example 5: Dark header for emphasis
#' tbl2pptx(results, "dark.pptx",
#'          dark_header = TRUE,
#'          caption = "Risk Factors")
#' 
#' # Example 6: With zebra stripes
#' tbl2pptx(results, "striped.pptx",
#'          zebra_stripes = TRUE)
#' 
#' # Example 7: Blank layout with custom positioning
#' tbl2pptx(results, "blank.pptx",
#'          layout = "Blank",
#'          left = 1,
#'          top = 1.5,
#'          width = 8)
#' 
#' # Example 8: Using company template
#' tbl2pptx(results, "branded.pptx",
#'          template = "company_template.pptx",
#'          layout = "Data Slide",
#'          caption = "Q4 Results")
#' 
#' # Example 9: Get flextable for customization
#' ft <- tbl2pptx(results, "base.pptx", return_ft = TRUE)
#' 
#' # Customize colors
#' ft <- flextable::color(ft, j = "p-value", color = "darkred")
#' ft <- flextable::bold(ft, i = ~ p_value < 0.01)
#' 
#' # Example 10: Presentation-optimized table
#' tbl2pptx(results, "presentation.pptx",
#'          caption = "Main Analysis Results",
#'          font_size = 11,
#'          condense_table = TRUE,
#'          zebra_stripes = TRUE,
#'          dark_header = TRUE,
#'          bold_significant = TRUE)
#' 
#' # Example 11: Descriptive statistics slide
#' desc <- desctbl(
#'     data = clintrial,
#'     strata = "treatment",
#'     vars = c("age", "sex", "bmi"),
#'     var_labels = clintrial_labels
#' )
#' 
#' tbl2pptx(desc, "baseline.pptx",
#'          caption = "Baseline Characteristics",
#'          font_size = 10)
#' 
#' # Example 12: Conference presentation style
#' tbl2pptx(results, "conference.pptx",
#'          caption = "Study Outcomes",
#'          font_family = "Calibri",
#'          font_size = 14,  # Large for big rooms
#'          dark_header = TRUE,
#'          condense_table = TRUE)
#' }
#'
#' @export
tbl2pptx <- function(table,
                     file,
                     caption = NULL,
                     font_size = 10,
                     font_family = "Arial",
                     format_headers = TRUE,
                     bold_significant = TRUE,
                     sig_threshold = 0.05,
                     indent_groups = FALSE,
                     condense_table = FALSE,
                     zebra_stripes = FALSE,
                     dark_header = FALSE,
                     width = NULL,
                     align = NULL,
                     template = NULL,
                     layout = "Title and Content",
                     master = "Office Theme",
                     left = 0.5,
                     top = 1.5,
                     return_ft = FALSE,
                     ...) {
    
    if (!requireNamespace("flextable", quietly = TRUE)) {
        stop("Package 'flextable' required. Install with: install.packages('flextable')")
    }
    if (!requireNamespace("officer", quietly = TRUE)) {
        stop("Package 'officer' required. Install with: install.packages('officer')")
    }
    
    if (!grepl("\\.pptx$", tolower(file))) {
        stop("File must have .pptx extension")
    }
    
                                        # Process table using shared function
                                        # For PowerPoint, default to landscape-style dimensions
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
        paper = "letter",
        orientation = "landscape",  # Usually better for slides
        width = width,
        align = align
    )
    
                                        # Extract flextable
    ft <- result$ft
    
                                        # Adjust width for PowerPoint if not specified
    if (is.null(width)) {
                                        # Standard slide is 10 inches wide, leave margins
        slide_width <- 9  # Leave 0.5 inch margins on each side
        ft <- flextable::width(ft, width = slide_width / ncol(table))
    }
    
                                        # Create or read PowerPoint presentation
    if (!is.null(template) && file.exists(template)) {
        pres <- officer::read_pptx(template, ...)
    } else {
        pres <- officer::read_pptx(...)
    }
    
                                        # Add a new slide
    pres <- officer::add_slide(pres, layout = layout, master = master)
    
                                        # Track actual position used (may be adjusted if caption is added)
    actual_top <- top
    
                                        # Add title if provided and layout supports it
    if (!is.null(caption)) {
        if (layout %in% c("Title and Content", "Title Only", "Two Content")) {
            tryCatch({
                pres <- officer::ph_with(pres, value = caption, 
                                         location = officer::ph_location_type(type = "title"))
            }, error = function(e) {
                                        # If title placeholder not available, add as text box
                pres <- officer::ph_with(pres, value = caption,
                                         location = officer::ph_location(left = left, top = top - 0.5,
                                                                         width = 9, height = 0.5))
                                        # Adjust table position down
                actual_top <- top + 0.5
            })
        } else {
                                        # For layouts without title placeholder, add as text box above table
            pres <- officer::ph_with(pres, value = caption,
                                     location = officer::ph_location(left = left, top = top - 0.5,
                                                                     width = 9, height = 0.5))
                                        # Adjust table position down
            actual_top <- top + 0.5
        }
    }
    
                                        # Add the table to the slide
    if (layout == "Title and Content" && is.null(caption)) {
                                        # Try to use content placeholder if available
        tryCatch({
            pres <- officer::ph_with(pres, value = ft,
                                     location = officer::ph_location_type(type = "body"))
        }, error = function(e) {
                                        # Fall back to absolute positioning
            pres <- officer::ph_with(pres, value = ft,
                                     location = officer::ph_location(left = left, top = actual_top))
        })
    } else {
                                        # Use absolute positioning
        pres <- officer::ph_with(pres, value = ft,
                                 location = officer::ph_location(left = left, top = actual_top))
    }
    
                                        # Save presentation
    print(pres, target = file)
    
    message(sprintf("Table exported to %s", file))
    
                                        # Return based on user preference
    if (return_ft) {
                                        # Return flextable directly for easy access
        return(ft)
    } else {
                                        # FIX: Create a proper list to hold attributes (same as tbl2docx)
                                        # Return invisibly with flextable and metadata as attributes
        result_obj <- list(
            file = file,
            caption = caption,
            layout = layout,
            master = master,
            template = template,
            position = list(left = left, top = actual_top)
        )
        attr(result_obj, "flextable") <- ft
        class(result_obj) <- c("tbl2pptx_result", "list")
        return(invisible(result_obj))
    }
}

#' Print method for tbl2pptx results
#' @keywords internal
#' @export
print.tbl2pptx_result <- function(x, ...) {
    cat("Table exported to:", x$file, "\n")
    if (!is.null(x$caption)) {
        cat("Caption:", x$caption, "\n")
    }
    cat("Layout:", x$layout, "\n")
    cat("Master:", x$master, "\n")
    if (!is.null(x$template)) {
        cat("Template:", x$template, "\n")
    }
    cat("Position: left =", x$position$left, ", top =", x$position$top, "\n")
    cat("Flextable object available via: attr(result, 'flextable')\n")
    invisible(x)
}
