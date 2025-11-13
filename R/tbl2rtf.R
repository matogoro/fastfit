#' Export Table to Rich Text Format (RTF)
#'
#' Converts a data frame, data.table, or matrix to a Rich Text Format (.rtf) 
#' document using the flextable and officer packages. Creates professional, 
#' widely-compatible tables with extensive formatting options. RTF files can be 
#' opened and edited in Microsoft Word, LibreOffice, WordPad, and many other word 
#' processors. Particularly useful for regulatory submissions, cross-platform 
#' compatibility, and when maximum editability is required.
#'
#' @param table A data.frame, data.table, or matrix to export. Can be output from 
#'   \code{\link{desctbl}}, \code{\link{fit}}, \code{\link{uscreen}}, 
#'   \code{\link{fastfit}}, \code{\link{compfit}}, or any tabular data.
#'   
#' @param file Character string specifying the output RTF filename. Must have 
#'   \code{.rtf} extension. Example: \code{"results.rtf"}, \code{"Table1.rtf"}.
#'   
#' @param caption Character string. Optional caption displayed above the table 
#'   in the RTF document. Default is \code{NULL}.
#'   
#' @param font_size Numeric. Base font size in points for table content. 
#'   Default is 8. Typical range: 8-12 points. Headers use slightly larger size.
#'   
#' @param font_family Character string. Font family name for the table. Must be 
#'   a font installed on the system. Default is \code{"Arial"}. Common options: 
#'   \code{"Times New Roman"}, \code{"Calibri"}, \code{"Courier New"}.
#'   
#' @param format_headers Logical. If \code{TRUE}, formats column headers by 
#'   italicizing statistical notation (\emph{n}, \emph{p}), converting underscores 
#'   to spaces, and improving readability. Default is \code{TRUE}.
#'   
#' @param bold_significant Logical. If \code{TRUE}, applies bold formatting to 
#'   p-values below the significance threshold. Makes significant results stand 
#'   out. Default is \code{TRUE}.
#'   
#' @param p_threshold Numeric. P-value threshold for bold formatting. Only 
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
#' @param ... Additional arguments (currently unused, reserved for future extensions).
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
#'   In both cases, creates a .rtf file at the specified location.
#'
#' @details
#' \strong{Package Requirements:}
#' 
#' This function requires:
#' \itemize{
#'   \item \strong{flextable} - For creating formatted tables
#'   \item \strong{officer} - For RTF document generation
#' }
#' 
#' Install if needed:
#' ```r
#' install.packages(c("flextable", "officer"))
#' ```
#' 
#' \strong{RTF Format Advantages:}
#' 
#' RTF (Rich Text Format) is a universal document format with several advantages:
#' \itemize{
#'   \item \strong{Maximum compatibility} - Opens in virtually all word processors
#'   \item \strong{Cross-platform} - Works on Windows, Mac, Linux without conversion
#'   \item \strong{Fully editable} - Native text format, not embedded objects
#'   \item \strong{Lightweight} - Smaller file sizes than DOCX
#'   \item \strong{Regulatory compliance} - Widely accepted for submissions (FDA, EMA)
#'   \item \strong{Long-term accessibility} - Simple text-based format
#'   \item \strong{Version control friendly} - Text-based, works with diff tools
#' }
#' 
#' Applications that can open RTF files:
#' \itemize{
#'   \item Microsoft Word (Windows, Mac)
#'   \item LibreOffice Writer
#'   \item Apache OpenOffice Writer
#'   \item WordPad (Windows built-in)
#'   \item TextEdit (Mac built-in)
#'   \item Google Docs (with import)
#'   \item Pages (Mac)
#'   \item Many other word processors
#' }
#' 
#' \strong{Regulatory Submissions:}
#' 
#' RTF is commonly used for regulatory submissions because:
#' \itemize{
#'   \item Accepted by FDA, EMA, PMDA, and other agencies
#'   \item Format specification is open and stable
#'   \item No proprietary dependencies
#'   \item Tables remain editable for reviewers
#'   \item Consistent rendering across platforms
#'   \item Can be validated programmatically
#' }
#' 
#' For clinical trial reporting (ICH E3 guidelines), RTF tables are often preferred 
#' over PDF because reviewers can extract data, add comments, and compare versions.
#' 
#' \strong{Output Features:}
#' 
#' The generated RTF document contains:
#' \itemize{
#'   \item Fully editable table (native RTF table, not image)
#'   \item Professional typography and spacing
#'   \item Proper page setup (size, orientation, margins)
#'   \item Caption (if provided) as separate paragraph above table
#'   \item All formatting preserved but editable
#'   \item Compatible with RTF 1.5 specification
#' }
#' 
#' \strong{Editability:}
#' 
#' Once created, you can edit the RTF table in any word processor:
#' \itemize{
#'   \item Modify cell contents
#'   \item Adjust column widths
#'   \item Change fonts and colors
#'   \item Add/remove rows or columns
#'   \item Apply different table styles
#'   \item Copy/paste into other documents
#'   \item Convert to other formats (DOCX, PDF, HTML)
#' }
#' 
#' This makes RTF ideal for creating tables that require review, approval, or 
#' modification by multiple stakeholders across different platforms.
#' 
#' \strong{Further Customization with Flextable:}
#' 
#' For programmatic customization beyond the built-in options, access the 
#' flextable object:
#' 
#' \emph{Method 1: Via attribute (default)}
#' ```r
#' result <- tbl2rtf(table, "output.rtf")
#' ft <- attr(result, "flextable")
#' 
#' # Customize flextable
#' ft <- flextable::bold(ft, i = 1, j = 1, part = "body")
#' ft <- flextable::color(ft, i = 2, j = 3, color = "red")
#' 
#' # Re-save if needed
#' flextable::save_as_rtf(ft, path = "customized.rtf")
#' ```
#' 
#' \emph{Method 2: Direct return}
#' ```r
#' ft <- tbl2rtf(table, "output.rtf", return_ft = TRUE)
#' 
#' # Customize immediately
#' ft <- flextable::bg(ft, bg = "yellow", part = "header")
#' ft <- flextable::autofit(ft)
#' 
#' # Save to new file
#' flextable::save_as_rtf(ft, path = "custom.rtf")
#' ```
#' 
#' \strong{Page Layout:}
#' 
#' The function automatically sets up the RTF document with:
#' \itemize{
#'   \item Specified paper size and orientation
#'   \item Standard margins (1 inch by default)
#'   \item Table positioned at document start
#'   \item Left-aligned table placement
#' }
#' 
#' For landscape orientation:
#' \itemize{
#'   \item Automatically swaps page dimensions
#'   \item Applies landscape property
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
#'   \item Can adjust individual column widths in word processor after creation
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
#'   \item Statistical notation: Italicized appropriately
#' }
#' 
#' \strong{Format Comparison:}
#' 
#' Choose the right export format for your needs:
#' 
#' \strong{RTF (\code{tbl2rtf}):}
#' \itemize{
#'   \item \strong{Best for:} Regulatory submissions, cross-platform compatibility
#'   \item \strong{Pros:} Universal compatibility, fully editable, lightweight
#'   \item \strong{Cons:} Limited advanced formatting compared to DOCX
#' }
#' 
#' \strong{DOCX (\code{tbl2docx}):}
#' \itemize{
#'   \item \strong{Best for:} Modern Word workflows, internal documents
#'   \item \strong{Pros:} Rich formatting, native Word format
#'   \item \strong{Cons:} Requires Office/compatible software
#' }
#' 
#' \strong{PDF (\code{tbl2pdf}):}
#' \itemize{
#'   \item \strong{Best for:} Final publications, fixed layout
#'   \item \strong{Pros:} Professional appearance, consistent rendering
#'   \item \strong{Cons:} Not editable without special tools
#' }
#' 
#' \strong{PPTX (\code{tbl2pptx}):}
#' \itemize{
#'   \item \strong{Best for:} Presentations, talks, meetings
#'   \item \strong{Pros:} Designed for projection, large fonts
#'   \item \strong{Cons:} Limited page layout options
#' }
#'
#' @seealso
#' \code{\link{tbl2docx}} for Word documents,
#' \code{\link{tbl2pptx}} for PowerPoint slides,
#' \code{\link{tbl2pdf}} for PDF output,
#' \code{\link{tbl2html}} for HTML tables,
#' \code{\link{tbl2tex}} for LaTeX output,
#' \code{\link[flextable]{flextable}} for the underlying table object,
#' \code{\link[flextable]{save_as_rtf}} for direct RTF export
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
#'     labels = clintrial_labels
#' )
#' 
#' # Example 1: Basic RTF export
#' tbl2rtf(results, "results.rtf")
#' 
#' # Example 2: With caption
#' tbl2rtf(results, "captioned.rtf",
#'         caption = "Table 1: Multivariable Logistic Regression Results")
#' 
#' # Example 3: Landscape orientation for wide tables
#' tbl2rtf(results, "wide.rtf",
#'         orientation = "landscape")
#' 
#' # Example 4: Custom font and size
#' tbl2rtf(results, "custom_font.rtf",
#'         font_family = "Times New Roman",
#'         font_size = 11)
#' 
#' # Example 5: Hierarchical display
#' tbl2rtf(results, "indented.rtf",
#'         indent_groups = TRUE)
#' 
#' # Example 6: Condensed table
#' tbl2rtf(results, "condensed.rtf",
#'         condense_table = TRUE)
#' 
#' # Example 7: With zebra stripes
#' tbl2rtf(results, "striped.rtf",
#'         zebra_stripes = TRUE)
#' 
#' # Example 8: Dark header style
#' tbl2rtf(results, "dark.rtf",
#'         dark_header = TRUE)
#' 
#' # Example 9: A4 paper for international submissions
#' tbl2rtf(results, "a4.rtf",
#'         paper = "a4")
#' 
#' # Example 10: Get flextable for customization
#' result <- tbl2rtf(results, "base.rtf")
#' ft <- attr(result, "flextable")
#' 
#' # Customize the flextable
#' ft <- flextable::bold(ft, i = 1, part = "body")
#' ft <- flextable::color(ft, j = "p-value", color = "blue")
#' 
#' # Re-save
#' flextable::save_as_rtf(ft, path = "customized.rtf")
#' 
#' # Example 11: Direct flextable return
#' ft <- tbl2rtf(results, "direct.rtf", return_ft = TRUE)
#' ft <- flextable::bg(ft, bg = "yellow", part = "header")
#' 
#' # Example 12: Regulatory submission table
#' tbl2rtf(results, "submission.rtf",
#'         caption = "Table 2: Adjusted Odds Ratios for Mortality",
#'         font_family = "Times New Roman",
#'         font_size = 10,
#'         indent_groups = TRUE,
#'         zebra_stripes = FALSE,
#'         bold_significant = TRUE)
#' 
#' # Example 13: Custom column alignment
#' tbl2rtf(results, "aligned.rtf",
#'         align = c("left", "left", "center", "right", "right"))
#' 
#' # Example 14: Disable significance bolding
#' tbl2rtf(results, "no_bold.rtf",
#'         bold_significant = FALSE)
#' 
#' # Example 15: Stricter significance threshold
#' tbl2rtf(results, "strict.rtf",
#'         bold_significant = TRUE,
#'         p_threshold = 0.01)
#' 
#' # Example 16: Descriptive statistics for baseline characteristics
#' desc <- desctbl(
#'     data = clintrial,
#'     strata = "treatment",
#'     vars = c("age", "sex", "bmi", "stage"),
#'     labels = clintrial_labels
#' )
#' 
#' tbl2rtf(desc, "baseline.rtf",
#'         caption = "Table 1: Baseline Patient Characteristics",
#'         zebra_stripes = TRUE)
#' 
#' # Example 17: Clinical trial efficacy table
#' tbl2rtf(results, "efficacy.rtf",
#'         caption = "Table 3: Primary Efficacy Analysis - Intent to Treat Population",
#'         font_family = "Courier New",  # Monospace for alignment
#'         paper = "letter",
#'         orientation = "landscape",
#'         condense_table = TRUE)
#' }
#'
#' @export
tbl2rtf <- function(table,
                    file,
                    caption = NULL,
                    font_size = 8,
                    font_family = "Arial",
                    format_headers = TRUE,
                    bold_significant = TRUE,
                    p_threshold = 0.05,
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
    
    if (!grepl("\\.rtf$", tolower(file))) {
        stop("File must have .rtf extension")
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
        p_threshold = p_threshold,
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
    
                                        # Define page dimensions based on paper size
    page_width <- if(paper == "a4") 8.27 else 8.5
    page_height <- if(paper == "a4") 11.69 else if(paper == "legal") 14 else 11
    
    if (orientation == "landscape") {
        temp <- page_width
        page_width <- page_height
        page_height <- temp
    }
    
                                        # Set page size properties for the flextable
                                        # RTF uses page_size() from officer
    ps <- officer::page_size(width = page_width, 
                             height = page_height, 
                             orient = orientation)
    
    pm <- officer::page_mar(bottom = 1, top = 1, right = 1, left = 1)
    
                                        # Create section properties for RTF
    sect_properties <- officer::prop_section(
        page_size = ps,
        page_margins = pm,
        type = "continuous"
    )
    
                                        # Add section properties to flextable
                                        # This helps maintain consistent page layout in RTF
    ft <- flextable::set_table_properties(ft, layout = "autofit")
    
                                        # Save as RTF with caption if provided
    if (!is.null(result$caption)) {
                                        # Save flextable to temporary RTF to get complete, valid RTF structure
        temp_rtf <- tempfile(fileext = ".rtf")
        flextable::save_as_rtf(ft, path = temp_rtf)
        
                                        # Read the entire RTF content
        rtf_text <- paste(readLines(temp_rtf, warn = FALSE), collapse = "\n")
        
                                        # Find the first \trowd (table row definition)
                                        # We'll insert the caption right before it
        trowd_pos <- regexpr("\\\\trowd", rtf_text)
        
        if (trowd_pos[1] > 0) {
                                        # Split at the table start
            before_table <- substr(rtf_text, 1, trowd_pos[1] - 1)
            table_and_after <- substr(rtf_text, trowd_pos[1], nchar(rtf_text))
            
                                        # Create caption paragraph with proper RTF formatting
                                        # Use \qc for centered, or \ql for left-aligned
            caption_rtf <- sprintf(
                "\\pard\\ql\\f0\\fs%d %s\\par\n\\par\n",
                font_size * 2,  # RTF uses half-points
                result$caption
            )
            
                                        # Combine: before_table + caption + table_and_after
            rtf_final <- paste0(before_table, caption_rtf, table_and_after)
            
        } else {
                                        # Fallback if we can't find \trowd
                                        # Try to insert after the last closing brace of header tables
            warning("Could not find table start, attempting fallback caption insertion")
            
                                        # Find the position after {\colortbl...} and {\fonttbl...}
                                        # Look for the last } before \viewkind or similar
            insert_pattern <- "\\}\\s*\\\\viewkind"
            insert_pos <- regexpr(insert_pattern, rtf_text)
            
            if (insert_pos[1] > 0) {
                                        # Insert after the }
                pos <- insert_pos[1]
                caption_rtf <- sprintf(
                    "\n\\pard\\ql\\f0\\fs%d %s\\par\n\\par\n",
                    font_size * 2,
                    result$caption
                )
                rtf_final <- paste0(
                    substr(rtf_text, 1, pos),
                    caption_rtf,
                    substr(rtf_text, pos + 1, nchar(rtf_text))
                )
            } else {
                                        # Last resort: just prepend to the whole thing after header
                warning("Caption insertion failed, table will be generated without caption")
                rtf_final <- rtf_text
            }
        }
        
                                        # Write final RTF
        writeLines(rtf_final, file)
        
                                        # Clean up temp file
        unlink(temp_rtf)
        
                                        # Write final RTF
        writeLines(rtf_final, file)
        
                                        # Clean up temp file
        unlink(temp_rtf)
        
    } else {
                                        # No caption - direct RTF export
        flextable::save_as_rtf(ft, path = file)
    }
    
    message(sprintf("Table exported to %s", file))
    
                                        # Return based on user preference
    if (return_ft) {
                                        # Return flextable directly for easy access
        return(ft)
    } else {
                                        # Create a proper list to hold attributes
                                        # Return invisibly with flextable as attribute for backward compatibility
        result_obj <- list(
            file = file,
            caption = result$caption
        )
        attr(result_obj, "flextable") <- ft
        class(result_obj) <- c("tbl2rtf_result", "list")
        return(invisible(result_obj))
    }
}

#' Print method for tbl2rtf results
#' @keywords internal
#' @export
print.tbl2rtf_result <- function(x, ...) {
    cat("Table exported to:", x$file, "\n")
    if (!is.null(x$caption)) {
        cat("Caption:", x$caption, "\n")
    }
    cat("Flextable object available via: attr(result, 'flextable')\n")
    invisible(x)
}
