#' Set Variable Labels for Data Table Columns
#'
#' Efficiently sets label attributes for multiple columns in a data.table object.
#' This function follows data.table conventions and modifies the data.table by 
#' reference for optimal performance, making it ideal for labeling variables in
#' the fastfit package workflow.
#'
#' @param dt A data.table object to which labels will be applied.
#' 
#' @param ... Named arguments where names correspond to column names in \code{dt} 
#'   and values are the character strings to use as labels. This is the primary
#'   method for setting labels inline. All arguments must be named.
#'   
#' @param .labels Optional named list or named vector of labels where names 
#'   correspond to column names in \code{dt}. If provided, this takes precedence 
#'   over \code{...} arguments. This is useful when you have a pre-defined label
#'   dictionary or when programmatically setting labels. Default is \code{NULL}.
#'   
#' @param .warn Logical. If \code{TRUE}, issues a warning when any column names 
#'   in the labels are not found in \code{dt}. Useful for catching typos or 
#'   detecting when expected columns are missing. Default is \code{TRUE}.
#'
#' @return
#' Returns the modified data.table invisibly (allowing for potential piping or 
#' chaining). The primary purpose is the side effect of setting label attributes
#' by reference.
#'
#' @details
#' This function uses \code{data.table::setattr()} internally for optimal 
#' performance with data.table objects. Labels are stored as "label" attributes 
#' on individual columns and can be retrieved using \code{attr(dt$column, "label")}
#' or with \code{sapply(dt, function(x) attr(x, "label"))}.
#' 
#' The function modifies the data.table by reference (in-place) and does not 
#' create copies, making it memory efficient even for very large datasets. This
#' follows standard data.table conventions for reference semantics.
#' 
#' Only columns that exist in \code{dt} will have their labels set. If a column
#' name appears in the labels but not in \code{dt}, it will be skipped (with a
#' warning if \code{.warn = TRUE}).
#' 
#' Labels are particularly useful in the fastfit package for:
#' \itemize{
#'   \item Creating publication-ready tables with descriptive variable names
#'   \item Generating figures with proper axis labels
#'   \item Documenting datasets with metadata
#'   \item Maintaining consistency across analyses
#' }
#'
#' @seealso
#' \code{\link[data.table]{setattr}} for the underlying attribute-setting function,
#' \code{\link[data.table]{setnames}} for renaming columns,
#' \code{\link{attr}} for retrieving attributes
#'
#' @examples
#' # Load example data
#' data(clintrial)
#' library(data.table)
#' 
#' # Convert to data.table if needed
#' dt <- as.data.table(clintrial)
#' 
#' # Example 1: Set labels using named arguments (inline method)
#' setlabels(dt,
#'   age = "Age (years)",
#'   sex = "Sex",
#'   bmi = "Body Mass Index (kg/m²)",
#'   treatment = "Treatment Group"
#' )
#' 
#' # Example 2: Check that labels were set
#' attr(dt$age, "label")  # "Age (years)"
#' attr(dt$bmi, "label")  # "Body Mass Index (kg/m²)"
#' 
#' # Example 3: View all labels at once
#' sapply(dt, function(x) attr(x, "label"))
#' 
#' # Example 4: Using the .labels parameter with a named list
#' # (useful when labels are defined separately)
#' my_labels <- list(
#'   age = "Age in years",
#'   sex = "Biological sex", 
#'   bmi = "BMI (kg/m²)",
#'   os_months = "Overall Survival Time (months)",
#'   os_status = "Death Event"
#' )
#' setlabels(dt, .labels = my_labels)
#' 
#' # Example 5: Using the package's built-in label dictionary
#' data(clintrial_labels)
#' setlabels(dt, .labels = clintrial_labels)
#' 
#' # Example 6: Suppress warnings for missing columns
#' # (useful when applying a standard label set to different datasets)
#' setlabels(dt, 
#'   age = "Age (years)",
#'   nonexistent_col = "This column doesn't exist",
#'   .warn = FALSE  # No warning will be issued
#' )
#' 
#' # Example 7: Using with a named vector instead of list
#' label_vec <- c(
#'   "age" = "Patient Age",
#'   "sex" = "Patient Sex"
#' )
#' setlabels(dt, .labels = label_vec)
#' 
#' # Example 8: Chaining with other data.table operations
#' dt2 <- copy(dt)  # Make a copy
#' setlabels(dt2, 
#'   age = "Age at enrollment",
#'   bmi = "Body Mass Index"
#' )[age > 50]  # Can immediately use the data.table
#' 
#' # Example 9: Programmatically setting labels
#' # Build labels dynamically
#' clinical_vars <- c("age", "bmi", "creatinine")
#' label_list <- setNames(
#'   as.list(paste(clinical_vars, "(baseline)")),
#'   clinical_vars
#' )
#' setlabels(dt, .labels = label_list)
#' 
#' # Example 10: Labels persist through operations
#' dt_subset <- dt[age > 40]
#' attr(dt_subset$age, "label")  # Label is preserved
#' 
#' # Example 11: Update existing labels
#' # Setting a label again will overwrite the previous one
#' setlabels(dt, age = "Age (years, revised)")
#' attr(dt$age, "label")  # Now shows revised label
#' 
#' # Example 12: Using labels in table output functions
#' # (conceptual example - actual implementation depends on other fastfit functions)
#' setlabels(dt, .labels = clintrial_labels)
#' # Now descriptive() or other fastfit functions can use these labels
#' # instead of raw variable names
#'
#' @import data.table
#' @export
setlabels <- function(dt, ..., .labels = NULL, .warn = TRUE) {
    ## Input validation
    if (!is.data.table(dt)) {
        stop("dt must be a data.table object")
    }
    
    ## Handle labels input
    if (!is.null(.labels)) {
        
        ## Use .labels if provided
        if (is.vector(.labels) && !is.list(.labels)) {
            labels <- as.list(.labels)
        } else {
            labels <- .labels
        }
        
        if (is.null(names(labels)) || any(names(labels) == "")) {
            stop(".labels must be a named list or vector")
        }
    } else {
        
        ## Use ... arguments
        labels <- list(...)
        if (length(labels) == 0) {
            warning("No labels provided")
            return(invisible(dt))
        }
        
        if (is.null(names(labels)) || any(names(labels) == "")) {
            stop("All arguments must be named")
        }
    }
    
    ## Check for missing variables
    vars <- names(labels)
    missing_vars <- setdiff(vars, names(dt))
    
    if (length(missing_vars) > 0 && .warn) {
        warning("Variables not found in data.table: ", 
                paste(missing_vars, collapse = ", "))
    }
    
    ## Set labels efficiently
    existing_vars <- intersect(vars, names(dt))
    for (i in seq_along(existing_vars)) {
        var <- existing_vars[i]
        setattr(dt[[var]], "label", labels[[var]])
    }
    
    invisible(dt)
}
