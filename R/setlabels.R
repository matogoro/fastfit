#' Set Variable Labels for data.table Columns
#'
#' @description 
#' Efficiently set label attributes for multiple columns in a data.table object.
#' This function follows data.table conventions and modifies the data.table by 
#' reference for optimal performance.
#'
#' @param dt A data.table object
#' @param ... Named arguments where names correspond to column names in \code{dt} 
#'   and values are the labels to assign
#' @param .labels Alternative way to pass labels as a named list or vector 
#'   (optional, overrides \code{...} if provided)
#' @param .warn Logical. Whether to warn about columns not found in \code{dt}. 
#'   Default is \code{TRUE}
#'
#' @details
#' This function uses \code{setattr()} internally for optimal performance with 
#' data.table objects. Labels are stored as "label" attributes and can be 
#' retrieved using \code{attr(dt$column, "label")}.
#' 
#' The function modifies the data.table by reference and does not create copies,
#' making it memory efficient for large datasets.
#'
#' @return
#' Returns the modified data.table invisibly (for potential chaining), though
#' the primary purpose is the side effect of setting labels by reference.
#'
#' @examples
#' library(data.table)
#' 
#' # Create sample data
#' dt <- data.table(
#'   age = c(25, 30, 35),
#'   income = c(50000, 75000, 90000),
#'   education = c("HS", "College", "Graduate")
#' )
#' 
#' # Set labels using named arguments
#' setlabels(dt,
#'   age = "Age in years",
#'   income = "Annual income in USD",
#'   education = "Education level"
#' )
#' 
#' # Check labels
#' attr(dt$age, "label")
#' sapply(dt, function(x) attr(x, "label"))
#' 
#' # Using .labels parameter
#' label_dict <- list(
#'   age = "Age in years", 
#'   income = "Annual income in USD"
#' )
#' setlabels(dt, .labels = label_dict)
#' 
#' # Suppress warnings for missing columns
#' setlabels(dt, 
#'   age = "Age in years",
#'   nonexistent_col = "This will warn",
#'   .warn = FALSE
#' )
#'
#' @seealso
#' \code{\link[data.table]{setattr}}, \code{\link[data.table]{setnames}}, 
#' \code{\link{attr}}
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
