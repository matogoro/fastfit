#' Detect if model is univariable or multivariable
#' @keywords internal
detect_model_type <- function(model) {
    ## Get number of non-intercept terms
    n_terms <- length(stats::coef(model))
    
    ## Account for intercept
    if ("(Intercept)" %in% names(stats::coef(model))) {
        n_terms <- n_terms - 1
    }
    
    ## For Cox models, no intercept to worry about
    if (inherits(model, c("coxph", "coxme", "clogit"))) {
        n_terms <- length(stats::coef(model))
    }
    
    ## Check for factor expansions - if model has xlevels, count base variables
    if (!is.null(model$xlevels)) {
        n_vars <- length(model$xlevels)
        ## Add any continuous variables (those not in xlevels)
        term_names <- names(stats::coef(model))
        term_names <- term_names[term_names != "(Intercept)"]
        for (term in term_names) {
            is_factor_term <- FALSE
            for (var in names(model$xlevels)) {
                if (grepl(paste0("^", var), term)) {
                    is_factor_term <- TRUE
                    break
                }
            }
            if (!is_factor_term) n_vars <- n_vars + 1
        }
        n_terms <- n_vars
    }
    
    return(ifelse(n_terms == 1, "Univariable", "Multivariable"))
}

#' Get readable model type name
#' @keywords internal
get_model_type_name <- function(model) {
    model_class <- class(model)[1]
    
    ## Remove wrapper classes
    if (model_class == "mmodel") {
        model_class <- class(model)[2]
    }
    
    ## Map to readable names
    type_map <- c(
        "glm" = "Logistic",  # Will be refined based on family
        "lm" = "Linear",
        "coxph" = "Cox PH",
        "clogit" = "Conditional Logistic",
        "coxme" = "Mixed Effects Cox",
        "glmer" = "Mixed Effects GLM",
        "lmer" = "Mixed Effects Linear"
    )
    
    ## For GLM, be more specific based on family
    if (model_class == "glm") {
        family <- model$family$family
        link <- model$family$link
        
        if (family == "binomial") {
            return("Logistic")
        } else if (family == "poisson") {
            return("Poisson")
        } else if (family == "gaussian") {
            return("Linear (GLM)")
        } else if (family == "Gamma") {
            return("Gamma")
        } else if (family == "quasibinomial") {
            return("Quasi-Binomial")
        } else if (family == "quasipoisson") {
            return("Quasi-Poisson")
        } else {
            return(paste0(stringr::str_to_title(family), " GLM"))
        }
    }
    
    return(type_map[model_class] %||% model_class)
}

#' Parse term into variable and group
#' @keywords internal
parse_term <- function(terms, xlevels = NULL) {
    result <- data.table::data.table(
                              variable = character(length(terms)),
                              group = character(length(terms))
                          )
    
    for (i in seq_along(terms)) {
        term <- terms[i]
        
        ## Check if it's a factor term
        matched <- FALSE
        if (!is.null(xlevels)) {
            for (var in names(xlevels)) {
                if (grepl(paste0("^", var), term)) {
                    ## Extract the level
                    level <- sub(paste0("^", var), "", term)
                    result$variable[i] <- var  ## lowercase
                    result$group[i] <- level    ## lowercase
                    matched <- TRUE
                    break
                }
            }
        }
        
        ## If not matched as factor, it's continuous or interaction
        if (!matched) {
            result$variable[i] <- term
            result$group[i] <- ""
        }
    }
    
    return(result)
}
