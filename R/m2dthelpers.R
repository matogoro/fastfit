#' Detect if model is univariable or multivariable
#' @keywords internal
detect_model_type <- function(model) {
    ## Get coefficient names once
    coef_names <- names(stats::coef(model))
    
    ## Remove intercept from count
    if ("(Intercept)" %in% coef_names) {
        term_names <- coef_names[coef_names != "(Intercept)"]
    } else {
        term_names <- coef_names
    }
    
    ## For models without factor variables, just count terms
    if (is.null(model$xlevels) || length(model$xlevels) == 0) {
        n_terms <- length(term_names)
        return(data.table::fifelse(n_terms == 1, "Univariable", "Multivariable"))
    }
    
    ## For models with factor variables, count unique base variables
    ## Start with number of factor variables
    n_vars <- length(model$xlevels)
    
    ## Create regex pattern to match all factor variable terms
    factor_pattern <- paste0("^(", paste(names(model$xlevels), collapse = "|"), ")")
    
    ## Find which terms are NOT factor terms (i.e., continuous variables)
    is_factor_term <- grepl(factor_pattern, term_names)
    
    ## Count unique continuous variables
    ## Each continuous term that's not a factor level is a separate variable
    continuous_terms <- term_names[!is_factor_term]
    n_continuous <- length(unique(continuous_terms))
    
    ## Total unique variables
    n_terms <- n_vars + n_continuous
    
    return(data.table::fifelse(n_terms == 1, "Univariable", "Multivariable"))
}

#' Get readable model type name
#' @keywords internal
get_model_type_name <- function(model) {
    model_class <- class(model)[1]
    
    ## Remove wrapper classes
    if (model_class == "mmodel") {
        model_class <- class(model)[2]
    }
    
    ## For GLM, be more specific based on family
    if (model_class == "glm") {
        family <- model$family$family
        
        ## OPTIMIZED: Use switch instead of multiple if-else
        return(switch(family,
                      binomial = "Logistic",
                      poisson = "Poisson",
                      gaussian = "Linear (GLM)",
                      Gamma = "Gamma",
                      quasibinomial = "Quasi-Binomial",
                      quasipoisson = "Quasi-Poisson",
                      paste0(stringr::str_to_title(family), " GLM")
                      ))
    }
    
    ## Map to readable names for other model types
    ## OPTIMIZED: Use switch instead of named vector lookup
    return(switch(model_class,
                  lm = "Linear",
                  coxph = "Cox PH",
                  clogit = "Conditional Logistic",
                  coxme = "Mixed Effects Cox",
                  glmer = "Mixed Effects GLM",
                  lmer = "Mixed Effects Linear",
                  model_class  # default
                  ))
}

#' Parse term into variable and group
#' @keywords internal
parse_term <- function(terms, xlevels = NULL) {
    n_terms <- length(terms)
    
    ## Initialize result vectors
    variable <- character(n_terms)
    group <- character(n_terms)
    
    if (!is.null(xlevels) && length(xlevels) > 0) {
        ## OPTIMIZED: Vectorized approach
        xlevel_names <- names(xlevels)
        
        ## For each factor variable, find all matching terms at once
        for (var in xlevel_names) {
            ## Find which terms start with this variable name
            pattern <- paste0("^", var)
            matches <- grepl(pattern, terms)
            
            if (any(matches)) {
                ## Extract levels for all matching terms at once
                variable[matches] <- var
                group[matches] <- sub(pattern, "", terms[matches])
            }
        }
        
        ## Any remaining unmatched terms are continuous variables
        unmatched <- variable == ""
        if (any(unmatched)) {
            variable[unmatched] <- terms[unmatched]
            ## group already initialized to "" for these
        }
    } else {
        ## No factor variables - all continuous
        variable <- terms
        ## group already initialized to ""
    }
    
    return(data.table::data.table(variable = variable, group = group))
}
