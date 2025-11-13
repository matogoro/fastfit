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
    
    ## Identify interaction terms (contain ":")
    is_interaction <- grepl(":", term_names, fixed = TRUE)
    interaction_terms <- term_names[is_interaction]
    main_terms <- term_names[!is_interaction]
    
    ## For models without factor variables, count unique base variables from main effects
    if (is.null(model$xlevels) || length(model$xlevels) == 0) {
        ## Count main effect terms
        n_terms <- length(main_terms)
        
        ## If we have interactions, extract unique base variables from them too
        if (length(interaction_terms) > 0) {
            ## Extract variables from interactions (split on ":")
            interaction_vars <- unique(unlist(strsplit(interaction_terms, ":", fixed = TRUE)))
            ## Count any interaction variables not already in main effects
            new_vars <- setdiff(interaction_vars, main_terms)
            n_terms <- n_terms + length(new_vars)
        }
        
        return(data.table::fifelse(n_terms == 1, "Univariable", "Multivariable"))
    }
    
    ## For models with factor variables, count unique base variables
    factor_pattern <- paste0("^(", paste(names(model$xlevels), collapse = "|"), ")")
    is_factor_term <- grepl(factor_pattern, main_terms)
    factor_vars_present <- unique(sapply(main_terms[is_factor_term], function(term) {
        for (var in names(model$xlevels)) {
            if (startsWith(term, var)) return(var)
        }
        return(NA_character_)
    }))
    factor_vars_present <- factor_vars_present[!is.na(factor_vars_present)]
    n_vars <- length(factor_vars_present)
    
    ## Count unique continuous variables IN MAIN EFFECTS
    continuous_terms <- main_terms[!is_factor_term]
    n_continuous <- length(unique(continuous_terms))
    
    ## Total unique variables (from main effects only)
    ## Interactions don't add new variables to the count
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
    
    ## Interactions should not be parsed as factor variables
    is_interaction <- grepl(":", terms, fixed = TRUE)
    
    if (!is.null(xlevels) && length(xlevels) > 0) {
        ## OPTIMIZED: Vectorized approach
        xlevel_names <- names(xlevels)
        
        ## For each factor variable, find all matching terms at once
        ## BUT SKIP INTERACTION TERMS
        for (var in xlevel_names) {
            ## Find which terms start with this variable name
            pattern <- paste0("^", var)
            matches <- grepl(pattern, terms) & !is_interaction  # Skip interactions
            
            if (any(matches)) {
                ## Extract levels for all matching terms at once
                variable[matches] <- var
                group[matches] <- sub(pattern, "", terms[matches])
            }
        }
        
        ## Any remaining unmatched terms (including interactions) are treated as-is
        unmatched <- variable == ""
        if (any(unmatched)) {
            variable[unmatched] <- terms[unmatched]
            ## group already initialized to "" for these
        }
    } else {
        ## No factor variables - all continuous (including interactions)
        variable <- terms
        ## group already initialized to ""
    }
    
    return(data.table::data.table(variable = variable, group = group))
}
