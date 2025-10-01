#' Convert Model Object to Data Table
#'
#' Extracts coefficients, confidence intervals, and comprehensive model statistics 
#' from fitted regression models and converts them to a standardized data.table 
#' format suitable for further analysis or publication.
#'
#' @param model A fitted model object. Supported classes include glm (generalized
#'   linear models), lm (linear models), coxph (Cox proportional hazards),
#'   clogit (conditional logistic), coxme (mixed effects Cox), and glmer
#'   (generalized linear mixed effects). Also accepts models wrapped with mmodel().
#' @param conf_level Numeric confidence level for confidence intervals. Must be
#'   between 0 and 1. Default is 0.95 (95\% CI).
#' @param variable_name Character string to identify the model source in the output.
#'   Useful when combining results from multiple models. Default is NULL, which
#'   uses "Multivariable" as the identifier.
#' @param keep_qc_stats Logical. If TRUE, includes model quality statistics such as
#'   AIC, BIC, R-squared, concordance, and model fit tests. These appear as 
#'   additional columns in the output. Default is TRUE.
#' @param terms_to_exclude Character vector of term names to exclude from output.
#'   Useful for removing intercepts or other unwanted parameters. Default is
#'   "(Intercept)".
#' @param add_reference_rows Logical. If TRUE, adds rows for reference categories
#'   of factor variables with appropriate labels and baseline values (OR/HR = 1,
#'   Estimate = 0). Default is TRUE.
#' @param reference_label Character string used to label reference category rows
#'   in the output. Default is "reference".
#'
#' @return A data.table containing extracted model information with the following
#'   standard columns:
#'   - variable: Model identifier from variable_name parameter
#'   - term: Predictor/coefficient name
#'   - n: Sample size
#'   - events: Number of events (for survival/logistic models)
#'   - coefficient: Raw coefficient estimate
#'   - se: Standard error
#'   - OR/HR/RR/Estimate: Effect estimate (type depends on model)
#'   - CI_lower, CI_upper: Confidence interval bounds
#'   - statistic: Test statistic (z or t value)
#'   - p_value: P-value for coefficient test
#'   - sig: Significance markers (*** p<0.001, ** p<0.01, * p<0.05)
#'   - sig_binary: Logical indicator for p<0.05
#'   
#'   Additional quality control columns when keep_qc_stats = TRUE vary by model type:
#'   - GLM: AIC, BIC, deviance, null_deviance, c_statistic (logistic)
#'   - LM: R2, adj_R2, sigma, df_residual
#'   - Cox: concordance, rsq, logtest_stat, wald_test, score_test
#'   
#'   The output includes attributes:
#'   - model_class: The class of the input model
#'   - model_family: The family for GLM models
#'   - conf_level: The confidence level used
#'
#' @export
m2dt <- function(model, 
                 conf_level = 0.95,
                 variable_name = NULL,
                 keep_qc_stats = TRUE,
                 terms_to_exclude = "(Intercept)",
                 add_reference_rows = TRUE,
                 reference_label = "reference") {
    
    model_class <- class(model)[1]
    
    ## Remove mmodel class if present to get underlying model type
    if (model_class == "mmodel") {
        model_class <- class(model)[2]
    }
    
    ## Initialize base data.table
    dt <- data.table::data.table()
    
    ## Extract based on model type
    if (model_class %in% c("glm", "lm")) {
        
        coef_summary <- stats::coef(summary(model))
        z_score <- stats::qnorm((1 + conf_level) / 2)
        
        ## Determine effect measure
        is_logistic <- model_class == "glm" && model$family$family == "binomial"
        is_poisson <- model_class == "glm" && model$family$family == "poisson"
        should_exp <- is_logistic || is_poisson || 
            (model_class == "glm" && model$family$link == "log")
        
        effect_name <- if (is_logistic) "OR" else if (is_poisson) "RR" else "Estimate"
        
        dt <- data.table::data.table(
                              variable = variable_name %||% "Multivariable",
                              term = rownames(coef_summary),
                              n = stats::nobs(model),
                              events = if (is_logistic && !is.null(model$y)) sum(model$y) else NA_integer_,
                              coefficient = coef_summary[, "Estimate"],
                              se = coef_summary[, "Std. Error"]
                          )
        
        ## Add effect estimates
        if (should_exp) {
            dt[, `:=`(
                effect = exp(coefficient),
                CI_lower = exp(coefficient - z_score * se),
                CI_upper = exp(coefficient + z_score * se)
            )]
        } else {
            dt[, `:=`(
                effect = coefficient,
                CI_lower = coefficient - z_score * se,
                CI_upper = coefficient + z_score * se
            )]
        }
        
        ## Rename effect column
        data.table::setnames(dt, "effect", effect_name)
        
        ## Add test statistics
        stat_col <- if ("z value" %in% colnames(coef_summary)) "z value" else "t value"
        dt[, `:=`(
            statistic = coef_summary[, stat_col],
            p_value = coef_summary[, ncol(coef_summary)]
        )]
        
        ## Add QC stats if requested
        if (keep_qc_stats) {
            if (model_class == "glm") {
                dt[, `:=`(
                    AIC = stats::AIC(model),
                    BIC = stats::BIC(model),
                    deviance = stats::deviance(model),
                    null_deviance = model$null.deviance,
                    df_residual = stats::df.residual(model)
                )]
                
                ## Add R-squared for non-binomial GLMs
                if (model$family$family != "binomial") {
                    dt[, R2 := 1 - (deviance / null_deviance)]
                }
                
                ## For binomial, add discrimination/calibration metrics
                if (is_logistic && keep_qc_stats) {
                    ## C-statistic (if pROC available)
                    if (requireNamespace("pROC", quietly = TRUE)) {
                        roc_obj <- pROC::roc(model$y, stats::fitted(model), quiet = TRUE)
                        dt[, c_statistic := as.numeric(pROC::auc(roc_obj))]
                    }
                    
                    ## Hosmer-Lemeshow test (if ResourceSelection available)
                    if (requireNamespace("ResourceSelection", quietly = TRUE)) {
                        hl <- ResourceSelection::hoslem.test(model$y, stats::fitted(model), g = 10)
                        dt[, `:=`(
                            hoslem_chi2 = hl$statistic,
                            hoslem_p = hl$p.value
                        )]
                    }
                }
            } else if (model_class == "lm") {
                summ <- summary(model)
                dt[, `:=`(
                    R2 = summ$r.squared,
                    adj_R2 = summ$adj.r.squared,
                    AIC = stats::AIC(model),
                    BIC = stats::BIC(model),
                    sigma = summ$sigma,
                    df_residual = stats::df.residual(model)
                )]
            }
        }
        
    } else if (model_class %in% c("coxph", "clogit")) {
        
        if (!requireNamespace("survival", quietly = TRUE))
            stop("Package 'survival' required")
        
        summ <- summary(model)
        coef_summary <- summ$coefficients
        conf_int <- stats::confint(model, level = conf_level)
        
        dt <- data.table::data.table(
                              variable = variable_name %||% "Multivariable",
                              term = rownames(coef_summary),
                              n = if (!is.null(model$n)) model$n[1] else summ$n,
                              events = if (!is.null(model$nevent)) model$nevent 
                                       else if (!is.null(model$n)) model$n[2] 
                                       else summ$nevent,
                              coefficient = coef_summary[, "coef"],
                              se = coef_summary[, "se(coef)"],
                              HR = coef_summary[, "exp(coef)"],
                              CI_lower = exp(conf_int[, 1]),
                              CI_upper = exp(conf_int[, 2]),
                              statistic = coef_summary[, "z"],
                              p_value = coef_summary[, "Pr(>|z|)"]
                          )
        
        ## Add QC stats
        if (keep_qc_stats) {
            dt[, `:=`(
                concordance = summ$concordance[1],
                concordance_se = summ$concordance[2],
                rsq = summ$rsq[1],
                rsq_max = summ$rsq[2],
                logtest_stat = summ$logtest[1],
                logtest_p = summ$logtest[3],
                wald_test = summ$waldtest[1],
                wald_p = summ$waldtest[3],
                score_test = summ$sctest[1],
                score_p = summ$sctest[3]
            )]
            
            ## Add AIC/BIC
            dt[, `:=`(
                AIC = stats::AIC(model),
                BIC = stats::extractAIC(model, k = log(n))[2]
            )]
        }
        
    } else if (model_class == "coxme") {
        
        if (!requireNamespace("coxme", quietly = TRUE))
            stop("Package 'coxme' required")
        
        coef_vals <- stats::coef(model)
        se_vals <- sqrt(diag(stats::vcov(model)))
        z_vals <- coef_vals / se_vals
        p_vals <- 2 * stats::pnorm(-abs(z_vals))
        z_score <- stats::qnorm((1 + conf_level) / 2)
        
        dt <- data.table::data.table(
                              variable = variable_name %||% "Multivariable", 
                              term = names(coef_vals),
                              n = model$n[1],
                              events = model$n[2],
                              coefficient = coef_vals,
                              se = se_vals,
                              HR = exp(coef_vals),
                              CI_lower = exp(coef_vals - z_score * se_vals),
                              CI_upper = exp(coef_vals + z_score * se_vals),
                              statistic = z_vals,
                              p_value = p_vals
                          )
        
        if (keep_qc_stats) {
            dt[, `:=`(
                loglik = stats::logLik(model),
                df = attr(stats::logLik(model), "df"),
                AIC = stats::AIC(model)
            )]
        }
        
    } else {
        stop("Model class '", model_class, "' not yet supported")
    }
    
    ## Remove excluded terms
    if (!is.null(terms_to_exclude)) {
        dt <- dt[!term %in% terms_to_exclude]
    }
    
    ## Add reference rows for factor variables
    if (add_reference_rows && !is.null(model$xlevels)) {
        
        ## Add reference column to existing data
        dt[, reference := ""]
        
        ## Create a list to hold the final ordered result
        final_dt <- list()
        
        ## Track which terms have been processed
        processed_terms <- character(0)
        
        ## Process each term in the original order
        for (i in seq_len(nrow(dt))) {
            current_term <- dt$term[i]
            
            ## Skip if already processed
            if (current_term %in% processed_terms) next
            
            ## Check if this term belongs to a factor variable
            factor_var <- NULL
            for (var in names(model$xlevels)) {
                if (grepl(paste0("^", var), current_term)) {
                    factor_var <- var
                    break
                }
            }
            
            if (!is.null(factor_var)) {
                ## This is a factor variable
                levels_order <- model$xlevels[[factor_var]]
                ref_level <- levels_order[1]
                
                ## Create reference row
                ref_row <- dt[1, ]
                ref_row[, `:=`(
                    term = paste0(factor_var, ref_level),
                    coefficient = 0,
                    se = NA_real_,
                    statistic = NA_real_,
                    p_value = NA_real_,
                    CI_lower = NA_real_,
                    CI_upper = NA_real_,
                    reference = reference_label
                )]
                
                ## Set effect estimates for reference
                if ("OR" %in% names(dt)) ref_row[, OR := 1]
                if ("HR" %in% names(dt)) ref_row[, HR := 1]
                if ("RR" %in% names(dt)) ref_row[, RR := 1]
                if ("Estimate" %in% names(dt)) ref_row[, Estimate := 0]
                
                ## Add reference row first
                final_dt[[length(final_dt) + 1]] <- ref_row
                
                ## Now add all other levels for this factor in order
                for (level in levels_order[-1]) {
                    level_term <- paste0(factor_var, level)
                    level_row <- dt[term == level_term]
                    if (nrow(level_row) > 0) {
                        final_dt[[length(final_dt) + 1]] <- level_row
                        processed_terms <- c(processed_terms, level_term)
                    }
                }
                
            } else {
                ## This is not a factor variable (e.g., continuous)
                final_dt[[length(final_dt) + 1]] <- dt[term == current_term]
                processed_terms <- c(processed_terms, current_term)
            }
        }
        
        ## Combine all rows
        dt <- rbindlist(final_dt, fill = TRUE)
    }

    ## Add significance markers
    dt[, `:=`(
        sig = data.table::fcase(
                              is.na(p_value), "",
                              p_value < 0.001, "***",
                              p_value < 0.01, "**",
                              p_value < 0.05, "*",
                              p_value < 0.10, ".",
                              default = ""
                          ),
        sig_binary = !is.na(p_value) & p_value < 0.05
    )]
    
    ## Set attributes for model info
    data.table::setattr(dt, "model_class", model_class)
    data.table::setattr(dt, "model_family", if (model_class == "glm") model$family$family else NA)
    data.table::setattr(dt, "conf_level", conf_level)

    ## Force data.table to finalize
    dt[]  # This forces data.table to update its internal state
    
    return(dt)
}
