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
#' @param keep_qc_stats Logical. If TRUE, includes model quality statistics such as
#'   AIC, BIC, R-squared, concordance, and model fit tests. These appear as 
#'   additional columns in the output. Default is TRUE.
#' @param include_intercept Logical. If TRUE, includes the model intercept in output.
#'   If FALSE, removes the intercept row from results. Useful for creating cleaner
#'   presentation tables. Default is TRUE for backward compatibility.
#' @param terms_to_exclude Character vector of term names to exclude from output.
#'   Useful for removing specific unwanted parameters. Default is NULL. Note: If
#'   include_intercept is FALSE, "(Intercept)" is automatically added to this list.
#' @param add_reference_rows Logical. If TRUE, adds rows for reference categories
#'   of factor variables with appropriate labels and baseline values (OR/HR = 1,
#'   Estimate = 0). Default is TRUE.
#' @param reference_label Character string used to label reference category rows
#'   in the output. Default is "reference".
#'
#' @return A data.table containing extracted model information with the following
#'   standard columns:
#'   - model_scope: Univariable (unadjusted) or multivariable (adjusted) model
#'   - model_type: Type of regression model
#'   - variable: Variable summarized
#'   - group: Group summarized
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
#' @keywords internal
#' @export
m2dt <- function(model, 
                 conf_level = 0.95,
                 keep_qc_stats = TRUE,
                 include_intercept = TRUE,
                 terms_to_exclude = NULL,
                 add_reference_rows = TRUE,
                 reference_label = "reference") {
    
                                        # ENHANCED: Handle intercept exclusion
                                        # If include_intercept is FALSE, add "(Intercept)" to exclusion list
    if (!include_intercept) {
        if (is.null(terms_to_exclude)) {
            terms_to_exclude <- "(Intercept)"
        } else {
            terms_to_exclude <- unique(c(terms_to_exclude, "(Intercept)"))
        }
    }
    
                                        # For backward compatibility: if terms_to_exclude is still NULL at this point,
                                        # set it to empty character vector (no exclusions)
    if (is.null(terms_to_exclude)) {
        terms_to_exclude <- character(0)
    }
    
    model_class <- class(model)[1]
    
    ## Remove mmodel class if present to get underlying model type
    if (model_class == "mmodel") {
        model_class <- class(model)[2]
    }

    ## Auto-detect model type if not specified
    model_scope <- detect_model_type(model)

    ## Get model type name
    model_type_name <- get_model_type_name(model)
    
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

        ## Events calculation for logistic regression
        events_value <- NA_integer_
        if (is_logistic && !is.null(model$y)) {
            if (is.factor(model$y)) {
                ## Convert factor to numeric (assumes binary factor with levels 0/1 or similar)
                events_value <- sum(as.numeric(model$y) == 2)  # Second level is typically "Yes" or "1"
            } else {
                events_value <- sum(model$y)
            }
        }
        
        dt <- data.table::data.table(
                              model_scope = model_scope %||% "Multivariable",
                              model_type = model_type_name,
                              term = rownames(coef_summary),
                              n = stats::nobs(model),
                              events = events_value,
                              coefficient = coef_summary[, "Estimate"],
                              se = coef_summary[, "Std. Error"]
                          )

        dt[, `:=`(
            ## Raw coefficients
            coef = coefficient,
            coef_lower = coefficient - z_score * se,
            coef_upper = coefficient + z_score * se,
            
            ## Exponentiated versions
            exp_coef = exp(coefficient),
            exp_lower = exp(coefficient - z_score * se),
            exp_upper = exp(coefficient + z_score * se)
        )]
        
        if (is_logistic) {
            dt[, `:=`(
                OR = exp_coef,
                CI_lower = exp_lower,
                CI_upper = exp_upper
            )]
        } else if (is_poisson) {
            dt[, `:=`(
                RR = exp_coef,
                CI_lower = exp_lower,
                CI_upper = exp_upper
            )]
        } else if (should_exp) {
            ## Other log-link models
            dt[, `:=`(
                Estimate = exp_coef,
                CI_lower = exp_lower,
                CI_upper = exp_upper
            )]
        } else {
            ## Linear models - use raw coefficients
            dt[, `:=`(
                Estimate = coef,
                CI_lower = coef_lower,
                CI_upper = coef_upper
            )]
        }
        
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
                              model_scope = model_scope %||% "Multivariable",
                              model_type = model_type_name,
                              term = rownames(coef_summary),
                              n = if (!is.null(model$n)) model$n[1] else summ$n,
                              events = if (!is.null(model$nevent)) model$nevent 
                                       else if (!is.null(model$n)) model$n[2] 
                                       else summ$nevent,
                              coefficient = coef_summary[, "coef"],
                              se = coef_summary[, "se(coef)"],
                                        # Store both versions
                              coef = coef_summary[, "coef"],
                              coef_lower = conf_int[, 1],
                              coef_upper = conf_int[, 2],
                              exp_coef = coef_summary[, "exp(coef)"],
                              exp_lower = exp(conf_int[, 1]),
                              exp_upper = exp(conf_int[, 2]),
                                        # Primary display columns
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
                likelihood_ratio_test = summ$logtest[1],
                likelihood_ratio_df = summ$logtest[2],
                likelihood_ratio_p = summ$logtest[3],
                wald_test = summ$waldtest[1],
                wald_df = summ$waldtest[2],
                wald_p = summ$waldtest[3],
                score_test = summ$sctest[1],
                score_df = summ$sctest[2],
                score_p = summ$sctest[3]
            )]
        }
        
    } else if (model_class %in% c("coxme", "lme", "lmer", "glmer")) {
        
        ## Mixed effects models
        summ <- summary(model)
        
        if (model_class == "coxme") {
            coef_vec <- stats::fixef(model)
            vcov_mat <- as.matrix(stats::vcov(model))
            se_vec <- sqrt(diag(vcov_mat))
            z_vec <- coef_vec / se_vec
            p_vec <- 2 * (1 - stats::pnorm(abs(z_vec)))
            
            dt <- data.table::data.table(
                                  model_scope = model_scope %||% "Multivariable",
                                  model_type = model_type_name,
                                  term = names(coef_vec),
                                  n = model$n[1],
                                  events = model$n[2],
                                  coefficient = coef_vec,
                                  se = se_vec,
                                  coef = coef_vec,
                                  coef_lower = coef_vec - stats::qnorm((1 + conf_level) / 2) * se_vec,
                                  coef_upper = coef_vec + stats::qnorm((1 + conf_level) / 2) * se_vec,
                                  exp_coef = exp(coef_vec),
                                  exp_lower = exp(coef_vec - stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  exp_upper = exp(coef_vec + stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  HR = exp(coef_vec),
                                  CI_lower = exp(coef_vec - stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  CI_upper = exp(coef_vec + stats::qnorm((1 + conf_level) / 2) * se_vec),
                                  statistic = z_vec,
                                  p_value = p_vec
                              )
            
        } else {
            ## lmer/glmer from lme4
            if (!requireNamespace("lme4", quietly = TRUE))
                stop("Package 'lme4' required")
            
            coef_summary <- stats::coef(summ)
            
            ## Determine if should exponentiate
            should_exp <- model_class == "glmer" && 
                (summ$family == "binomial" || summ$link == "log")
            
            dt <- data.table::data.table(
                                  model_scope = model_scope %||% "Multivariable",
                                  model_type = model_type_name,
                                  term = rownames(coef_summary),
                                  n = stats::nobs(model),
                                  events = NA_integer_,
                                  coefficient = coef_summary[, "Estimate"],
                                  se = coef_summary[, "Std. Error"]
                              )
            
            z_score <- stats::qnorm((1 + conf_level) / 2)
            dt[, `:=`(
                coef = coefficient,
                coef_lower = coefficient - z_score * se,
                coef_upper = coefficient + z_score * se,
                exp_coef = if (should_exp) exp(coefficient) else coefficient,
                exp_lower = if (should_exp) exp(coefficient - z_score * se) else coefficient - z_score * se,
                exp_upper = if (should_exp) exp(coefficient + z_score * se) else coefficient + z_score * se
            )]
            
            ## Add appropriate effect column
            if (model_class == "glmer" && summ$family == "binomial") {
                dt[, `:=`(
                    OR = exp_coef,
                    CI_lower = exp_lower,
                    CI_upper = exp_upper
                )]
            } else if (should_exp) {
                dt[, `:=`(
                    RR = exp_coef,
                    CI_lower = exp_lower,
                    CI_upper = exp_upper
                )]
            } else {
                dt[, `:=`(
                    Estimate = coef,
                    CI_lower = coef_lower,
                    CI_upper = coef_upper
                )]
            }
            
            ## Add test statistics
            if ("t value" %in% colnames(coef_summary)) {
                dt[, `:=`(
                    statistic = coef_summary[, "t value"],
                    p_value = if ("Pr(>|t|)" %in% colnames(coef_summary)) 
                                  coef_summary[, "Pr(>|t|)"] else NA_real_
                )]
            } else if ("z value" %in% colnames(coef_summary)) {
                dt[, `:=`(
                    statistic = coef_summary[, "z value"],
                    p_value = if ("Pr(>|z|)" %in% colnames(coef_summary)) 
                                  coef_summary[, "Pr(>|z|)"] else NA_real_
                )]
            }
        }
        
    } else {
        stop("Unsupported model class: ", model_class)
    }
    
                                        # ENHANCED: Filter out excluded terms before any further processing
                                        # This ensures intercepts are removed early in the pipeline
    if (length(terms_to_exclude) > 0) {
        dt <- dt[!term %in% terms_to_exclude]
    }

    ## Process terms to extract variable and group information
    if (!("variable" %in% names(dt))) {
        dt[, `:=`(variable = "", group = "")]
        
        for (i in seq_len(nrow(dt))) {
            term_str <- dt$term[i]
            
                                        # NOTE: We no longer need to check terms_to_exclude here
                                        # because we've already filtered them out above
            
            ## Check against xlevels if available
            if (!is.null(model$xlevels)) {
                matched <- FALSE
                for (var in names(model$xlevels)) {
                    if (grepl(paste0("^", var), term_str)) {
                        level <- gsub(paste0("^", var), "", term_str)
                        dt[i, `:=`(variable = var, group = level)]
                        matched <- TRUE
                        break
                    }
                }
                if (!matched) {
                    dt[i, `:=`(variable = term_str, group = "")]
                }
            } else {
                dt[i, `:=`(variable = term_str, group = "")]
            }
        }
        
        ## Add n_group and events_group for categorical variables
        if (!is.null(model$xlevels)) {
            
            ## Try to get the source data
            data_source <- NULL
            if (!is.null(model$data)) {
                data_source <- model$data
            } else if (!is.null(model$model)) {
                data_source <- model$model
            }
            
            for (var in names(model$xlevels)) {
                levels <- model$xlevels[[var]]
                
                ## Calculate n for each level if data available
                if (!is.null(data_source) && var %in% names(data_source)) {
                    for (level in levels) {
                        level_n <- sum(data_source[[var]] == level, na.rm = TRUE)
                        dt[variable == var & group == level, n_group := level_n]
                        
                        ## For binomial GLM, also calculate events
                        if (model_class == "glm" && model$family$family == "binomial") {
                            outcome_var <- all.vars(model$formula)[1]
                            if (outcome_var %in% names(data_source)) {
                                level_data <- data_source[data_source[[var]] == level & !is.na(data_source[[var]]), ]
                                        # FIX: Handle factor outcomes
                                if (is.factor(level_data[[outcome_var]])) {
                                    level_events <- sum(as.numeric(level_data[[outcome_var]]) == 2, na.rm = TRUE)
                                } else {
                                    level_events <- sum(level_data[[outcome_var]], na.rm = TRUE)
                                }
                                dt[variable == var & group == level, events_group := level_events]
                            }
                        } else if (model_class %in% c("coxph", "clogit")) {
                            ## For survival, extract from Surv() formula
                            outcome_str <- as.character(model$formula)[2]
                            if (grepl("^Surv\\(", outcome_str)) {
                                surv_expr <- gsub("Surv\\(|\\)", "", outcome_str)
                                surv_parts <- trimws(strsplit(surv_expr, ",")[[1]])
                                event_var <- surv_parts[2]
                                
                                if (event_var %in% names(data_source)) {
                                    level_data <- data_source[data_source[[var]] == level & !is.na(data_source[[var]]), ]
                                    level_events <- sum(level_data[[event_var]], na.rm = TRUE)
                                    dt[variable == var & group == level, events_group := level_events]
                                }
                            }
                        }
                    }
                }
            }
        }
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
                
                ## Calculate n and events for reference level
                ref_n_group <- NA_real_
                ref_events_group <- NA_real_
                
                ## Get the data source (same as used for other levels)
                data_source <- NULL
                if (!is.null(model$data)) {
                    data_source <- model$data
                } else if (!is.null(model$model) && factor_var %in% names(model$model)) {
                    data_source <- model$model
                }
                
                if (!is.null(data_source) && factor_var %in% names(data_source)) {
                    ## Calculate reference level counts
                    ref_n_group <- sum(data_source[[factor_var]] == ref_level, na.rm = TRUE)
                    
                    if (model_class == "glm" && model$family$family == "binomial") {
                        ## Get outcome variable and calculate events
                        outcome_var <- all.vars(model$formula)[1]
                        if (outcome_var %in% names(data_source)) {
                            ref_data <- data_source[data_source[[factor_var]] == ref_level & !is.na(data_source[[factor_var]]), ]
                                        # FIX: Handle factor outcomes
                            if (is.factor(ref_data[[outcome_var]])) {
                                ref_events_group <- sum(as.numeric(ref_data[[outcome_var]]) == 2, na.rm = TRUE)
                            } else {
                                ref_events_group <- sum(ref_data[[outcome_var]], na.rm = TRUE)
                            }
                        }
                    } else if (model_class %in% c("coxph", "clogit")) {
                        outcome_str <- as.character(model$formula)[2]
                        if (grepl("^Surv\\(", outcome_str)) {
                            surv_expr <- gsub("Surv\\(|\\)", "", outcome_str)
                            surv_parts <- trimws(strsplit(surv_expr, ",")[[1]])
                            event_var <- surv_parts[2]
                            
                            if (event_var %in% names(data_source)) {
                                ref_data <- data_source[data_source[[factor_var]] == ref_level & !is.na(data_source[[factor_var]]), ]
                                ref_events_group <- sum(ref_data[[event_var]], na.rm = TRUE)
                            }
                        }
                    }
                }
                
                ## Create reference row with correct counts
                ref_row <- dt[1, ]
                ref_row[, `:=`(
                    term = paste0(factor_var, ref_level),
                    variable = factor_var,
                    group = ref_level,
                    n = n,
                    n_group = ref_n_group,
                    events = events,
                    events_group = ref_events_group,
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
    
    ## Reorder columns to put variable and level early
    col_order <- c("model_scope", "model_type", "term", "variable", "group")

    ## Add count columns if they exist
    for (col in c("n", "n_group", "events", "events_group")) {
        if (col %in% names(dt)) {
            col_order <- c(col_order, col)
        }
    }

    ## Add remaining columns
    other_cols <- setdiff(names(dt), col_order)
    setcolorder(dt, c(col_order, other_cols))

    ## Remove term column entirely
    dt[, term := NULL]
    
    ## Set attributes for model info
    data.table::setattr(dt, "model_class", model_class)
    data.table::setattr(dt, "model_family", if (model_class == "glm") model$family$family else NA)
    data.table::setattr(dt, "conf_level", conf_level)

    ## Force data.table to finalize
    dt[]
    
    return(dt)
}
