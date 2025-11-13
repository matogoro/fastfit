#' Convert Model to Data Table
#'
#' Extracts coefficients, confidence intervals, and comprehensive model statistics 
#' from fitted regression models and converts them to a standardized data.table 
#' format suitable for further analysis or publication. This is a core utility 
#' function in the fastfit package for working with regression results.
#'
#' @param model A fitted model object. Supported classes include:
#'   \itemize{
#'     \item \code{glm} - Generalized linear models (logistic, Poisson, etc.)
#'     \item \code{lm} - Linear models
#'     \item \code{coxph} - Cox proportional hazards models
#'     \item \code{clogit} - Conditional logistic regression
#'     \item \code{coxme} - Mixed effects Cox models
#'     \item \code{glmer} - Generalized linear mixed effects models
#'   }
#'   Also accepts models wrapped with \code{mmodel()}.
#'   
#' @param conf_level Numeric confidence level for confidence intervals. Must be
#'   between 0 and 1. Default is 0.95 (95\% CI).
#'   
#' @param keep_qc_stats Logical. If \code{TRUE}, includes model quality statistics 
#'   such as AIC, BIC, R-squared, concordance, and model fit tests. These appear 
#'   as additional columns in the output. Default is \code{TRUE}.
#'   
#' @param include_intercept Logical. If \code{TRUE}, includes the model intercept 
#'   in output. If \code{FALSE}, removes the intercept row from results. Useful 
#'   for creating cleaner presentation tables. Default is \code{TRUE}.
#'   
#' @param terms_to_exclude Character vector of term names to exclude from output.
#'   Useful for removing specific unwanted parameters (e.g., nuisance variables,
#'   spline terms). Default is \code{NULL}. Note: If \code{include_intercept = FALSE}, 
#'   "(Intercept)" is automatically added to this list.
#'   
#' @param add_reference_rows Logical. If \code{TRUE}, adds rows for reference 
#'   categories of factor variables with appropriate labels and baseline values 
#'   (OR/HR = 1, Estimate = 0). This makes tables more complete and easier to 
#'   interpret. Default is \code{TRUE}.
#'   
#' @param reference_label Character string used to label reference category rows
#'   in the output. Appears in the \code{reference} column. Default is \code{"reference"}.
#'
#' @return A \code{data.table} containing extracted model information with the 
#'   following standard columns:
#'   \describe{
#'     \item{model_scope}{Character. Either "Univariable" (unadjusted model with 
#'       single predictor) or "Multivariable" (adjusted model with multiple predictors)}
#'     \item{model_type}{Character. Type of regression (e.g., "Logistic", "Linear", 
#'       "Cox PH", "Poisson")}
#'     \item{variable}{Character. Variable name (for factor variables, the base 
#'       variable name without the level)}
#'     \item{group}{Character. Group/level name for factor variables; empty string 
#'       for continuous variables}
#'     \item{n}{Integer. Total sample size used in the model}
#'     \item{n_group}{Integer. Sample size for this specific variable level 
#'       (factor variables only)}
#'     \item{events}{Integer. Total number of events in the model (for survival 
#'       and logistic models)}
#'     \item{events_group}{Integer. Number of events for this specific variable 
#'       level (for survival and logistic models with factor variables)}
#'     \item{coefficient}{Numeric. Raw regression coefficient (log odds, log hazard, 
#'       etc.)}
#'     \item{se}{Numeric. Standard error of the coefficient}
#'     \item{OR/HR/RR/Estimate}{Numeric. Effect estimate - column name depends on 
#'       model type:
#'       \itemize{
#'         \item \code{OR} for logistic regression (odds ratio)
#'         \item \code{HR} for Cox models (hazard ratio)
#'         \item \code{RR} for Poisson regression (rate/risk ratio)
#'         \item \code{Estimate} for linear models or other GLMs
#'       }}
#'     \item{CI_lower}{Numeric. Lower bound of confidence interval for effect estimate}
#'     \item{CI_upper}{Numeric. Upper bound of confidence interval for effect estimate}
#'     \item{statistic}{Numeric. Test statistic (z-value for GLM/Cox, t-value for LM)}
#'     \item{p_value}{Numeric. P-value for coefficient test}
#'     \item{sig}{Character. Significance markers: "***" (p<0.001), "**" (p<0.01), 
#'       "*" (p<0.05), "." (p<0.10), "" (p≥0.10)}
#'     \item{sig_binary}{Logical. Binary indicator: \code{TRUE} if p<0.05, 
#'       \code{FALSE} otherwise}
#'     \item{reference}{Character. Contains \code{reference_label} for reference 
#'       category rows when \code{add_reference_rows = TRUE}, empty string otherwise}
#'   }
#'   
#'   Additional quality control columns when \code{keep_qc_stats = TRUE} (vary by 
#'   model type):
#'   \describe{
#'     \item{GLM models}{AIC, BIC, deviance, null_deviance, df_residual. For logistic 
#'       regression: c_statistic (concordance/AUC)}
#'     \item{LM models}{R2 (R-squared), adj_R2 (adjusted R-squared), sigma (residual 
#'       standard error), df_residual}
#'     \item{Cox models}{concordance (C-index), rsq (R-squared approximation), 
#'       logtest_stat, logtest_p (likelihood ratio test), wald_test, wald_p 
#'       (Wald test), score_test, score_p (score test)}
#'   }
#'   
#'   The output includes the following attributes accessible via \code{attr()}:
#'   \describe{
#'     \item{model_class}{Character. The class of the input model object}
#'     \item{model_family}{Character. The family for GLM models (e.g., "binomial", 
#'       "gaussian"); \code{NA} for non-GLM models}
#'     \item{conf_level}{Numeric. The confidence level used for intervals}
#'   }
#'
#' @details
#' This function automatically detects whether a model is univariable or multivariable 
#' by counting the number of unique predictor variables (accounting for factor 
#' variable expansion). It handles factor variables intelligently by:
#' \itemize{
#'   \item Parsing factor terms into base variable names and levels
#'   \item Computing group-specific sample sizes and event counts
#'   \item Optionally adding reference category rows with OR/HR = 1
#'   \item Maintaining factor level ordering from the original model
#' }
#' 
#' For logistic regression (\code{glm} with \code{family = binomial}), the function 
#' returns odds ratios (OR). For Cox models, it returns hazard ratios (HR). For 
#' Poisson regression, it returns rate/risk ratios (RR). For linear models and 
#' other GLMs, it returns the raw coefficient estimates.
#' 
#' The function calculates confidence intervals using the normal approximation 
#' (z-score method) for all model types, which is appropriate for large samples.
#' 
#' When \code{add_reference_rows = TRUE}, the function adds rows for reference 
#' categories of factor variables with:
#' \itemize{
#'   \item Effect estimate = 1 (for OR/HR/RR) or 0 (for Estimate)
#'   \item Coefficient = 0
#'   \item Standard error, statistic, p-value, and CI bounds = NA
#'   \item Group-specific sample sizes and event counts calculated from the data
#'   \item The \code{reference} column populated with \code{reference_label}
#' }
#'
#' @seealso 
#' \code{\link{glm}}, \code{\link{lm}}, \code{\link[survival]{coxph}}, 
#' \code{\link[survival]{clogit}}, \code{\link[lme4]{glmer}}
#' 
#' @examples
#' # Load example data
#' data(clintrial)
#' 
#' # Example 1: Simple logistic regression
#' model1 <- glm(os_status ~ age + sex, 
#'               data = clintrial, 
#'               family = binomial)
#' result1 <- m2dt(model1)
#' print(result1)
#' 
#' # Example 2: Remove intercept from output
#' result2 <- m2dt(model1, include_intercept = FALSE)
#' print(result2)
#' 
#' # Example 3: Logistic regression with factor variable
#' model3 <- glm(os_status ~ age + treatment, 
#'               data = clintrial, 
#'               family = binomial)
#' result3 <- m2dt(model3, add_reference_rows = TRUE)
#' # Note: reference category for treatment is included with OR = 1
#' print(result3)
#' 
#' # Example 4: Custom confidence level and reference label
#' result4 <- m2dt(model3, 
#'                 conf_level = 0.90,
#'                 reference_label = "ref",
#'                 add_reference_rows = TRUE)
#' print(result4)
#' 
#' # Example 5: Cox proportional hazards model
#' library(survival)
#' model5 <- coxph(Surv(os_months, os_status) ~ age + sex + treatment, 
#'                 data = clintrial)
#' result5 <- m2dt(model5)
#' # Returns hazard ratios (HR) instead of odds ratios
#' print(result5)
#' 
#' # Example 6: Linear model with QC statistics
#' model6 <- lm(bmi ~ age + sex + smoking, data = clintrial)
#' result6 <- m2dt(model6, keep_qc_stats = TRUE)
#' # Includes R-squared, AIC, BIC, etc.
#' print(result6)
#' 
#' # Example 7: Exclude specific terms
#' model7 <- glm(os_status ~ age + sex + treatment + site, 
#'               data = clintrial, 
#'               family = binomial)
#' result7 <- m2dt(model7, 
#'                 terms_to_exclude = c("siteB", "siteC"),
#'                 include_intercept = FALSE)
#' print(result7)
#' 
#' # Example 8: Access model attributes
#' result8 <- m2dt(model1)
#' attr(result8, "model_class")    # "glm"
#' attr(result8, "model_family")   # "binomial"
#' attr(result8, "conf_level")     # 0.95
#' 
#' # Example 9: Univariable vs Multivariable detection
#' uni_model <- glm(os_status ~ age, 
#'                  data = clintrial, 
#'                  family = binomial)
#' uni_result <- m2dt(uni_model)
#' print(uni_result$model_scope)  # "Univariable"
#' 
#' multi_model <- glm(os_status ~ age + sex + treatment, 
#'                    data = clintrial, 
#'                    family = binomial)
#' multi_result <- m2dt(multi_model)
#' print(multi_result$model_scope)  # "Multivariable"
#' 
#' # Example 10: Poisson regression for count data
#' model10 <- glm(los_days ~ age + treatment, 
#'                data = clintrial, 
#'                family = poisson)
#' result10 <- m2dt(model10)
#' # Returns rate ratios (RR)
#' print(result10)
#' 
#' @export
m2dt <- function(model, 
                 conf_level = 0.95,
                 keep_qc_stats = TRUE,
                 include_intercept = TRUE,
                 terms_to_exclude = NULL,
                 add_reference_rows = TRUE,
                 reference_label = "reference") {
    
    ## Handle intercept exclusion
    if (!include_intercept) {
        if (is.null(terms_to_exclude)) {
            terms_to_exclude <- "(Intercept)"
        } else {
            terms_to_exclude <- unique(c(terms_to_exclude, "(Intercept)"))
        }
    }
    
    ## If terms_to_exclude is still NULL at this point, set it to empty character vector (no exclusions)
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
            coef_vec <- coxme::fixef(model)
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
    
                                        # Filter out excluded terms before any further processing
    if (length(terms_to_exclude) > 0) {
        dt <- dt[!term %in% terms_to_exclude]
    }

    ## Process terms to extract variable and group information - OPTIMIZED
    if (!("variable" %in% names(dt))) {
        ## Use vectorized parse_term helper function
        parsed <- parse_term(dt$term, model$xlevels)
        dt[, `:=`(variable = parsed$variable, group = parsed$group)]
        
        ## Add n_group and events_group for categorical variables - OPTIMIZED
        if (!is.null(model$xlevels)) {
            
            ## Initialize columns
            dt[, `:=`(n_group = NA_real_, events_group = NA_real_)]
            
            ## Get data source
            data_source <- NULL
            if (!is.null(model$data)) {
                data_source <- model$data
            } else if (!is.null(model$model)) {
                data_source <- model$model
            }
            
            if (!is.null(data_source)) {
                ## Convert to data.table for efficient operations
                data_dt <- data.table::as.data.table(data_source)
                
                ## Get outcome variable name once
                outcome_var <- NULL
                event_var <- NULL
                
                if (model_class == "glm" && model$family$family == "binomial") {
                    outcome_var <- all.vars(model$formula)[1]
                } else if (model_class %in% c("coxph", "clogit")) {
                    outcome_str <- as.character(model$formula)[2]
                    if (grepl("^Surv\\(", outcome_str)) {
                        surv_expr <- gsub("Surv\\(|\\)", "", outcome_str)
                        surv_parts <- trimws(strsplit(surv_expr, ",")[[1]])
                        event_var <- surv_parts[2]
                    }
                }
                
                ## Process all factor variables at once - VECTORIZED
                for (var in names(model$xlevels)) {
                    if (var %in% names(data_dt)) {
                        ## Vectorized count calculation using data.table
                        if (!is.null(outcome_var) && outcome_var %in% names(data_dt)) {
                            ## For binomial: count n and events by level
                            outcome_col <- data_dt[[outcome_var]]
                            
                            ## Handle factor outcomes
                            if (is.factor(outcome_col)) {
                                data_dt[, .events_calc := as.numeric(get(outcome_var)) == 2]
                            } else {
                                data_dt[, .events_calc := get(outcome_var)]
                            }
                            
                            ## Aggregate in one pass
                            counts <- data_dt[!is.na(get(var)), .(
                                                                    n_group = .N,
                                                                    events_group = sum(.events_calc, na.rm = TRUE)
                                                                ), by = var]
                            
                            ## Clean up temporary column
                            data_dt[, .events_calc := NULL]
                            
                            ## Set names for joining
                            data.table::setnames(counts, var, "group")
                            counts[, variable := var]
                            
                            ## Join back to main dt using data.table update join
                            dt[counts, `:=`(
                                           n_group = i.n_group,
                                           events_group = i.events_group
                                       ), on = .(variable, group)]
                            
                        } else if (!is.null(event_var) && event_var %in% names(data_dt)) {
                            ## For survival: count n and events by level
                            counts <- data_dt[!is.na(get(var)), .(
                                                                    n_group = .N,
                                                                    events_group = sum(get(event_var), na.rm = TRUE)
                                                                ), by = var]
                            
                            data.table::setnames(counts, var, "group")
                            counts[, variable := var]
                            
                            dt[counts, `:=`(
                                           n_group = i.n_group,
                                           events_group = i.events_group
                                       ), on = .(variable, group)]
                            
                        } else {
                            ## Just count n by level
                            counts <- data_dt[!is.na(get(var)), .N, by = var]
                            data.table::setnames(counts, c("group", "n_group"))
                            counts[, variable := var]
                            
                            dt[counts, n_group := i.n_group, on = .(variable, group)]
                        }
                    }
                }
            }
        }
    }

    ## Add reference rows for factor variables while maintaining original order
    if (add_reference_rows && !is.null(model$xlevels)) {
        
        ## Add reference column to existing data
        dt[, reference := ""]
        
        ## Pre-calculate all reference level counts at once for efficiency
        ref_counts <- list()
        
        data_source <- NULL
        if (!is.null(model$data)) {
            data_source <- model$data
        } else if (!is.null(model$model)) {
            data_source <- model$model
        }
        
        if (!is.null(data_source)) {
            data_dt <- data.table::as.data.table(data_source)
            
            ## Get outcome/event info once
            outcome_var <- NULL
            event_var <- NULL
            
            if (model_class == "glm" && model$family$family == "binomial") {
                outcome_var <- all.vars(model$formula)[1]
            } else if (model_class %in% c("coxph", "clogit")) {
                outcome_str <- as.character(model$formula)[2]
                if (grepl("^Surv\\(", outcome_str)) {
                    surv_expr <- gsub("Surv\\(|\\)", "", outcome_str)
                    surv_parts <- trimws(strsplit(surv_expr, ",")[[1]])
                    event_var <- surv_parts[2]
                }
            }
            
            ## Calculate reference counts for all factors at once
            for (var in names(model$xlevels)) {
                if (var %in% names(data_dt)) {
                    ref_level <- model$xlevels[[var]][1]
                    
                    if (!is.null(outcome_var) && outcome_var %in% names(data_dt)) {
                        outcome_col <- data_dt[[outcome_var]]
                        if (is.factor(outcome_col)) {
                            data_dt[, .events_calc := as.numeric(get(outcome_var)) == 2]
                        } else {
                            data_dt[, .events_calc := get(outcome_var)]
                        }
                        
                        ref_counts[[var]] <- data_dt[get(var) == ref_level & !is.na(get(var)), .(
                                                                                                   n_group = .N,
                                                                                                   events_group = sum(.events_calc, na.rm = TRUE)
                                                                                               )]
                        
                        data_dt[, .events_calc := NULL]
                        
                    } else if (!is.null(event_var) && event_var %in% names(data_dt)) {
                        ref_counts[[var]] <- data_dt[get(var) == ref_level & !is.na(get(var)), .(
                                                                                                   n_group = .N,
                                                                                                   events_group = sum(get(event_var), na.rm = TRUE)
                                                                                               )]
                        
                    } else {
                        ref_counts[[var]] <- data_dt[get(var) == ref_level & !is.na(get(var)), .(
                                                                                                   n_group = .N
                                                                                               )]
                    }
                }
            }
        }
        
        ## Build final table while maintaining original order
        ## Create a list to hold the final ordered result
        final_rows <- list()
        
        ## Track which terms have been processed
        processed_terms <- character(0)
        
        ## Process each term in the original order
        for (i in seq_len(nrow(dt))) {
            current_term <- dt$term[i]
            
            ## Skip if already processed
            if (current_term %in% processed_terms) next
            
            ## Interaction terms contain ":" and should be treated as continuous-like terms
            if (grepl(":", current_term, fixed = TRUE)) {
                ## This is an interaction term - add it as-is
                final_rows[[length(final_rows) + 1]] <- dt[term == current_term]
                processed_terms <- c(processed_terms, current_term)
                next
            }
            
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
                
                ## Get reference counts from pre-calculated values
                ref_n_group <- NA_real_
                ref_events_group <- NA_real_
                if (!is.null(ref_counts[[factor_var]])) {
                    if (nrow(ref_counts[[factor_var]]) > 0) {
                        ref_n_group <- ref_counts[[factor_var]]$n_group[1]
                        if ("events_group" %in% names(ref_counts[[factor_var]])) {
                            ref_events_group <- ref_counts[[factor_var]]$events_group[1]
                        }
                    }
                }
                
                ## Create reference row with correct counts
                ref_row <- dt[1, ]
                ref_row[, `:=`(
                    term = paste0(factor_var, ref_level),
                    variable = factor_var,
                    group = ref_level,
                    n_group = ref_n_group,
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
                final_rows[[length(final_rows) + 1]] <- ref_row
                
                ## Now add all other levels for this factor in order
                for (level in levels_order[-1]) {
                    level_term <- paste0(factor_var, level)
                    level_row <- dt[term == level_term]
                    if (nrow(level_row) > 0) {
                        final_rows[[length(final_rows) + 1]] <- level_row
                        processed_terms <- c(processed_terms, level_term)
                    }
                }
                
            } else {
                ## This is not a factor variable (e.g., continuous)
                final_rows[[length(final_rows) + 1]] <- dt[term == current_term]
                processed_terms <- c(processed_terms, current_term)
            }
        }
        
        ## Combine all rows
        dt <- data.table::rbindlist(final_rows, fill = TRUE)
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
    data.table::setcolorder(dt, c(col_order, other_cols))

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
