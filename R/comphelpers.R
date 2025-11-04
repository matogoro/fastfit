#' Check model convergence
#' @keywords internal
check_convergence <- function(model) {
    
    if (inherits(model, "glm")) {
        if (!model$converged) return("No")
        if (any(abs(coef(model)) > 10)) return("Suspect")
        return("Yes")
        
    } else if (inherits(model, c("coxph", "coxme"))) {
        
        ## Cox models don't always have an iter field
        ## Check for actual convergence issues
        if (!is.null(model$info) && !is.null(model$info$convergence)) {
            if (!model$info$convergence) return("No")
        }
        
        ## Check for extreme coefficients as a proxy for issues
        if (any(abs(coef(model)) > 10, na.rm = TRUE)) return("Suspect")
        
        ## Check if model failed (no coefficients)
        if (length(coef(model)) == 0) return("Failed")
        
        return("Yes")
        
    } else {
        
        return("Yes")
        
    }
}

#' Extract comprehensive model metrics based on academic consensus
#' @keywords internal
extract_model_metrics <- function(model, raw_data, model_type) {
    
    metrics <- list(
        n = raw_data$n[1],
        events = if ("events" %chin% names(raw_data)) raw_data$events[1] else NA_integer_,
        predictors = length(coef(model)),
        aic = if ("AIC" %chin% names(raw_data)) raw_data$AIC[1] else AIC(model),
        bic = if ("BIC" %chin% names(raw_data)) raw_data$BIC[1] else BIC(model),
        ## Initialize all possible metrics as NA to avoid missing list elements
        pseudo_r2 = NA_real_,
        concordance = NA_real_,
        global_p = NA_real_,
        null_deviance = NA_real_,
        residual_deviance = NA_real_,
        deviance_ratio = NA_real_,
        mcfadden_r2 = NA_real_,
        nagelkerke_r2 = NA_real_,
        tjur_r2 = NA_real_,
        c_statistic = NA_real_,
        brier_score = NA_real_,
        hoslem_p = NA_real_
    )
    
    ## GLM-specific metrics
    if (model_type == "glm") {
        ## Deviance metrics
        metrics$null_deviance <- model$null.deviance
        metrics$residual_deviance <- model$deviance
        metrics$deviance_ratio <- 1 - (model$deviance / model$null.deviance)
        
        if (model$family$family == "binomial") {
            ## Multiple pseudo R-squared measures
            n <- length(model$y)
            
            ## McFadden's R²
            null_model <- glm(model$y ~ 1, family = model$family)
            metrics$mcfadden_r2 <- as.numeric(1 - (logLik(model)/logLik(null_model)))
            metrics$pseudo_r2 <- metrics$mcfadden_r2  # Set generic pseudo_r2
            
            ## Nagelkerke's R² (Cox-Snell adjusted)
            cox_snell <- 1 - exp((model$deviance - model$null.deviance)/n)
            metrics$nagelkerke_r2 <- cox_snell / (1 - exp(-model$null.deviance/n))
            
            ## Tjur's R² (coefficient of discrimination)
            pred_probs <- fitted(model)
            metrics$tjur_r2 <- mean(pred_probs[model$y == 1]) - mean(pred_probs[model$y == 0])
            
            ## C-statistic (concordance)
            if (requireNamespace("pROC", quietly = TRUE)) {
                roc_obj <- pROC::roc(model$y, fitted(model), quiet = TRUE)
                metrics$c_statistic <- as.numeric(pROC::auc(roc_obj))
                metrics$concordance <- metrics$c_statistic  ## Set generic concordance
                
                ## Optimal threshold metrics
                coords <- pROC::coords(roc_obj, "best", ret = c("threshold", "sensitivity", "specificity"))
                metrics$optimal_threshold <- coords$threshold
                metrics$sensitivity <- coords$sensitivity
                metrics$specificity <- coords$specificity
            }
            
            ## Brier score
            metrics$brier_score <- mean((fitted(model) - model$y)^2)
            
            ## Hosmer-Lemeshow test (if ResourceSelection available)
            if (requireNamespace("ResourceSelection", quietly = TRUE) && n >= 40) {
                tryCatch({
                    hl <- ResourceSelection::hoslem.test(model$y, fitted(model), g = 10)
                    metrics$hoslem_p <- hl$p.value
                }, error = function(e) {
                    metrics$hoslem_p <- NA_real_
                })
            }
            
            ## Global LR test
            lr_stat <- model$null.deviance - model$deviance
            df <- model$df.null - model$df.residual
            metrics$global_p <- pchisq(lr_stat, df, lower.tail = FALSE)
        }
        
    } else if (model_type == "coxph") {
        
        ## Initialize Cox-specific metrics
        metrics$c_index <- NA_real_
        metrics$rsq <- NA_real_
        metrics$lr_test_p <- NA_real_
        
        ## Cache summary to avoid repeated calls
        summ <- summary(model)
        
        ## Concordance
        if (!is.null(model$concordance)) {
            metrics$c_index <- model$concordance["concordance"]
            metrics$concordance <- metrics$c_index
            
        } else if (!is.null(summ$concordance)) {
            metrics$c_index <- summ$concordance["C"]
            metrics$concordance <- metrics$c_index
        }
        
        ## R-squared (different location in Cox models)
        if (!is.null(summ$rsq)) {
            metrics$rsq <- summ$rsq[1]
            metrics$pseudo_r2 <- metrics$rsq
        }
        
        ## Global tests
        if (!is.null(summ$logtest)) {
            metrics$lr_test_p <- summ$logtest["pvalue"]
            metrics$global_p <- metrics$lr_test_p
            
        } else if (!is.null(summ$waldtest)) {
            metrics$global_p <- summ$waldtest["pvalue"]
        }
    }
    
    return(metrics)
}

#' Build comprehensive comparison table
#' @keywords internal
build_comparison_table <- function(comparison, model_type) {
    
    ## Select columns based on model type and academic consensus
    if (model_type == "glm") {
        ## Focus on discrimination, calibration, and information criteria
        key_cols <- c("Model", "N", "Events", "Predictors", "Converged",
                      "AIC", "BIC", "C-statistic", "Brier Score",
                      "McFadden R2", "Nagelkerke R2", 
                      "Hoslem p", "Global p")
        
        ## Rename for display - single call is more efficient
        old_names <- c("c_statistic", "brier_score", "mcfadden_r2", 
                       "nagelkerke_r2", "hoslem_p", "global_p")
        new_names <- c("C-statistic", "Brier Score", "McFadden R2", 
                       "Nagelkerke R2", "Hoslem p", "Global p")
        setnames(comparison, old_names, new_names, skip_absent = TRUE)
        
    } else if (model_type == "coxph") {
        key_cols <- c("Model", "N", "Events", "Predictors", "Converged",
                      "AIC", "BIC", "C-index", "R2", "R2 max",
                      "PH test p", "Global p")
        
        old_names <- c("c_index", "rsq", "rsq_max", "ph_global_p", "lr_test_p")
        new_names <- c("C-index", "R2", "R2 max", "PH test p", "Global p")
        setnames(comparison, old_names, new_names, skip_absent = TRUE)
        
    } else if (model_type == "lm") {
        key_cols <- c("Model", "N", "Predictors", "Converged",
                      "AIC", "BIC", "R2", "Adj R2", "RMSE",
                      "F-stat", "Global p")
        
        old_names <- c("r_squared", "adj_r_squared", "rmse", "f_statistic", "global_p")
        new_names <- c("R2", "Adj R2", "RMSE", "F-stat", "Global p")
        setnames(comparison, old_names, new_names, skip_absent = TRUE)
    }
    
    ## Keep only relevant columns that exist
    key_cols <- intersect(key_cols, names(comparison))
    comparison <- comparison[, ..key_cols]
    
    ## Format numeric columns appropriately
    format_model_comparison(comparison)
    
    return(comparison)
}

#' Format model comparison table
#' @keywords internal
format_model_comparison <- function(comparison) {
    
    ## Round numeric columns to appropriate precision
    numeric_cols <- names(comparison)[sapply(comparison, is.numeric)]
    
    for (col in numeric_cols) {
        if (col %chin% c("AIC", "BIC")) {
            comparison[, (col) := round(get(col), 1)]
        } else if (col %chin% c("C-statistic", "C-index", "Concordance", 
                                "McFadden R2", "Nagelkerke R2", "R2", "Adj R2",
                                "Brier Score", "RMSE")) {
            comparison[, (col) := round(get(col), 3)]
        }
    }
    
    return(comparison)
}

#' Calculate FastFit scores for model comparison
#' @keywords internal
calculate_model_scores <- function(comparison, model_type, scoring_weights = NULL) {
    
    ## Define default weights based on model type
    default_weights <- list(
        glm = list(convergence = 0.15, aic = 0.25, concordance = 0.40, 
                   pseudo_r2 = 0.15, brier = 0.05),
        coxph = list(convergence = 0.15, aic = 0.30, concordance = 0.40, 
                     global_p = 0.15),
        lm = list(convergence = 0.15, aic = 0.25, pseudo_r2 = 0.45, 
                  rmse = 0.15)
    )
    
    ## Use provided weights or defaults
    if (is.null(scoring_weights)) {
        weights <- default_weights[[model_type]]
        if (is.null(weights)) {
            weights <- list(convergence = 0.20, aic = 0.40, bic = 0.40)
        }
    } else {
        weights <- scoring_weights
    }
    
    ## Validate weights sum to 1
    weight_sum <- sum(unlist(weights))
    if (abs(weight_sum - 1) > 0.01) {
        warning("Scoring weights do not sum to 1, normalizing...")
        weights <- lapply(weights, function(x) x / weight_sum)
    }
    
    ## Initialize scores with pre-allocation
    n_models <- nrow(comparison)
    scores <- list(
        conv_score = numeric(n_models),
        aic_score = numeric(n_models),
        concordance_score = numeric(n_models),
        pseudo_r2_score = numeric(n_models),
        brier_score = numeric(n_models),
        global_score = numeric(n_models),
        rmse_score = numeric(n_models),
        bic_score = numeric(n_models),
        total = numeric(n_models)
    )
    
    ## Convergence score (universal)
    scores$conv_score <- data.table::fcase(
                                         comparison$Converged == "Yes", 100,
                                         comparison$Converged == "Suspect", 70,
                                         comparison$Converged == "No", 30,
                                         default = 0
                                     )
    
    ## AIC score (universal, lower is better)
    if (!all(is.na(comparison$AIC))) {
        aic_values <- comparison$AIC[!is.na(comparison$AIC)]
        if (length(aic_values) > 1) {
            aic_best <- min(aic_values)
            aic_worst <- max(aic_values)
            if (aic_worst - aic_best > 0) {
                scores$aic_score <- data.table::fifelse(
                                                    is.na(comparison$AIC), 0,
                                                    100 * (1 - (comparison$AIC - aic_best)/(aic_worst - aic_best))
                )
            } else {
                scores$aic_score <- rep(100, n_models)
            }
        } else {
            scores$aic_score <- data.table::fifelse(is.na(comparison$AIC), 0, 100)
        }
        scores$aic_score[is.na(scores$aic_score)] <- 0
    } else {
        scores$aic_score <- rep(50, n_models)
    }
    
    ## Model-specific scores
    if (model_type == "glm") {
        ## Concordance/C-statistic using fcase for better performance
        if ("Concordance" %chin% names(comparison) && !all(is.na(comparison$Concordance))) {
            scores$concordance_score <- data.table::fcase(
                                                        is.na(comparison$Concordance), 0,
                                                        comparison$Concordance <= 0.5, 0,
                                                        comparison$Concordance <= 0.6, 40 * (comparison$Concordance - 0.5)/0.1,
                                                        comparison$Concordance <= 0.7, 40 + 20 * (comparison$Concordance - 0.6)/0.1,
                                                        comparison$Concordance <= 0.8, 60 + 20 * (comparison$Concordance - 0.7)/0.1,
                                                        comparison$Concordance <= 0.9, 80 + 10 * (comparison$Concordance - 0.8)/0.1,
                                                        default = 90 + 10 * (comparison$Concordance - 0.9)/0.1
                                                    )
            scores$concordance_score <- pmin(100, scores$concordance_score)
        } else {
            scores$concordance_score <- rep(50, n_models)
        }
        
        ## Pseudo-R^2
        if ("Pseudo-R^2" %chin% names(comparison) && !all(is.na(comparison$`Pseudo-R^2`))) {
            ## McFadden's R2 rarely exceeds 0.4 for good models
            scores$pseudo_r2_score <- data.table::fifelse(
                is.na(comparison$`Pseudo-R^2`), 0,
                pmin(100, comparison$`Pseudo-R^2` * 250)  # 0.4 -> 100 points
            )
        } else {
            scores$pseudo_r2_score <- rep(50, n_models)
        }

        ## Brier score - only if column exists and weight is non-zero
        if ("Brier Score" %chin% names(comparison) && "brier" %chin% names(weights) && weights$brier > 0) {
            ## Lower is better (0 = perfect, 0.25 = no skill)
            scores$brier_score <- data.table::fifelse(
                is.na(comparison$`Brier Score`), 50,
                100 * (1 - comparison$`Brier Score`/0.25)
            )
        } else {
            scores$brier_score <- rep(50, n_models)
            weights$brier <- 0  # Set weight to 0 if not using
        }
        
        ## Normalize weights if Brier is excluded
        if (weights$brier == 0) {
            remaining_weights <- weights[names(weights) != "brier"]
            weight_sum <- sum(unlist(remaining_weights))
            remaining_weights <- lapply(remaining_weights, function(x) x / weight_sum)
            weights[names(remaining_weights)] <- remaining_weights
        }
        
        ## Calculate weighted total
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$concordance_score * weights$concordance +
            scores$pseudo_r2_score * weights$pseudo_r2
        
    } else if (model_type == "coxph") {
        ## Similar adjustments for Cox models
        if ("Concordance" %chin% names(comparison) && !all(is.na(comparison$Concordance))) {
            scores$concordance_score <- data.table::fcase(
                                                        is.na(comparison$Concordance), 0,
                                                        comparison$Concordance <= 0.5, 0,
                                                        default = pmin(100, 200 * (comparison$Concordance - 0.5))
                                                    )
        } else {
            scores$concordance_score <- rep(50, n_models)
        }
        
        ## Global p-value score
        if ("Global p" %chin% names(comparison)) {
            global_numeric <- suppressWarnings(as.numeric(gsub("< ", "", comparison$`Global p`)))
            scores$global_score <- data.table::fcase(
                                                   is.na(global_numeric), 50,
                                                   global_numeric > 0.05, 50,
                                                   default = 50 + 50 * (0.05 - global_numeric)/0.05
                                               )
        } else {
            scores$global_score <- rep(50, n_models)
        }
        
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$concordance_score * weights$concordance +
            scores$global_score * weights$global_p
        
    } else if (model_type == "lm") {
        ## R^2 scoring
        if ("Pseudo-R^2" %chin% names(comparison) && !all(is.na(comparison$`Pseudo-R^2`))) {
            scores$pseudo_r2_score <- comparison$`Pseudo-R^2` * 100
            scores$pseudo_r2_score[is.na(scores$pseudo_r2_score)] <- 0
        } else {
            scores$pseudo_r2_score <- rep(50, n_models)
        }
        
        ## RMSE scoring would go here if column exists
        scores$rmse_score <- rep(50, n_models)  # Placeholder
        
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$pseudo_r2_score * weights$pseudo_r2 +
            scores$rmse_score * weights$rmse
        
    } else {
        ## Generic fallback
        scores$bic_score <- rep(50, n_models)
        if (!all(is.na(comparison$BIC))) {
            bic_values <- comparison$BIC[!is.na(comparison$BIC)]
            if (length(bic_values) > 1) {
                bic_best <- min(bic_values)
                bic_worst <- max(bic_values)
                if (bic_worst - bic_best > 0) {
                    scores$bic_score <- 100 * (1 - (comparison$BIC - bic_best)/(bic_worst - bic_best))
                } else {
                    scores$bic_score <- rep(100, n_models)
                }
            } else {
                scores$bic_score <- rep(100, n_models)
            }
            scores$bic_score[is.na(scores$bic_score)] <- 0
        }
        
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$bic_score * weights$bic
    }
    
    ## Add final score to comparison table using := for efficiency
    comparison[, `FastFit Score` := round(scores$total, 1)]
    
    ## Sort by score (highest first)
    setorder(comparison, -`FastFit Score`)
    
    ## Store detailed scores as attribute
    setattr(comparison, "detailed_scores", scores)
    setattr(comparison, "weights", weights)
    
    return(comparison)
}

#' Combine coefficient tables from multiple models
#' @keywords internal
combine_coefficient_tables <- function(coef_list, model_names) {
    if (length(coef_list) == 0) return(NULL)
    
    ## Add model identifier to each table
    for (i in seq_along(coef_list)) {
        if (!is.null(coef_list[[i]])) {
            coef_list[[i]]$Model <- model_names[i]
        }
    }
    
    ## Combine all tables
    combined <- data.table::rbindlist(coef_list, fill = TRUE)
    
    ## Reorder columns to put Model first
    if ("Model" %in% names(combined)) {
        cols <- names(combined)
        new_order <- c("Model", setdiff(cols, "Model"))
        combined <- combined[, ..new_order]
    }
    
    return(combined)
}

#' Print method showing scoring methodology
#' @keywords internal
#' @export
print.compfit_result <- function(x, ...) {
    cat("\nModel Comparison Results\n")
    cat("Outcome: ", attr(x, "outcome"), "\n", sep = "")
    cat("Model Type: ", attr(x, "model_type"), "\n", sep = "")
    
    ## Show scoring weights
    weights <- attr(x, "weights")
    if (!is.null(weights)) {
        cat("\nFastFit Score Weights:\n")
        for (metric in names(weights)) {
            ## Formatting for metric names
            display_name <- switch(metric,
                                   "convergence" = "Convergence",
                                   "aic" = "AIC",
                                   "bic" = "BIC",
                                   "concordance" = "Concordance",
                                   "c_stat" = "C-statistic",
                                   "c_index" = "C-index",
                                   "pseudo_r2" = "Pseudo-R^2",
                                   "global_p" = "Global p-value",
                                   "rmse" = "RMSE",
                                   "brier" = "Brier score",
                                   "adj_r2" = "Adjusted R²",
                                   metric  # fallback
                                   )
            cat(sprintf("  %s: %.0f%%\n", display_name, weights[[metric]] * 100))
        }
    }
    
    ## Identify best model
    if (nrow(x) > 0) {
        best_model <- x$Model[1]  # Already sorted by score
        cat("\nRecommended Model: ", best_model, " (FastFit Score: ", x$`FastFit Score`[1], ")\n", sep = "")
    }
    
    cat("\nModels ranked by selection score:\n")
    NextMethod("print", x)

    cat("\nFastFit Score interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor\n")
    
    if (!is.null(attr(x, "coefficients"))) {
        cat("\nNote: Coefficient comparison available via attr(result, 'coefficients')\n")
    }
    
    invisible(x)
}
