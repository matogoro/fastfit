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
        events = if ("events" %in% names(raw_data)) raw_data$events[1] else NA_integer_,
        predictors = length(coef(model)),
        aic = if ("AIC" %in% names(raw_data)) raw_data$AIC[1] else AIC(model),
        bic = if ("BIC" %in% names(raw_data)) raw_data$BIC[1] else BIC(model),
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
        
        ## Concordance
        if (!is.null(model$concordance)) {
            metrics$c_index <- model$concordance["concordance"]
            metrics$concordance <- metrics$c_index
            
        } else if (!is.null(summary(model)$concordance)) {
            metrics$c_index <- summary(model)$concordance["C"]
            metrics$concordance <- metrics$c_index
        }
        
        ## R-squared (different location in Cox models)
        summ <- summary(model)
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
        
        ## Rename for display
        setnames(comparison, 
                 c("c_statistic", "brier_score", "mcfadden_r2", 
                   "nagelkerke_r2", "hoslem_p", "global_p"),
                 c("C-statistic", "Brier Score", "McFadden R2", 
                   "Nagelkerke R2", "Hoslem p", "Global p"),
                 skip_absent = TRUE)
        
    } else if (model_type == "coxph") {
        key_cols <- c("Model", "N", "Events", "Predictors", "Converged",
                      "AIC", "BIC", "C-index", "R2", "R2 max",
                      "PH test p", "Global p")
        
        setnames(comparison,
                 c("c_index", "rsq", "rsq_max", "ph_global_p", "lr_test_p"),
                 c("C-index", "R2", "R2 max", "PH test p", "Global p"),
                 skip_absent = TRUE)
        
    } else if (model_type == "lm") {
        key_cols <- c("Model", "N", "Predictors", "Converged",
                      "AIC", "BIC", "R2", "Adj R2", "RMSE",
                      "F-stat", "Global p")
        
        setnames(comparison,
                 c("r_squared", "adj_r_squared", "rmse", "f_statistic", "global_p"),
                 c("R2", "Adj R2", "RMSE", "F-stat", "Global p"),
                 skip_absent = TRUE)
    }
    
    ## Keep only relevant columns that exist
    key_cols <- intersect(key_cols, names(comparison))
    comparison <- comparison[, ..key_cols]
    
    ## Format numeric columns appropriately
    format_model_comparison(comparison)
    
    return(comparison)
}

#' Format comparison table for display
#' @keywords internal  
format_model_comparison <- function(dt) {
    ## Round appropriately based on metric type
    
    ## 3 decimal places for proportions/probabilities
    for (col in c("C-statistic", "C-index", "Brier Score", "R2", "Adj R2", 
                  "McFadden R2", "Nagelkerke R2", "Tjur R2", "R2 max")) {
        if (col %in% names(dt)) {
            dt[[col]] <- round(dt[[col]], 3)
        }
    }
    
    ## 1 decimal for information criteria
    for (col in c("AIC", "BIC", "AICc")) {
        if (col %in% names(dt)) {
            dt[[col]] <- round(dt[[col]], 1)
        }
    }
    
    ## Format p-values
    for (col in grep("\\sp$|^Global p$|^Hoslem p$|^PH test p$", names(dt), value = TRUE)) {
        dt[[col]] <- format_pvalue(dt[[col]], 3)
    }
    
    return(dt)
}

#' Calculate model selection scores based on multiple criteria
#' @keywords internal
calculate_model_scores <- function(comparison, model_type, custom_weights = NULL) {
    
    ## Initialize scoring data.table
    scores <- data.table::data.table(Model = comparison$Model)
    
    ## Define default weights based on model type
    default_weights <- list(
        glm = list(
            convergence = 0.15,
            aic = 0.25,
            concordance = 0.30,
            pseudo_r2 = 0.20,
            brier = 0.10
        ),
        coxph = list(
            convergence = 0.20,
            aic = 0.30,
            concordance = 0.30,
            global_p = 0.20
        ),
        lm = list(
            convergence = 0.15,
            aic = 0.25,
            pseudo_r2 = 0.35,
            rmse = 0.25
        ),
        generic = list(
            convergence = 0.25,
            aic = 0.50,
            bic = 0.25
        )
    )
    
    ## Use custom weights if provided, otherwise use defaults
    if (!is.null(custom_weights)) {
        
        ## Validate custom weights
        if (!is.list(custom_weights)) {
            stop("scoring_weights must be a named list")
        }
        
        ## Check that weights sum to approximately 1
        weight_sum <- sum(unlist(custom_weights))
        if (abs(weight_sum - 1) > 0.01) {
            warning(sprintf("Scoring weights sum to %.2f, not 1.0. Normalizing...", weight_sum))
            custom_weights <- lapply(custom_weights, function(x) x / weight_sum)
        }
        
        ## Use custom weights, filling in any missing with defaults
        base_weights <- default_weights[[model_type]] %||% default_weights$generic
        weights <- modifyList(base_weights, custom_weights)
    } else {
        weights <- default_weights[[model_type]] %||% default_weights$generic
    }
    
    ## Score 1: Convergence (0-100)
    scores$conv_score <- ifelse(comparison$Converged == "Yes", 100,
                         ifelse(comparison$Converged == "Suspect", 50, 0))
    
    ## Score 2: AIC - lower is better
    if (!all(is.na(comparison$AIC))) {
        aic_values <- comparison$AIC[!is.na(comparison$AIC)]
        if (length(aic_values) > 1) {
            ## Use relative scoring - best model gets 100, others proportionally less
            aic_best <- min(aic_values)
            aic_worst <- max(aic_values)
            ## Avoid division by zero
            if (aic_worst - aic_best > 0) {
                scores$aic_score <- 100 * (1 - (comparison$AIC - aic_best)/(aic_worst - aic_best))
            } else {
                scores$aic_score <- 100  ## All models have same AIC
            }
        } else {
            scores$aic_score <- 100  ## Single model
        }
        scores$aic_score[is.na(scores$aic_score)] <- 0
    } else {
        scores$aic_score <- 50
    }
    
    ## Model-specific scores
    if (model_type == "glm") {
        ## Concordance/C-statistic
        ## 0.5 = no discrimination (0 points)
        ## 0.6 = poor (40 points)
        ## 0.7 = acceptable (60 points)
        ## 0.8 = good (80 points)
        ## 0.9+ = excellent (90-100 points)
        if ("Concordance" %in% names(comparison) && !all(is.na(comparison$Concordance))) {
            scores$concordance_score <- ifelse(
                is.na(comparison$Concordance), 0,
                                        ifelse(comparison$Concordance <= 0.5, 0,
                                        ifelse(comparison$Concordance <= 0.6, 40 * (comparison$Concordance - 0.5)/0.1,
                                        ifelse(comparison$Concordance <= 0.7, 40 + 20 * (comparison$Concordance - 0.6)/0.1,
                                        ifelse(comparison$Concordance <= 0.8, 60 + 20 * (comparison$Concordance - 0.7)/0.1,
                                        ifelse(comparison$Concordance <= 0.9, 80 + 10 * (comparison$Concordance - 0.8)/0.1,
                                               90 + 10 * (comparison$Concordance - 0.9)/0.1)))))
            )
            scores$concordance_score <- pmin(100, scores$concordance_score)  # Cap at 100
        } else {
            scores$concordance_score <- 50
        }
        
        ## Pseudo-R^2
        ## 0 to 0.5 is typical range
        if ("Pseudo-R^2" %in% names(comparison) && !all(is.na(comparison$`Pseudo-R^2`))) {
            ## McFadden's R2 rarely exceeds 0.4 for good models
            scores$pseudo_r2_score <- ifelse(
                is.na(comparison$`Pseudo-R^2`), 0,
                pmin(100, comparison$`Pseudo-R^2` * 250)  # 0.4 -> 100 points
            )
        } else {
            scores$pseudo_r2_score <- 50
        }

        ## Brier score - only if column exists and weight is non-zero
        if ("Brier Score" %in% names(comparison) && "brier" %in% names(weights) && weights$brier > 0) {
            ## Lower is better (0 = perfect, 0.25 = no skill)
            scores$brier_score <- ifelse(
                is.na(comparison$`Brier Score`), 50,
                100 * (1 - comparison$`Brier Score`/0.25)
            )
        } else {
            scores$brier_score <- 50
            weights$brier <- 0  # Set weight to 0 if not using
        }
        
        ## Normalize weights if Brier is excluded
        if (weights$brier == 0) {
            remaining_weights <- weights[names(weights) != "brier"]
            weight_sum <- sum(unlist(remaining_weights))
            for (w in names(remaining_weights)) {
                weights[[w]] <- weights[[w]] / weight_sum
            }
        }
        
        ## Calculate weighted total
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$concordance_score * weights$concordance +
            scores$pseudo_r2_score * weights$pseudo_r2
        
    } else if (model_type == "coxph") {
        ## Similar adjustments for Cox models
        if ("Concordance" %in% names(comparison) && !all(is.na(comparison$Concordance))) {
            scores$concordance_score <- ifelse(
                is.na(comparison$Concordance), 0,
                                        ifelse(comparison$Concordance <= 0.5, 0,
                                               pmin(100, 200 * (comparison$Concordance - 0.5)))
            )
        } else {
            scores$concordance_score <- 50
        }
        
        ## Global p-value score
        if ("Global p" %in% names(comparison)) {
            global_numeric <- suppressWarnings(as.numeric(gsub("< ", "", comparison$`Global p`)))
            scores$global_score <- ifelse(
                is.na(global_numeric), 50,
                                   ifelse(global_numeric > 0.05, 50,  # Not significant
                                          50 + 50 * (0.05 - global_numeric)/0.05)  # More significant = higher score
            )
        } else {
            scores$global_score <- 50
        }
        
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$concordance_score * weights$concordance +
            scores$global_score * weights$global_p
        
    } else if (model_type == "lm") {
        ## R^2 scoring
        if ("Pseudo-R^2" %in% names(comparison) && !all(is.na(comparison$`Pseudo-R^2`))) {
            scores$pseudo_r2_score <- comparison$`Pseudo-R^2` * 100
            scores$pseudo_r2_score[is.na(scores$pseudo_r2_score)] <- 0
        } else {
            scores$pseudo_r2_score <- 50
        }
        
        ## RMSE scoring would go here if column exists
        scores$rmse_score <- 50  # Placeholder
        
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$pseudo_r2_score * weights$pseudo_r2 +
            scores$rmse_score * weights$rmse
        
    } else {
        ## Generic fallback
        scores$bic_score <- 50
        if (!all(is.na(comparison$BIC))) {
            bic_values <- comparison$BIC[!is.na(comparison$BIC)]
            if (length(bic_values) > 1) {
                bic_best <- min(bic_values)
                bic_worst <- max(bic_values)
                if (bic_worst - bic_best > 0) {
                    scores$bic_score <- 100 * (1 - (comparison$BIC - bic_best)/(bic_worst - bic_best))
                } else {
                    scores$bic_score <- 100
                }
            } else {
                scores$bic_score <- 100
            }
            scores$bic_score[is.na(scores$bic_score)] <- 0
        }
        
        scores$total <- scores$conv_score * weights$convergence +
            scores$aic_score * weights$aic +
            scores$bic_score * weights$bic
    }
    
    ## Add final score to comparison table
    comparison$`Composite Score` <- round(scores$total, 1)
    
    ## Sort by score (highest first)
    setorder(comparison, -`Composite Score`)
    
    ## Store detailed scores as attribute
    setattr(comparison, "detailed_scores", scores)
    setattr(comparison, "weights", weights)
    
    return(comparison)
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
    ## In print.compfit_result function:
    if (!is.null(weights)) {
        cat("\nComposite Score Weights:\n")
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
        cat("\nRecommended Model: ", best_model, " (Composite Score: ", x$`Composite Score`[1], ")\n", sep = "")
    }
    
    cat("\nModels ranked by selection score:\n")
    NextMethod("print", x)

    cat("\nScore interpretation: 85+ Excellent, 75-84 Very Good, 65-74 Good, 55-64 Fair, < 55 Poor\n")
    
    if (!is.null(attr(x, "coefficients"))) {
        cat("\nNote: Coefficient comparison available via attr(result, 'coefficients')\n")
    }
    
    invisible(x)
}
