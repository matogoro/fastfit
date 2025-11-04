#' Package Imports
#'
#' This file declares imports from base R and stats packages to avoid
#' "no visible global function definition" warnings in R CMD check.
#'
#' @name fastfit-imports
#' @importFrom stats sd var median quantile IQR
#' @importFrom stats logLik AIC BIC coef confint vcov
#' @importFrom stats chisq.test fisher.test t.test wilcox.test kruskal.test
#' @importFrom stats lm glm anova
#' @importFrom stats family binomial gaussian poisson
#' @importFrom stats predict fitted residuals
#' @importFrom stats formula terms model.frame model.matrix
#' @importFrom stats na.omit complete.cases
#' @importFrom stats cor pnorm qnorm
#' @importFrom stats nobs
#' @importFrom grDevices axisTicks
#' @importFrom utils head tail str
#' @keywords internal
NULL
