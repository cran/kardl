## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)

## ----srr-tags, eval = FALSE, echo = FALSE-------------------------------------
# #' srr tags for the kardl package
# #'
# #' This package follows structured statistical reporting rules using the
# #' \code{srrstats} framework. The standards governing implementation and
# #' documentation are summarized below.
# #'
# #' @srrstats {G1.0} The package documentation and accompanying materials cite
# #'           the ARDL bounds-testing and nonlinear ARDL literature on which the
# #'           implemented estimators, asymmetry decomposition, bounds tests, and
# #'           dynamic multipliers are based.
# #' @srrstats {G1.1} The package is described as an implementation and extension
# #'           of ARDL and NARDL workflows, supporting mixed symmetric and
# #'           asymmetric regressors, flexible lag selection, and dynamic
# #'           multiplier methods.
# #' @srrstats {G1.2} The README, NEWS file, and package website describe the
# #'           development status, recent changes, and future maintenance plans
# #'           of the package.
# #' @srrstats {G1.3} Key statistical concepts such as ARDL, NARDL, ECM, bounds
# #'           testing, short-run asymmetry, long-run asymmetry, and dynamic
# #'           multipliers are defined in function documentation and vignettes.
# #' @srrstats {G1.4} All user-facing functions and S3 methods are documented
# #'           using \code{roxygen2}, and the generated Rd files are maintained
# #'           as part of the package documentation workflow.
# #' @srrstats {G1.4a} Internal helper functions are documented where necessary
# #'           and are excluded from the public help index using \code{@noRd}
# #'           or marked as internal.

## ----install-cran, eval=FALSE-------------------------------------------------
# install.packages("kardl")
# library(kardl)

## ----install, eval=FALSE------------------------------------------------------
# # Install required packages
# install.packages(c(
#   "stats", "msm", "lmtest", "nlWaldTest", "car", "strucchange",
#   "utils", "ggplot2"
# ))
# # Install kardl from GitHub
# install.packages("devtools")
# devtools::install_github("karamelikli/kardl")

## ----load---------------------------------------------------------------------
library(kardl)

## ----data-prepare-------------------------------------------------------------
# Define the model formula
my_formula <- CPI ~ ER + PPI + asymmetric(ER + PPI) + deterministic(covid) +
  trend

## ----eval=FALSE---------------------------------------------------------------
# same_formula <- y ~ asymmetric(x1) +
#   sasymmetric(x2 + x3) +
#   lasymmetric(x4 + x5) +
#   deterministic(dummy1) + trend
# same_formula <- y ~ asymmetric(x1) +
#   sasymmetric(x2 + x3) +
#   lasymmetric(x4 + x5) +
#   deterministic(dummy1) + trend
# same_formula <- y ~ asym(x1) + sasym(x2 + x3) + lasym(x4 + x5) +
#   det(dummy1) + trend
# same_formula <- y ~ a(x1) + s(x2 + x3) + l(x4 + x5) + d(dummy1) + trend

## ----model-grid---------------------------------------------------------------
# Set model options
kardl_set(criterion = "BIC", different_asym_lag = TRUE, data = imf_example_data)
# Estimate model with grid mode
kardl_model <- kardl(
  data = imf_example_data, formula = my_formula,
  maxlag = 4, mode = "grid"
)
# View results
kardl_model

## ----model-grid-custom-summary------------------------------------------------
# Display model summary
summary(kardl_model)

## ----model-user-defined-------------------------------------------------------
kardl_model2 <- kardl(
  data = imf_example_data, my_formula,
  mode = c(2, 1, 1, 3, 0)
)
# View results
kardl_extract(kardl_model2, "opt_lag")

## ----model-user-defined-summary-----------------------------------------------
# Display model summary
summary(kardl_model2)

## ----model-all-vars-----------------------------------------------------------
kardl_set(data = imf_example_data)
kardl(formula = CPI ~ . + deterministic(covid), mode = "grid")

## ----lag-criteria-------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
# Convert lag_criteria to a data frame
lag_criteria <- as.data.frame(kardl_extract(kardl_model, "lag_criteria"))
colnames(lag_criteria) <- c("lag", "AIC", "BIC", "AICc", "HQ")
lag_criteria <- lag_criteria |> mutate(across(c(AIC, BIC, HQ), as.numeric))

# Pivot to long format
lag_criteria_long <- lag_criteria |>
  select(-AICc) |>
  pivot_longer(
    cols = c(AIC, BIC, HQ),
    names_to = "Criteria",
    values_to = "Value"
  )

# Find minimum values
min_values <- lag_criteria_long |>
  group_by(Criteria) |>
  slice_min(order_by = Value) |>
  ungroup()

# Plot
ggplot(
  lag_criteria_long,
  aes(x = lag, y = Value, color = Criteria, group = Criteria)
) +
  geom_line() +
  geom_point(
    data = min_values, aes(x = lag, y = Value),
    color = "red", size = 3, shape = 8
  ) +
  geom_text(
    data = min_values, aes(x = lag, y = Value, label = lag),
    vjust = 1.5, color = "black", size = 3.5
  ) +
  labs(
    title = "Lag Criteria Comparison",
    x = "Lag Configuration",
    y = "Criteria Value"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## ----ecm-estimation-----------------------------------------------------------
ecm_model <- ecm(
  data = imf_example_data, formula = my_formula,
  maxlag = 4, mode = "grid_custom"
)
# View results
summary(ecm_model)

## ----long-run-----------------------------------------------------------------
# Long-run coefficients
my_long <- kardl_longrun(kardl_model)
my_long

## ----long-run-summary---------------------------------------------------------
# Summary of long-run coefficients
summary(my_long)

## ----asymmetry-test-----------------------------------------------------------
ast <- imf_example_data |>
  kardl(
    CPI ~ ER + PPI + asymmetric(ER + PPI) +
      deterministic(covid) + trend,
    mode = c(1, 2, 3, 0, 1),
    data = _
  ) |>
  symmetrytest()
ast

## ----asymmetry-test-summary---------------------------------------------------
# Summary of symmetry test
summary(ast)

## ----pssf---------------------------------------------------------------------
test_result <- kardl_model |> pssf(case = 3, signif_level = "0.05")
test_result

## ----pssf-summary-------------------------------------------------------------
summary(test_result)

## ----psst---------------------------------------------------------------------
test_result <- kardl_model |> psst(case = 3, signif_level = "0.05")
test_result

## ----psst-summary-------------------------------------------------------------
summary(test_result)

## ----narayan------------------------------------------------------------------
test_result <- kardl_model |> narayan(case = 3, signif_level = "0.05")
test_result

## ----narayan-summary----------------------------------------------------------
summary(test_result)

## ----dynamic-multipliers------------------------------------------------------
multipliers <- kardl_model |> mplier()
# View multipliers of the model
head(kardl_extract(multipliers, "multipliers"))
# View long-run multipliers
kardl_extract(multipliers, "omega")
# View short-run multipliers
head(kardl_extract(multipliers, "lambda"))

## ----plot-multipliers---------------------------------------------------------
plot(multipliers, variables = c("ER", "PPI"))

## ----bootstrap-multipliers----------------------------------------------------
bootstrap_results <- kardl_model |>
  bootstrap(horizon = 12, replications = 10)
# View bootstrap summary
summary(bootstrap_results)

## ----plot-bootstrap-multipliers-----------------------------------------------
plot(bootstrap_results, variables = "ER")

## ----asym-custom--------------------------------------------------------------
# Set custom prefixes and suffixes
kardl_reset()
kardl_set(asym_prefix = c("asyP_", "asyN_"), asym_suffix = c("_PP", "_NN"))
kardl_custom <- kardl(data = imf_example_data, my_formula)
kardl_custom

