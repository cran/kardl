## ----setup, include=FALSE-----------------------------------------------------

knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)


## ----install-cran, eval=FALSE-------------------------------------------------
# 
# install.packages("kardl")
# library(kardl)
# 

## ----install, eval=FALSE------------------------------------------------------
# 
# # Install required packages
# install.packages(c("stats", "msm", "lmtest", "nlWaldTest", "car", "strucchange", "utils"))
# # Install kardl from GitHub
# install.packages("devtools")
# devtools::install_github("karamelikli/kardl")
# 

## ----load---------------------------------------------------------------------

library(kardl)


## ----data-prepare-------------------------------------------------------------

# Define the model formula
MyFormula <- CPI ~ ER + PPI + asym(ER + PPI) + deterministic(covid) + trend


## ----model-grid---------------------------------------------------------------

# Set model options
kardl_set(criterion = "BIC", differentAsymLag = TRUE, data=imf_example_data)
# Estimate model with grid mode
kardl_model <- kardl(data=imf_example_data,model= MyFormula, maxlag = 4, mode = "grid")
# View results
kardl_model
# Display model summary
summary(kardl_model)


## ----model-user-defined-------------------------------------------------------

kardl_model2 <- kardl(data=imf_example_data, MyFormula, mode = c(2, 1, 1, 3, 0))
# View results
kardl_model2$properLag


## ----model-all-vars-----------------------------------------------------------

kardl_set(data=imf_example_data)
kardl(model =  CPI ~ . + deterministic(covid), mode = "grid")


## ----lag-criteria-------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
# Convert LagCriteria to a data frame
LagCriteria <- as.data.frame(kardl_model[["LagCriteria"]])
colnames(LagCriteria) <- c("lag", "AIC", "BIC", "AICc", "HQ")
LagCriteria <- LagCriteria %>% mutate(across(c(AIC, BIC, HQ), as.numeric))

# Pivot to long format
LagCriteria_long <- LagCriteria %>%
  select(-AICc) %>%
  pivot_longer(cols = c(AIC, BIC, HQ), names_to = "Criteria", values_to = "Value")

# Find minimum values
min_values <- LagCriteria_long %>%
  group_by(Criteria) %>%
  slice_min(order_by = Value) %>%
  ungroup()

# Plot
ggplot(LagCriteria_long, aes(x = lag, y = Value, color = Criteria, group = Criteria)) +
  geom_line() +
  geom_point(data = min_values, aes(x = lag, y = Value), color = "red", size = 3, shape = 8) +
  geom_text(data = min_values, aes(x = lag, y = Value, label = lag), vjust = 1.5, color = "black", size = 3.5) +
  labs(title = "Lag Criteria Comparison", x = "Lag Configuration", y = "Criteria Value") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


## ----long-run-----------------------------------------------------------------

# Long-run coefficients
mylong <- kardl_longrun(kardl_model)
mylong


## ----asymmetry-test-----------------------------------------------------------

ast <- imf_example_data %>% kardl(CPI ~ ER + PPI + asym(ER + PPI) + deterministic(covid) + trend, mode = c(1, 2, 3, 0, 1)) %>% asymmetrytest()
ast
# View long-run hypotheses
ast$Lhypotheses
# Detailed results
summary(ast)


## ----pssf---------------------------------------------------------------------

A <- kardl_model %>% pssf(case = 3, signif_level = "0.05")
cat(paste0("The F statistic = ", A$statistic, " where k = ", A$k, "."))
cat(paste0("\nWe found '", A$Cont, "' at ", A$siglvl, "."))
A$criticalValues
summary(A)


## ----psst---------------------------------------------------------------------

A <- kardl_model %>% psst(case = 3, signif_level = "0.05")
cat(paste0("The t statistic = ", A$statistic, " where k = ", A$k, "."))
cat(paste0("\nWe found '", A$Cont, "' at ", A$siglvl, "."))
A$criticalValues
summary(A)


## ----narayan------------------------------------------------------------------

A <- kardl_model %>% narayan(case = 3, signif_level = "0.05")
cat(paste0("The F statistic = ", A$statistic, " where k = ", A$k, "."))
cat(paste0("\nWe found '", A$Cont, "' at ", A$siglvl, "."))
A$criticalValues
summary(A)


## ----banerjee-----------------------------------------------------------------

A <- kardl_model %>% banerjee(signif_level = "0.05")
cat(paste0("The ECM parameter = ", A$coef, ", k = ", A$k, ", t statistic = ", A$statistic, "."))
cat(paste0("\nWe found '", A$Cont, "' at ", A$siglvl, "."))
A$criticalValues
summary(A)


## ----recmt--------------------------------------------------------------------

recmt_model <- imf_example_data %>% recmt(MyFormula, mode = "grid_custom", case = 3)
recmt_model
# View results
summary(recmt_model)


## ----arch-test----------------------------------------------------------------

arch_result <- archtest(kardl_model$finalModel$model$residuals, q = 2)
summary(arch_result)


## ----asym-custom--------------------------------------------------------------

# Set custom prefixes and suffixes
kardl_reset()
kardl_set(AsymPrefix = c("asyP_", "asyN_"), AsymSuffix = c("_PP", "_NN"))
kardl_custom <- kardl(data=imf_example_data, MyFormula, mode = "grid_custom")
kardl_custom$properLag


