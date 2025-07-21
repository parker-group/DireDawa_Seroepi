
################################################################################
# Mixed Effects Modeling of IgG and IgM Antibody Levels by Virus (DENV, ZIKV, CHIKV)
# 
# Description:
#   - Fits linear mixed models (random intercept for HH_ID) to examine associations 
#     among quantitative IgG and IgM levels across dengue, Zika, and chikungunya.
#   - Exports model estimates, confidence intervals, and p-values to Excel.
#   - Calculates intra-class correlation coefficients (ICCs) to quantify household clustering.
#   - Visualizes pairwise Spearman correlations between antibody types in a heatmap.
#
# Author: DM Parker
# Last updated: 2025-07-21
# Source data: DireDawa_SeroEpi24.csv
################################################################################

##############################################################
#### mixed models for quant IgG and IgM across each other ####
##############################################################


#load data
dat1 <- read.csv("data/DireDawa_SeroEpi24.csv", header=T, sep=",")

# Convert HH_ID to a factor
dat1$HH_ID <- as.factor(dat1$HH_ID)

# Load necessary packages
library(lme4)
library(lmerTest)
library(broom)
library(openxlsx)

# Fit the IgG models
igg_model_dengue <- lmer(RDT.DENIgG ~ RDT.ZIKIgG + RDT.CHKIgG + (1 | HH_ID), data = dat1)
igg_model_zika <- lmer(RDT.ZIKIgG ~ RDT.DENIgG + RDT.CHKIgG + (1 | HH_ID), data = dat1)
igg_model_chik <- lmer(RDT.CHKIgG ~ RDT.DENIgG + RDT.ZIKIgG + (1 | HH_ID), data = dat1)

# Fit the IgM models
igm_model_dengue <- lmer(RDT.DENIgM ~ RDT.ZIKIgM + RDT.CHKIgM + RDT.ZIKIgG + RDT.CHKIgG + (1 | HH_ID), data = dat1)
igm_model_zika <- lmer(RDT.ZIKIgM ~ RDT.DENIgM + RDT.CHKIgM + RDT.DENIgG + RDT.CHKIgG + (1 | HH_ID), data = dat1)
igm_model_chik <- lmer(RDT.CHKIgM ~ RDT.DENIgM + RDT.ZIKIgM + RDT.DENIgG + RDT.ZIKIgG + (1 | HH_ID), data = dat1)


# Extract fixed effects estimates, standard errors, and p-values for IgG models
coef_dengue_igg <- summary(igg_model_dengue)$coefficients
coef_zika_igg <- summary(igg_model_zika)$coefficients
coef_chik_igg <- summary(igg_model_chik)$coefficients

# Extract confidence intervals for IgG models
ci_dengue_igg <- confint(igg_model_dengue, parm = "beta_", method = "Wald")
ci_zika_igg <- confint(igg_model_zika, parm = "beta_", method = "Wald")
ci_chik_igg <- confint(igg_model_chik, parm = "beta_", method = "Wald")

# Combine the IgG results into a single data frame
table_data_igg <- data.frame(
  Model = rep(c("Dengue IgG", "Zika IgG", "Chikungunya IgG"), each = nrow(coef_dengue_igg)),
  Predictor = rep(rownames(coef_dengue_igg), times = 3),
  Estimate = c(coef_dengue_igg[, "Estimate"], coef_zika_igg[, "Estimate"], coef_chik_igg[, "Estimate"]),
  Std.Error = c(coef_dengue_igg[, "Std. Error"], coef_zika_igg[, "Std. Error"], coef_chik_igg[, "Std. Error"]),
  CI.Lower = c(ci_dengue_igg[, 1], ci_zika_igg[, 1], ci_chik_igg[, 1]),
  CI.Upper = c(ci_dengue_igg[, 2], ci_zika_igg[, 2], ci_chik_igg[, 2]),
  p.value = c(coef_dengue_igg[, "Pr(>|t|)"], coef_zika_igg[, "Pr(>|t|)"], coef_chik_igg[, "Pr(>|t|)"])
)

# Format the IgG table
table_data_igg$Estimate <- round(table_data_igg$Estimate, 3)
table_data_igg$Std.Error <- round(table_data_igg$Std.Error, 3)
table_data_igg$CI.Lower <- round(table_data_igg$CI.Lower, 3)
table_data_igg$CI.Upper <- round(table_data_igg$CI.Upper, 3)
table_data_igg$p.value <- formatC(table_data_igg$p.value, format = "e", digits = 2)
table_data_igg$CI <- paste0("(", round(table_data_igg$CI.Lower, 3), ", ", round(table_data_igg$CI.Upper, 3), ")")



# Extract fixed effects estimates, standard errors, and p-values for IgM models
coef_dengue_igm <- summary(igm_model_dengue)$coefficients
coef_zika_igm <- summary(igm_model_zika)$coefficients
coef_chik_igm <- summary(igm_model_chik)$coefficients

# Extract confidence intervals for IgM models
ci_dengue_igm <- confint(igm_model_dengue, parm = "beta_", method = "Wald")
ci_zika_igm <- confint(igm_model_zika, parm = "beta_", method = "Wald")
ci_chik_igm <- confint(igm_model_chik, parm = "beta_", method = "Wald")

# Combine the IgM results into a single data frame
table_data_igm <- data.frame(
  Model = rep(c("Dengue IgM", "Zika IgM", "Chikungunya IgM"), each = nrow(coef_dengue_igm)),
  Predictor = rep(rownames(coef_dengue_igm), times = 3),
  Estimate = c(coef_dengue_igm[, "Estimate"], coef_zika_igm[, "Estimate"], coef_chik_igm[, "Estimate"]),
  Std.Error = c(coef_dengue_igm[, "Std. Error"], coef_zika_igm[, "Std. Error"], coef_chik_igm[, "Std. Error"]),
  CI.Lower = c(ci_dengue_igm[, 1], ci_zika_igm[, 1], ci_chik_igm[, 1]),
  CI.Upper = c(ci_dengue_igm[, 2], ci_zika_igm[, 2], ci_chik_igm[, 2]),
  p.value = c(coef_dengue_igm[, "Pr(>|t|)"], coef_zika_igm[, "Pr(>|t|)"], coef_chik_igm[, "Pr(>|t|)"])
)

# Format the IgM table
table_data_igm$Estimate <- round(table_data_igm$Estimate, 3)
table_data_igm$Std.Error <- round(table_data_igm$Std.Error, 3)
table_data_igm$CI.Lower <- round(table_data_igm$CI.Lower, 3)
table_data_igm$CI.Upper <- round(table_data_igm$CI.Upper, 3)
table_data_igm$p.value <- formatC(table_data_igm$p.value, format = "e", digits = 2)
table_data_igm$CI <- paste0("(", round(table_data_igm$CI.Lower, 3), ", ", round(table_data_igm$CI.Upper, 3), ")")


# Create a workbook object
wb <- createWorkbook()

# Add sheets for IgG and IgM results
addWorksheet(wb, "IgG Results")
addWorksheet(wb, "IgM Results")

# Write the tables to their respective sheets
writeData(wb, "IgG Results", table_data_igg)
writeData(wb, "IgM Results", table_data_igm)

# Save the workbook to the specified directory
## this is Table B and Table C in supplemental materials
saveWorkbook(wb, "tables/Mixed_Model_Results_IgG_IgM.xlsx", overwrite = TRUE)



########################################################
##### correlation matrix ###############################
########################################################

# Load necessary libraries
library(ggplot2)
library(reshape2)
library(dplyr)

# Extract the relevant columns for IgG and IgM values
igg_igm_data <- dat1[, c("RDT.DENIgG", "RDT.ZIKIgG", "RDT.CHKIgG", 
                         "RDT.DENIgM", "RDT.ZIKIgM", "RDT.CHKIgM")]

# Convert all columns to numeric
igg_igm_data <- igg_igm_data %>%
  mutate(across(everything(), as.numeric))
  
  
# Calculate the correlation matrix using Spearman's rho
cor_matrix <- cor(igg_igm_data, use = "complete.obs", method = "spearman")

# Rename the labels for better readability
colnames(cor_matrix) <- rownames(cor_matrix) <- c("DENGV IgG", "ZIKV IgG", "CHIKV IgG", 
                                                  "DENGV IgM", "ZIKV IgM", "CHIKV IgM")

# Melt the correlation matrix for visualization
cor_matrix_melted <- melt(cor_matrix)

# Filter to keep only the lower triangle of the correlation matrix
cor_matrix_melted <- cor_matrix_melted %>%
  filter(as.numeric(Var1) >= as.numeric(Var2))

# Plot the heatmap with points and labels, without axis labels
scatter1 <- ggplot(cor_matrix_melted, aes(Var1, Var2, fill = value, label = round(value, 2))) +
  geom_tile(color = "white") +
  geom_text(color = "black", size = 5) +  # Add correlation values as text
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1, 1), space = "Lab", 
                       name = "Spearman\nCorrelation") +
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),  # Remove x-axis title
        axis.title.y = element_blank()) +  # Remove y-axis title
  coord_fixed() +
  ggtitle("Spearman Correlation Matrix of IgG and IgM Values")



##save the figure as a high resoltion PDF
ggsave("figures/Fig3.pdf", 
       plot = scatter1, 
       width = 6, height = 6, dpi = 300, device = "pdf")


#####################################################################
##Household clustering analysis (looking at ICC)
### note that this depends on the linear mixed models above
#####################################################################

# Extract random effects variances and residual variance
var_dengue_igg <- as.data.frame(VarCorr(igg_model_dengue))
var_zika_igg <- as.data.frame(VarCorr(igg_model_zika))
var_chik_igg <- as.data.frame(VarCorr(igg_model_chik))

var_dengue_igm <- as.data.frame(VarCorr(igm_model_dengue))
var_zika_igm <- as.data.frame(VarCorr(igm_model_zika))
var_chik_igm <- as.data.frame(VarCorr(igm_model_chik))

# Calculate ICC for each model
icc_dengue_igg <- var_dengue_igg[1, "vcov"] / (var_dengue_igg[1, "vcov"] + attr(VarCorr(igg_model_dengue), "sc")^2)
icc_zika_igg <- var_zika_igg[1, "vcov"] / (var_zika_igg[1, "vcov"] + attr(VarCorr(igg_model_zika), "sc")^2)
icc_chik_igg <- var_chik_igg[1, "vcov"] / (var_chik_igg[1, "vcov"] + attr(VarCorr(igg_model_chik), "sc")^2)

icc_dengue_igm <- var_dengue_igm[1, "vcov"] / (var_dengue_igm[1, "vcov"] + attr(VarCorr(igm_model_dengue), "sc")^2)
icc_zika_igm <- var_zika_igm[1, "vcov"] / (var_zika_igm[1, "vcov"] + attr(VarCorr(igm_model_zika), "sc")^2)
icc_chik_igm <- var_chik_igm[1, "vcov"] / (var_chik_igm[1, "vcov"] + attr(VarCorr(igm_model_chik), "sc")^2)

# Display the ICC for each model
icc_dengue_igg
icc_zika_igg
icc_chik_igg

icc_dengue_igm
icc_zika_igm
icc_chik_igm

icc_table <- data.frame(
  Model = c("Dengue IgG", "Zika IgG", "Chikungunya IgG", "Dengue IgM", "Zika IgM", "Chikungunya IgM"),
  ICC = round(c(icc_dengue_igg, icc_zika_igg, icc_chik_igg, icc_dengue_igm, icc_zika_igm, icc_chik_igm), 3)
)