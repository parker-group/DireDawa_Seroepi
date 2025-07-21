
################################################################################
# Title: GEE-Based Seroprevalence Estimation for Arboviruses in Dire Dawa
#
# Description:
# This script estimates overall IgG and IgM seroprevalence for dengue, 
# chikungunya, and Zika using Generalized Estimating Equations (GEE) to 
# account for household-level clustering. Seroprevalence estimates are 
# obtained as the intercept of null GEE models and converted from the logit 
# scale to probabilities. Bootstrapping (n = 1000) is used to compute empirical 
# 95% confidence intervals for each seroprevalence estimate, using individual-level resampling. 
# Household-level resampling is optionally included as a sensitivity check.
# Results are visualized using forest plots and summarized in tabular format.
#
# Input:
# - A .csv file with binary IgG and IgM outcome variables for each virus
# - Household ID variable (HH_ID) for clustering
#
# Output:
# - Forest plots of overall seroprevalence with 95% CIs
# - A results table with formatted point estimates and confidence intervals
#
# Author: D.M. Parker
# Last Updated: July 2025
################################################################################


#load data
dat1 <- read.csv("data/DireDawa_SeroEpi24.csv", header=T, sep=",")


##########################################################################################
###  using GEEs to calculate the seropositivity accounting for HH clustering  ####
##########################################################################################

library(geepack)
library(boot)
library(ggplot2)
library(patchwork)
library(knitr)


# Ensure HH_ID is treated as a factor
dat1$HH_ID <- as.factor(dat1$HH_ID)


#first calculating an overall seroprevalence to IgG and IgM for all 3 viruses
#using GEEs to account for clustering of the data
# using bootstrapping to empirically estimate the confidence intervals around the 
## seroprevalence estimates (which are the intercept in this case, from null models)


# Remove missing values for each outcome
dat1dg <- dat1[!is.na(dat1$DengG), ]
dat1dm <- dat1[!is.na(dat1$DengM), ]
dat1cg <- dat1[!is.na(dat1$CHIKG), ]
dat1cm <- dat1[!is.na(dat1$CHIKM), ]
dat1zg <- dat1[!is.na(dat1$ZIKG), ]
dat1zm <- dat1[!is.na(dat1$ZIKM), ]

# Function to fit the GEE model and return the intercept estimate for Dengue IgG
gee_intercept_dengue_igg <- function(data, indices) {
  d <- data[indices, ]
  gee_model <- geeglm(DengG ~ 1, family = binomial(link = "logit"), data = d, id = HH_ID, corstr = "exchangeable")
  return(coef(gee_model)["(Intercept)"])
}

# Function to fit the GEE model and return the intercept estimate for Dengue IgM
gee_intercept_dengue_igm <- function(data, indices) {
  d <- data[indices, ]
  gee_model <- geeglm(DengM ~ 1, family = binomial(link = "logit"), data = d, id = HH_ID, corstr = "exchangeable")
  return(coef(gee_model)["(Intercept)"])
}

# Function to fit the GEE model and return the intercept estimate for Chikungunya IgG
gee_intercept_chik_igg <- function(data, indices) {
  d <- data[indices, ]
  gee_model <- geeglm(CHIKG ~ 1, family = binomial(link = "logit"), data = d, id = HH_ID, corstr = "exchangeable")
  return(coef(gee_model)["(Intercept)"])
}

# Function to fit the GEE model and return the intercept estimate for Chikungunya IgM
gee_intercept_chik_igm <- function(data, indices) {
  d <- data[indices, ]
  gee_model <- geeglm(CHIKM ~ 1, family = binomial(link = "logit"), data = d, id = HH_ID, corstr = "exchangeable")
  return(coef(gee_model)["(Intercept)"])
}

# Function to fit the GEE model and return the intercept estimate for Zika IgG
gee_intercept_zika_igg <- function(data, indices) {
  d <- data[indices, ]
  gee_model <- geeglm(ZIKG ~ 1, family = binomial(link = "logit"), data = d, id = HH_ID, corstr = "exchangeable")
  return(coef(gee_model)["(Intercept)"])
}

# Function to fit the GEE model and return the intercept estimate for Zika IgM
gee_intercept_zika_igm <- function(data, indices) {
  d <- data[indices, ]
  gee_model <- geeglm(ZIKM ~ 1, family = binomial(link = "logit"), data = d, id = HH_ID, corstr = "exchangeable")
  return(coef(gee_model)["(Intercept)"])
}

# set seed for reproducible boostrap results
set.seed(123)  

## run bootstrapping w/ 1000 resamples to estimate the intercept (logit-scale prevalence) and its variability
# Perform bootstrapping for Dengue IgG
boot_results_dengue_igg <- boot(data = dat1dg, statistic = gee_intercept_dengue_igg, R = 1000)
print(boot_results_dengue_igg)

# Perform bootstrapping for Dengue IgM
boot_results_dengue_igm <- boot(data = dat1dm, statistic = gee_intercept_dengue_igm, R = 1000)
print(boot_results_dengue_igm)

# Perform bootstrapping for Chikungunya IgG
boot_results_chik_igg <- boot(data = dat1cg, statistic = gee_intercept_chik_igg, R = 1000)
print(boot_results_chik_igg)

# Perform bootstrapping for Chikungunya IgM
boot_results_chik_igm <- boot(data = dat1cm, statistic = gee_intercept_chik_igm, R = 1000)
print(boot_results_chik_igm)

# Perform bootstrapping for Zika IgG
boot_results_zika_igg <- boot(data = dat1zg, statistic = gee_intercept_zika_igg, R = 1000)
print(boot_results_zika_igg)

# Perform bootstrapping for Zika IgM
boot_results_zika_igm <- boot(data = dat1zm, statistic = gee_intercept_zika_igm, R = 1000)
print(boot_results_zika_igm)







### Convert bootstrap percentiles from logit scale to probability scale
calculate_ci <- function(boot_results) {
  boot_ci <- boot.ci(boot_results, type = "perc")
  logit_ci <- boot_ci$percent[4:5]
  ci_lower_logit <- logit_ci[1]
  ci_upper_logit <- logit_ci[2]
  ci_lower_prob <- plogis(ci_lower_logit)
  ci_upper_prob <- plogis(ci_upper_logit)
  return(c(ci_lower_prob, ci_upper_prob))
}

# Helper function to summarize prevalence and CIs to console
print_seroprevalence <- function(boot_results, name) {
  # Extract the original intercept estimate
  intercept_estimate <- boot_results$t0
  
  # Convert the log-odds to probability
  seroprevalence <- plogis(intercept_estimate)
  
  # Calculate confidence intervals
  ci <- calculate_ci(boot_results)
  
  # Print the seroprevalence estimate and confidence intervals
  cat(paste(name, "Seroprevalence Estimate:", seroprevalence, "\n"))
  cat(paste(name, "95% Confidence Interval:", ci[1], "-", ci[2], "\n"))
}



########################################################################################################################
########################################################################################################################

# Optional: Use household-level bootstrap for CI sensitivity analysis
## Note: Point estimates remain the same as original GEE. This only affects CI.
household_bootstrap <- function(data, outcome_var, id_var = "HH_ID") {
  hh_ids <- unique(data[[id_var]])
  boot_estimates <- numeric(1000)

  for (i in 1:1000) {
    sampled_hhs <- sample(hh_ids, length(hh_ids), replace = TRUE)
    boot_data <- do.call(rbind, lapply(sampled_hhs, function(hh) data[data[[id_var]] == hh, ]))

    # Fit null GEE model
    model <- tryCatch({
      geeglm(reformulate("1", outcome_var), family = binomial, data = boot_data,
             id = boot_data[[id_var]], corstr = "exchangeable")
    }, error = function(e) NULL)

    # Store estimate if successful
    if (!is.null(model)) {
      boot_estimates[i] <- coef(model)[["(Intercept)"]]
    } else {
      boot_estimates[i] <- NA
    }
  }

  # Drop failed iterations
  boot_estimates <- na.omit(boot_estimates)

  # Convert to probability scale
  ci <- plogis(quantile(boot_estimates, c(0.025, 0.975)))
  return(ci)
}

########################################################################################################################
########################################################################################################################



# Calculate and print seroprevalence estimates and confidence intervals for Dengue IgG
print_seroprevalence(boot_results_dengue_igg, "Dengue IgG")

# Calculate and print seroprevalence estimates and confidence intervals for Dengue IgM
print_seroprevalence(boot_results_dengue_igm, "Dengue IgM")

# Calculate and print seroprevalence estimates and confidence intervals for Chikungunya IgG
print_seroprevalence(boot_results_chik_igg, "Chikungunya IgG")

# Calculate and print seroprevalence estimates and confidence intervals for Chikungunya IgM
print_seroprevalence(boot_results_chik_igm, "Chikungunya IgM")

# Calculate and print seroprevalence estimates and confidence intervals for Zika IgG
print_seroprevalence(boot_results_zika_igg, "Zika IgG")

# Calculate and print seroprevalence estimates and confidence intervals for Zika IgM
print_seroprevalence(boot_results_zika_igm, "Zika IgM")




##########################################################################################
###  Create forest plots and tabular output of seroprevalence results with CIs  #########
##########################################################################################


# Assemble final seroprev estimates and CIs into a dataframe for plotting and table output
results <- data.frame(
  Disease = c("Dengue", "Dengue", "Chikungunya", "Chikungunya", "Zika", "Zika"),
  Antibody = c("IgG", "IgM", "IgG", "IgM", "IgG", "IgM"),
  Estimate = c(plogis(boot_results_dengue_igg$t0), 
               plogis(boot_results_dengue_igm$t0),
               plogis(boot_results_chik_igg$t0),
               plogis(boot_results_chik_igm$t0),
               plogis(boot_results_zika_igg$t0),
               plogis(boot_results_zika_igm$t0)),
  CI_Lower = c(calculate_ci(boot_results_dengue_igg)[1],
               calculate_ci(boot_results_dengue_igm)[1],
               calculate_ci(boot_results_chik_igg)[1],
               calculate_ci(boot_results_chik_igm)[1],
               calculate_ci(boot_results_zika_igg)[1],
               calculate_ci(boot_results_zika_igm)[1]),
  CI_Upper = c(calculate_ci(boot_results_dengue_igg)[2],
               calculate_ci(boot_results_dengue_igm)[2],
               calculate_ci(boot_results_chik_igg)[2],
               calculate_ci(boot_results_chik_igm)[2],
               calculate_ci(boot_results_zika_igg)[2],
               calculate_ci(boot_results_zika_igm)[2])
)

# Reorder the data frame for plotting
results$Disease <- factor(results$Disease, levels = c("Dengue", "Chikungunya", "Zika"))
results$Antibody <- factor(results$Antibody, levels = c("IgG", "IgM"))

# Create a forest plot for IgG
results_igg <- subset(results, Antibody == "IgG")
forest_plot_igg <- ggplot(results_igg, aes(x = Disease, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "IgG (95% CIs)",
       x = "Disease",
       y = "Seroprevalence") +
  ylim(0, 1)

# Create a forest plot for IgM
results_igm <- subset(results, Antibody == "IgM")
forest_plot_igm <- ggplot(results_igm, aes(x = Disease, y = Estimate)) +
  geom_point() +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), width = 0.2) +
  coord_flip() +
  theme_minimal() +
  labs(title = "IgM (95% CIs)",
       x = "",
       y = "Seroprevalence") +
  ylim(0, 1)

# Combine the plots side-by-side
combined_plot <- forest_plot_igg + forest_plot_igm + plot_layout(ncol = 2)

# Print the combined plot
print(combined_plot)

###  Create a simple table of the results as well
# Format the estimates and confidence intervals
results$Estimate <- round(results$Estimate, 3)
results$CI_Lower <- round(results$CI_Lower, 3)
results$CI_Upper <- round(results$CI_Upper, 3)


# Create a formatted summary table (this is Table 2) with the results
kable(results, col.names = c("Disease", "Antibody", "Estimate", "CI Lower", "CI Upper"), 
      caption = "Seroprevalence Estimates with 95% Confidence Intervals")


#############################################################################################################
#######   optional: compare to estimates w/ bootstrapping done at house level ###############################
#############################################################################################################

# calculating CIs from bootstrapping at house level 
ci_dengue_igg_hh <- household_bootstrap(dat1dg, "DengG")
cat("Dengue IgG (HH bootstrap) CI:", round(ci_dengue_igg_hh[1], 3), "-", round(ci_dengue_igg_hh[2], 3), "\n")

ci_dengue_igm_hh <- household_bootstrap(dat1dm, "DengM")
cat("Dengue IgM (HH bootstrap) CI:", round(ci_dengue_igm_hh[1], 3), "-", round(ci_dengue_igm_hh[2], 3), "\n")

ci_chik_igg_hh <- household_bootstrap(dat1cg, "CHIKG")
cat("Chikungunya IgG (HH bootstrap) CI:", round(ci_chik_igg_hh[1], 3), "-", round(ci_chik_igg_hh[2], 3), "\n")

ci_chik_igm_hh <- household_bootstrap(dat1cm, "CHIKM")
cat("Chikungunya IgM (HH bootstrap) CI:", round(ci_chik_igm_hh[1], 3), "-", round(ci_chik_igm_hh[2], 3), "\n")

ci_zika_igg_hh <- household_bootstrap(dat1zg, "ZIKG")
cat("Zika IgG (HH bootstrap) CI:", round(ci_zika_igg_hh[1], 3), "-", round(ci_zika_igg_hh[2], 3), "\n")

ci_zika_igm_hh <- household_bootstrap(dat1zm, "ZIKM")
cat("Zika IgM (HH bootstrap) CI:", round(ci_zika_igm_hh[1], 3), "-", round(ci_zika_igm_hh[2], 3), "\n")





# Compile into a table for side-by-side comparison
hh_bootstrap_results <- data.frame(
  Disease = c("Dengue", "Dengue", "Chikungunya", "Chikungunya", "Zika", "Zika"),
  Antibody = c("IgG", "IgM", "IgG", "IgM", "IgG", "IgM"),
  Estimate = results$Estimate,  # Original point estimates from GEE
  CI_Lower_HH = round(c(ci_dengue_igg_hh[1], ci_dengue_igm_hh[1],
                        ci_chik_igg_hh[1], ci_chik_igm_hh[1],
                        ci_zika_igg_hh[1], ci_zika_igm_hh[1]), 3),
  CI_Upper_HH = round(c(ci_dengue_igg_hh[2], ci_dengue_igm_hh[2],
                        ci_chik_igg_hh[2], ci_chik_igm_hh[2],
                        ci_zika_igg_hh[2], ci_zika_igm_hh[2]), 3)
)

# Create the formatted table
kable(hh_bootstrap_results, col.names = c("Disease", "Antibody", "Estimate", "CI Lower (HH)", "CI Upper (HH)"),
      caption = "Seroprevalence Estimates with Household-Level Bootstrap 95% Confidence Intervals")