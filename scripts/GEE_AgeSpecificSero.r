#####################################################################################################
## Script: Age-stratified seroprevalence curves using GEE (IgG and IgM)
## Description: Estimates predicted probability of seropositivity by age, using GEEs with splines.
## Outcome: Side-by-side plots for IgG and IgM, stratified by virus, with household-level clustering
## Data: Dire Dawa serosurvey, cleaned file: DireDawa_SeroEpi24.csv
## Author: DM Parker
## Last updated: July 2025
#####################################################################################################

# Load necessary libraries
library(splines)
library(geepack)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(geepack)



#load data
dat1 <- read.csv("data/DireDawa_SeroEpi24.csv", header=T, sep=",")

####################################################################################################
###  using GEEs to calculate the predicted seropositivity by age accounting for HH clustering  #####
####################################################################################################


# Ensure HH_ID is treated as a factor
dat1$HH_ID <- as.factor(dat1$HH_ID)

# Function to fit GEE model with splines and calculate confidence intervals
fit_gee_model <- function(data, response_var, age_var, id_var, disease_name, df = 3) {
  # Define the formula using the response variable and age variable with splines
  formula <- as.formula(paste(response_var, "~ ns(", age_var, ", df =", df, ")"))
  # Fit the GEE model
  gee_model <- geeglm(formula, family = binomial(link = "logit"), data = data, id = data[[id_var]], corstr = "exchangeable")
  summary(gee_model)
  
  # Define a range of ages for prediction
  age_range <- seq(min(data[[age_var]]), max(data[[age_var]]), length.out = 100)
  
  # Create a new data frame for prediction
  new_data <- data.frame(D1.AGEX = age_range)
  
  # Predict the linear predictor
  linear_predictor <- predict(gee_model, newdata = new_data, type = "link")
  
  # Extract coefficients and variance-covariance matrix
  beta <- coef(gee_model)
  vcov_matrix <- vcov(gee_model)
  
  # Create design matrix for new data
  design_matrix <- model.matrix(~ ns(D1.AGEX, df = 3), data = new_data)
  
  # Calculate standard errors of the linear predictor
  se_fit <- sqrt(rowSums((design_matrix %*% vcov_matrix) * design_matrix))
  
  # Calculate the predicted probabilities using the logistic function
  predicted_prob <- plogis(linear_predictor)
  
  # Calculate 95% confidence intervals
  z_value <- qnorm(0.975)
  ci_lower <- plogis(linear_predictor - z_value * se_fit)
  ci_upper <- plogis(linear_predictor + z_value * se_fit)
  
  # Create a data frame for plotting
  plot_data <- data.frame(Age = age_range, PredictedProb = predicted_prob, CI_Lower = ci_lower, CI_Upper = ci_upper, Disease = disease_name)
  return(plot_data)
}

###############################
## for IgG first
###############################

# Prepare data and fit models for Dengue, Chikungunya, and Zika
# `disease_name` is passed to label each disease
dat1dg <- dat1[!is.na(dat1$DengG), ]
plot_data_dengue <- fit_gee_model(dat1dg, "DengG", "D1.AGEX", "HH_ID", "Dengue")

dat1cg <- dat1[!is.na(dat1$CHIKG), ]
plot_data_chikungunya <- fit_gee_model(dat1cg, "CHIKG", "D1.AGEX", "HH_ID", "Chikungunya")

dat1zg <- dat1[!is.na(dat1$ZIKG), ]
plot_data_zika <- fit_gee_model(dat1zg, "ZIKG", "D1.AGEX", "HH_ID", "Zika")

# Combine the data
combined_plot_data <- bind_rows(plot_data_dengue, plot_data_chikungunya, plot_data_zika)


# Plot the combined data for IgG with legend turned off
plot_igg <- ggplot(combined_plot_data, aes(x = Age, y = PredictedProb, color = Disease, fill = Disease)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Age (years)", y = "Predicted Probability of Being Seropositive") +
  ggtitle("Predicted Probability of IgG Seropositivity by Age") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 80, by = 10)) +
  scale_color_manual(
    values = c("Dengue" = "skyblue", "Zika" = "lightgreen", "Chikungunya" = "lightcoral"),
    labels = c("Dengue" = "DENV", "Zika" = "ZIKV", "Chikungunya" = "CHIKV")
  ) +
  scale_fill_manual(
    values = c("Dengue" = "skyblue", "Zika" = "lightgreen", "Chikungunya" = "lightcoral"),
    labels = c("Dengue" = "DENV", "Zika" = "ZIKV", "Chikungunya" = "CHIKV")
  ) +
  theme(
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.minor.x = element_blank()
  )


###############################
## then for IgM
###############################


# Prepare data and fit models for Dengue, Chikungunya, and Zika
# `disease_name` is passed to label each disease
dat1dm <- dat1[!is.na(dat1$DengM), ]
plot_data_dengueM <- fit_gee_model(dat1dm, "DengM", "D1.AGEX", "HH_ID", "Dengue")

dat1cm <- dat1[!is.na(dat1$CHIKM), ]
plot_data_chikungunyaM <- fit_gee_model(dat1cm, "CHIKM", "D1.AGEX", "HH_ID", "Chikungunya")

dat1zm <- dat1[!is.na(dat1$ZIKM), ]
plot_data_zikaM <- fit_gee_model(dat1zm, "ZIKM", "D1.AGEX", "HH_ID", "Zika")

# Combine the data
combined_plot_data <- bind_rows(plot_data_dengueM, plot_data_chikungunyaM, plot_data_zikaM)

# Plot the combined data for IgM
plot_igm <- ggplot(combined_plot_data, aes(x = Age, y = PredictedProb, color = Disease, fill = Disease)) +
  geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.2, color = NA) +
  geom_line(linewidth = 1) +
  labs(x = "Age (years)", y = "Predicted Probability of Being Seropositive") +
  ggtitle("Predicted Probability of IgM Seropositivity by Age") +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 80, by = 10)) +
  scale_color_manual(
    values = c("Dengue" = "skyblue", "Zika" = "lightgreen", "Chikungunya" = "lightcoral"),
    labels = c("Dengue" = "DENV", "Zika" = "ZIKV", "Chikungunya" = "CHIKV")
  ) +
  scale_fill_manual(
    values = c("Dengue" = "skyblue", "Zika" = "lightgreen", "Chikungunya" = "lightcoral"),
    labels = c("Dengue" = "DENV", "Zika" = "ZIKV", "Chikungunya" = "CHIKV")
  ) +
  theme(
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.minor.x = element_blank()
  )
#  theme(legend.position = "none")  # Turn off the legend for IgG


##############################################################
## plot age specific prevalence for IgG and IgM side by side
##############################################################

# Arrange the two plots side by side
grid.arrange(plot_igg, plot_igm, ncol = 2)

##save the figure as a high resoltion PDF
ggsave("figures/Fig2.pdf", 
       plot = grid.arrange(plot_igg, plot_igm, ncol = 2), 
       width = 12, height = 3.5, dpi = 300, device = "pdf")
	   

