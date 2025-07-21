
# Dire Dawa Seroprevalence Analysis by Age Group (Dengue, Zika, Chikungunya)
# Author: Daniel M. Parker
# Date: July 2025
# Description: Generates bar plots and summary tables of IgG/IgM seroprevalence by age group and by virus

# Load required libraries
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(writexl)

## load data from Dire Dawa
dat1 <- read.csv("data/DireDawa_PNTD/DireDawa_SeroEpi24.csv", header=T, sep=",")




##################
### here combined into a single plot but using spline functions for the predicted probabilities
# Remove missing values for each outcome
dat1dg <- dat1[!is.na(dat1$DengG), ]
dat1dm <- dat1[!is.na(dat1$DengM), ]
dat1cg <- dat1[!is.na(dat1$CHIKG), ]
dat1cm <- dat1[!is.na(dat1$CHIKM), ]
dat1zg <- dat1[!is.na(dat1$ZIKG), ]
dat1zm <- dat1[!is.na(dat1$ZIKM), ]


################################################################################################################
#### bar plots of age-specific seropositivity with CIs calculated using the Wilson score interval method #######
################################################################################################################

# Load necessary libraries
library(dplyr)
library(ggplot2)

# Create age groups up to 69 in 5-year intervals, and a single group for ages 70 and above
dat1 <- dat1 %>%
  mutate(AgeGroup = cut(D1.AGEX, breaks = c(seq(0, 70, by = 5), Inf), labels = c(paste(seq(0, 65, by = 5), seq(4, 69, by = 5), sep = "-"), "70+"), right = FALSE, include.lowest = TRUE))

dat1 <- dat1 %>%
  mutate(AgeGroup = case_when(
    D1.AGEX < 5 ~ "0-4",
    D1.AGEX >= 5 & D1.AGEX < 10 ~ "5-9",
    D1.AGEX >= 10 & D1.AGEX < 20 ~ "10-19",
    D1.AGEX >= 20 & D1.AGEX < 30 ~ "20-29",
    D1.AGEX >= 30 & D1.AGEX < 40 ~ "30-39",
    D1.AGEX >= 40 & D1.AGEX < 50 ~ "40-49",
    D1.AGEX >= 50 & D1.AGEX < 60 ~ "50-59",
    D1.AGEX >= 60 & D1.AGEX < 70 ~ "60-69",
    D1.AGEX >= 70 ~ "70+"
  ))

dat1$AgeGroup <- factor(dat1$AgeGroup, levels = c(
  "0-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70+"
))



# View the first few rows to check the new AgeGroup variable
head(dat1)

# Function to calculate proportions and Wilson score confidence intervals
calc_seroprevalence <- function(data, response_var) {
  data %>%
    group_by(AgeGroup) %>%
    summarize(
      Seropositive = sum(get(response_var), na.rm = TRUE),
      Total = n(),
      Seroprevalence = Seropositive / Total,
      CI_Lower = {
        p_hat <- Seropositive / Total
        z <- qnorm(0.975)
        denominator <- 1 + z^2 / Total
        numerator <- p_hat + z^2 / (2 * Total) - z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * Total)) / Total)
        numerator / denominator
      },
      CI_Upper = {
        p_hat <- Seropositive / Total
        z <- qnorm(0.975)
        denominator <- 1 + z^2 / Total
        numerator <- p_hat + z^2 / (2 * Total) + z * sqrt((p_hat * (1 - p_hat) + z^2 / (4 * Total)) / Total)
        numerator / denominator
      }
    )
}


# Calculate IgG seroprevalence and confidence intervals for Dengue, Zika, and Chikungunya IgG
dengue_seroprev_igg <- calc_seroprevalence(dat1, "DengG")
dengue_seroprev_igg <- dengue_seroprev_igg %>% mutate(Disease = "Dengue")

zika_seroprev_igg <- calc_seroprevalence(dat1, "ZIKG")
zika_seroprev_igg <- zika_seroprev_igg %>% mutate(Disease = "Zika")

chikungunya_seroprev_igg <- calc_seroprevalence(dat1, "CHIKG")
chikungunya_seroprev_igg <- chikungunya_seroprev_igg %>% mutate(Disease = "Chikungunya")

# Combine the data for all IgGs
combined_seroprev_igg <- bind_rows(dengue_seroprev_igg, zika_seroprev_igg, chikungunya_seroprev_igg)


# Calculate IgM seroprevalence and confidence intervals for Dengue, Zika, and Chikungunya IgM
dengue_seroprev_igm <- calc_seroprevalence(dat1, "DengM")
dengue_seroprev_igm <- dengue_seroprev_igm %>% mutate(Disease = "Dengue")

zika_seroprev_igm <- calc_seroprevalence(dat1, "ZIKM")
zika_seroprev_igm <- zika_seroprev_igm %>% mutate(Disease = "Zika")

chikungunya_seroprev_igm <- calc_seroprevalence(dat1, "CHIKM")
chikungunya_seroprev_igm <- chikungunya_seroprev_igm %>% mutate(Disease = "Chikungunya")

# Combine the data for all IgMs
combined_seroprev_igm <- bind_rows(dengue_seroprev_igm, zika_seroprev_igm, chikungunya_seroprev_igm)


# Create the combined bar plot with confidence intervals for IgG
plot_igg <- ggplot(combined_seroprev_igg, aes(x = AgeGroup, y = Seroprevalence, fill = Disease)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), position = position_dodge(width = 0.8), width = 0.2) +
  labs(y = "Proportion positive (CI)", title = "IgG Seroprevalence", fill=NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_fill_manual(values = c("Dengue" = "skyblue", "Zika" = "lightgreen", "Chikungunya" = "lightcoral")) +
  theme(
    axis.title.x = element_blank(),     # Remove main x-axis title
    axis.text.x = element_text(angle = 90, hjust = 1, size = 11)  # Keep and format tick labels
  )

# Create the combined bar plot with confidence intervals for IgM
plot_igm <- ggplot(combined_seroprev_igm, aes(x = AgeGroup, y = Seroprevalence, fill = Disease)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), position = position_dodge(width = 0.8), width = 0.2) +
  labs(x = "Age Group", y = "Proportion positive (CI)", title = "IgM Seroprevalence", fill=NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_fill_manual(values = c("Dengue" = "skyblue", "Zika" = "lightgreen", "Chikungunya" = "lightcoral"))+
  theme(
#    axis.title.x = element_blank(),     # Remove main x-axis title
    axis.text.x = element_text(angle = 90, hjust = 1, size = 11)  # Keep and format tick labels
  )
  
# Combine the two plots into one figure

# Open PNG device
png("figures/Fig1.png", width = 1700, height = 1700, res = 300)

# Draw the combined plots
grid.arrange(plot_igg, plot_igm, ncol = 1)

# Close the device
dev.off()


##save the figure as a high resoltion PDF
ggsave("figures/Fig1.pdf", 
       plot = grid.arrange(plot_igg, plot_igm, ncol = 1), 
       width = 6, height = 6, dpi = 300, device = "pdf")


#############################################################################################
#####  Generating a table of people tested by age group and text positives by age group #####
#############################################################################################
# Summarize number of IgG positives by AgeGroup for each disease
igg_table <- dat1 %>%
  group_by(AgeGroup) %>%
  summarize(
    Dengue_IgG_Pos = sum(DengG, na.rm = TRUE),
    Zika_IgG_Pos = sum(ZIKG, na.rm = TRUE),
    Chikungunya_IgG_Pos = sum(CHIKG, na.rm = TRUE),
    Total_Tested = n()
  )

# View table
print(igg_table)

igm_table <- dat1 %>%
  group_by(AgeGroup) %>%
  summarize(
    Dengue_IgM_Pos = sum(DengM, na.rm = TRUE),
    Zika_IgM_Pos = sum(ZIKM, na.rm = TRUE),
    Chikungunya_IgM_Pos = sum(CHIKM, na.rm = TRUE),
    Total_Tested = n()
  )

print(igm_table)


#export an excel file with these tables
write_xlsx(
  list(
    IgG_Summary = igg_table,
    IgM_Summary = igm_table
  ),
  path = "tables/Table1.xlsx"
)


