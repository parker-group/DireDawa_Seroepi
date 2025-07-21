
# Addis Ababa Seroprevalence Analysis by Age Group (Dengue, Zika, Chikungunya)
# Author: Daniel M. Parker
# Date: July 2025
# Description: Generates bar plots and summary tables of IgG/IgM seroprevalence by age group and by virus

library(geepack)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(writexl)

# Load data
dat_addis <- read.csv("data/AddisAbaba_Seroepi24.csv", header = TRUE, sep = ",")

# Assign age groups
dat_addis <- dat_addis %>%
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
  )) %>%
  mutate(AgeGroup = factor(AgeGroup, levels = c("0-4", "5-9", "10-19", "20-29", "30-39", "40-49", "50-59", "60-69", "70+")))

# Function for Wilson score confidence intervals
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

# IgG
dengue_igg <- calc_seroprevalence(dat_addis, "DengG") %>% mutate(Disease = "Dengue")
zika_igg <- calc_seroprevalence(dat_addis, "ZIKG") %>% mutate(Disease = "Zika")
chik_igg <- calc_seroprevalence(dat_addis, "CHIKG") %>% mutate(Disease = "Chikungunya")
igg_combined <- bind_rows(dengue_igg, zika_igg, chik_igg)

# IgM
dengue_igm <- calc_seroprevalence(dat_addis, "DengM") %>% mutate(Disease = "Dengue")
zika_igm <- calc_seroprevalence(dat_addis, "ZIKM") %>% mutate(Disease = "Zika")
chik_igm <- calc_seroprevalence(dat_addis, "CHIKM") %>% mutate(Disease = "Chikungunya")
igm_combined <- bind_rows(dengue_igm, zika_igm, chik_igm)

# Plot IgG
plot_igg_addis <- ggplot(igg_combined, aes(x = AgeGroup, y = Seroprevalence, fill = Disease)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), position = position_dodge(0.8), width = 0.2) +
  labs(y = "Proportion positive (CI)", title = "IgG Seroprevalence (Addis Ababa)", fill=NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_fill_manual(values = c("Dengue" = "skyblue", "Zika" = "lightgreen", "Chikungunya" = "lightcoral")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, size = 11))

# Plot IgM
plot_igm_addis <- ggplot(igm_combined, aes(x = AgeGroup, y = Seroprevalence, fill = Disease)) +
  geom_bar(stat = "identity", position = position_dodge(0.8), color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = CI_Lower, ymax = CI_Upper), position = position_dodge(0.8), width = 0.2) +
  labs(x = "Age Group", y = "Proportion positive (CI)", title = "IgM Seroprevalence (Addis Ababa)", fill=NULL) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent_format()) +
  scale_fill_manual(values = c("Dengue" = "skyblue", "Zika" = "lightgreen", "Chikungunya" = "lightcoral")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11))

# Save as PNG
png("figures/Addis_AgeSeroprev.png", width = 1700, height = 1700, res = 300)
grid.arrange(plot_igg_addis, plot_igm_addis, ncol = 1)
dev.off()

# Output tables
igg_table_addis <- dat_addis %>%
  group_by(AgeGroup) %>%
  summarize(
    Dengue_IgG_Pos = sum(DengG, na.rm = TRUE),
    Zika_IgG_Pos = sum(ZIKG, na.rm = TRUE),
    Chikungunya_IgG_Pos = sum(CHIKG, na.rm = TRUE),
    Total_Tested = n()
  )

igm_table_addis <- dat_addis %>%
  group_by(AgeGroup) %>%
  summarize(
    Dengue_IgM_Pos = sum(DengM, na.rm = TRUE),
    Zika_IgM_Pos = sum(ZIKM, na.rm = TRUE),
    Chikungunya_IgM_Pos = sum(CHIKM, na.rm = TRUE),
    Total_Tested = n()
  )

write_xlsx(
  list(
    IgG_Summary_Addis = igg_table_addis,
    IgM_Summary_Addis = igm_table_addis
  ),
  path = "tables/SuppleTableA.xlsx"
)
