#Creator: Erik Suring (erik.suring@odfw.oregon.gov)
#
#Six ODFW Life Cycle Monitoring sites formed the basis of a Marine Survival Index
#used to forecast Coho Salmon marine survival as part of PFMC Amendment 13
#fisheries management. In 2017 monitoring was cut at the NF Nehalem site. In 2025
#monitoring was cut at the WF Smith and Winchester sites. This code 
#evaluates the marine survival forecasts fit to alternative Marine Survival Indexes
#in MSI_forecast.r.
#

require(tidyverse)
require(gt)
require(scales)

#Coho Salmon marine survival forecast output from MSI_forecast.r
MSI_forecast <- read_csv("Forecast_eval.csv")
#OPIH Jacks/smolt from PMFC Preseason I, Table C-2
OPIH_jacks <- read_csv("OPIH_jacks.csv")

#A13 matrix MS categories from PMFC Preseason I, Table A-3 & A-4
MS_breaks <- c(0, 0.02, 0.045, 0.08, Inf)
OPIH_breaks <- c(0, 0.0008, 0.0014, 0.0040, Inf)
Break_labels <- c("ExLow", "Low", "Medium", "High")

#MSI_table <- MSI_forecast%>%
#  pivot_wider(names_from = MSI_type, values_from = Ensemble_forecast, id_cols=Forecast_year)
#MSI_table$CurrentMSI_cat <- cut(MSI_table$CurrentMSI,
#                                breaks = MS_breaks,
#                                labels = Break_labels,
#                                right = TRUE)

#Categorize the forecast survival
MSI_forecast$MSI_category <- cut(MSI_forecast$Ensemble_forecast,
                                 breaks = MS_breaks,
                                 labels = Break_labels,
                                 right = TRUE)
OPIH_jacks$MSI_category <- cut(OPIH_jacks$OPIH_Jacks_Smolt,
                                breaks = OPIH_breaks,
                                labels = Break_labels,
                                right = TRUE)
OPIH_jacks$MSI_type <- "OPIH_jacks"

#Add OPIH jacks to MSI_forecast table
OPIH_jacks <- OPIH_jacks %>%
  rename_at("OPIH_Jacks_Smolt", ~"Ensemble_forecast") %>%
  rename_at("ReturnYear", ~"Forecast_year")
OPIH_jacks <- OPIH_jacks %>%
  dplyr::select(-c("Jacks_Total_OPIH", "Smolts_Normal", "Delayed", "Smolts_Total")) %>%
  filter(between(Forecast_year, 2021, 2025))
MSI_table <- bind_rows(MSI_forecast, OPIH_jacks)

# Create a summary column for the table
MSI_table <- MSI_table %>%
  mutate(Forecast_Info = paste0("**", MSI_category, "**<br>", round(Ensemble_forecast * 100, 1), "%"))

# Pivot the data
MSI_table_pivot <- MSI_table %>%
  select(Forecast_year, MSI_type, Forecast_Info) %>%
  pivot_wider(names_from = MSI_type, values_from = Forecast_Info)

# Create the results table
Forecast_results_table <- MSI_table_pivot %>%
  gt(rowname_col = "Forecast_year") %>%
  tab_header(
    title = "Marine Survival Forecasts",
  ) %>%
  cols_label(
    Forecast_year = "Forecast Year",
    CurrentMSI = "Five Site MSI",
    InterpolatedMSI = "Interpolated MSI",
    ThreeSiteMSI = "Three Site MSI",
    FixedMSI = "Fixed Model",
    OPIH_jacks = "OPIH Indicator"
  ) %>%
  fmt_markdown(columns = everything()) %>%
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_body(columns = everything())
  ) %>%
  tab_options(
    table.font.size = "12px"
  )

# Color the cells that doesn't match the current MSI
for (col in names(MSI_table_pivot)[-1]) {
  Forecast_results_table <- Forecast_results_table %>%
    tab_style(
      style = cell_fill(color = "red", alpha = 0.2),
      locations = cells_body(
        columns = col,
        rows = grepl("\\*\\*High\\*\\*", MSI_table_pivot[[col]])
      )
    )
}

#Display and save the results table
print(Forecast_results_table)
gtsave(data = Forecast_results_table,
       filename = "Forecast_results.png",
       expand = 5,  # Optional: adds padding around the table
       vwidth = 7 * 96,  # Convert inches to pixels (assuming 96 DPI)
       vheight = 4 * 96
      )


MSI_forecast$MSI_type <- factor(MSI_forecast$MSI_type,
                                levels = c("CurrentMSI", "InterpolatedMSI", "ThreeSiteMSI", "FixedMSI"),
                                labels = c("Five site MSI (2017)", "Interpolated MSI", "Three site MSI", "Fixed model"))

Forecast_results_plot <- ggplot(data = MSI_forecast, aes(x = Forecast_year, y = Ensemble_forecast, color = MSI_type)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Five site MSI (2017)" = "black",
      "Three site MSI" = "blue",
      "Interpolated MSI" = "orange",
      "Fixed model" = "green"
    )
  ) +
  scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
  coord_cartesian(ylim = c(0.055, 0.08)) +
  xlab("Forecast Year") +
  ylab("Forecast Marine Survival") +
  #ggtitle("MSI Alternatives") +
  theme_bw() +
  theme(legend.position = c(0.07, 0.95),
        legend.justification = c(0, 1),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
  )
print(Forecast_results_plot)       
ggsave("Forecast_results_plot.png", plot = Forecast_results_plot +
         theme(plot.margin = margin(10, 30, 10, 10)) +
         coord_cartesian(clip = "off"),
       width = 7, height = 4, units = "in", dpi = 300)


