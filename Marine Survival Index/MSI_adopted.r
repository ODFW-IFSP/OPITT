#Creator: Erik Suring (erik.suring@odfw.oregon.gov)
#
#Six ODFW Life Cycle Monitoring sites formed the basis of a Marine Survival Index
#used to forecast Coho Salmon marine survival as part of PFMC Amendment 13
#fisheries management. In 2017 monitoring was cut at the NF Nehalem site. In 2025
#monitoring was cut at the WF Smith and Winchester sites. This code uses the
#remaining LCM sites to interpolate the missing marine survival estimates.
#
#MARSS inspiration forked from Washington Department of Fish & Wildlife, Fish Program
#https://github.com/wdfw-fp/OPI-H-Forecast-Evaluation-2023/tree/2024_forecast
#

# Load libraries
library(tidyverse)
library(MARSS)

#Coho Salmon marine survival estimate data from the ODFW LCM project
MSI_data<-read_csv("MSI.csv")
#Convert ratios to logit
MSI_data$MSAdj = log(MSI_data$MSAdj/(1-MSI_data$MSAdj))

MSI_table<-MSI_data%>%
  arrange(ReturnYear)%>%
  pivot_wider(names_from = Site,values_from = MSAdj,id_cols=ReturnYear)

#Original method: average of six LCM sites
MSI_Original <- MSI_data%>%
  group_by(ReturnYear)%>%
  filter(between(ReturnYear, 1999, 2017))%>%
  summarise(MSAdjOG=mean(MSAdj))
MSI_Original$MSAdjOG <- (exp(MSI_Original$MSAdjOG)/(exp(MSI_Original$MSAdjOG)+1))
MSI_ESU <- MSI_Original

#Adopted alternative 2025: use MARSS to interpolate missing values
marss_mat<-MSI_table%>%
  dplyr::select(-c(ReturnYear))%>%
  as.matrix()%>%
  t()

model=list(Q= "equalvarcov")
fit=MARSS(marss_mat, model=model)

fitted<-t(fit$states)
colnames(fitted)<-gsub("X.","mle_",colnames(fitted))

MSI_fitted<-MSI_table%>%
  bind_cols(fitted)%>%
  pivot_longer(cols=c(everything(), -ReturnYear),names_to = "Site")%>%
  mutate(type=ifelse(grepl("mle_",Site),"mle","obs"))%>%
  mutate(Site=gsub("mle_","",Site))%>%
  pivot_wider(names_from = type,values_from = value)%>%
  mutate(MSAdj=ifelse(!is.na(obs),obs,mle))

MSI_fitted <- MSI_fitted%>%
  group_by(ReturnYear)%>%
  summarise(MSAdjMARSS=mean(MSAdj))
MSI_fitted$MSAdjMARSS <- (exp(MSI_fitted$MSAdjMARSS)/(exp(MSI_fitted$MSAdjMARSS)+1))

MSI_ESU <- right_join(MSI_ESU, MSI_fitted)

#Alternative Check: Use a weighted average of the three remaining sites
MSI_Three_Weight <- MSI_table%>%
  dplyr::select(-c("NF Nehalem Total", "West Fork Smith", "Winchester Creek"))%>%
  mutate(Cascade=Cascade*(2.5/6))%>%
  mutate(SiletzMill=SiletzMill*(2.5/6))%>%
  mutate(YaquinaMill=YaquinaMill*(1/6))
#SiletzMill needs more weight in 1999 and 2001 as Cascade is NA
MSI_Three_Weight$SiletzMill[3] <- MSI_Three_Weight$SiletzMill[3]*2
MSI_Three_Weight$SiletzMill[1] <- MSI_Three_Weight$SiletzMill[1]*2

MSI_Three_Weight <- MSI_Three_Weight%>%
  mutate(MSAdjThree = rowSums(MSI_Three_Weight[ , c(2,3,4)], na.rm=TRUE))
MSI_Three <- MSI_Three_Weight%>%
  dplyr::select(-c("Cascade", "SiletzMill", "YaquinaMill"))
MSI_Three$MSAdjThree <- (exp(MSI_Three$MSAdjThree)/(exp(MSI_Three$MSAdjThree)+1))
MSI_ESU <- left_join(MSI_ESU, MSI_Three)

#Print graph of MSI input to be used in the marine survival forecast
MSI_graph <- MSI_ESU
MSI_graph_long <- MSI_graph %>%
  select(ReturnYear, MSAdjMARSS, MSAdjThree, MSAdjOG)
MSI_graph_long <- pivot_longer(MSI_graph_long,
                              cols = c(MSAdjMARSS, MSAdjThree, MSAdjOG),
                              names_to = "Metric",
                              values_to = "Value")
MSI_graph_long$Metric <- factor(MSI_graph_long$Metric,
                                levels = c("MSAdjOG", "MSAdjThree", "MSAdjMARSS"),
                                labels = c("Six site MSI (2013)", "Three site MSI", "Interpolated MSI (2025)"))

MSI_alt_plot <- ggplot(data = MSI_graph_long, aes(x = ReturnYear, y = Value, color = Metric)) +
  geom_line(linewidth = 1.5) +
  geom_point(size = 2) +
  scale_color_manual(
    name = NULL,
    values = c(
      "Six site MSI (2013)" = "#000000",
      "Three site MSI" = "#0072B2",
      "Interpolated MSI (2025)" = "#D55E00"
    )
  ) +
  scale_y_continuous(
    breaks = c(.02, 0.045, 0.08),
    labels = scales::label_percent(accuracy = 0.1)(c(.02, 0.045, 0.08))
    ) +
  xlab("Return Year") +
  ylab("Adjusted Marine Survival") +
  #ggtitle("MSI Alternatives") +
  theme_bw() +
  theme(legend.position = c(0.05, 0.95),
        legend.justification = c(0, 1),
        axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12)
  )
print(MSI_alt_plot)       
ggsave("MSI_plot.png", plot = MSI_alt_plot +
         theme(plot.margin = margin(10, 30, 10, 10)) +
         coord_cartesian(clip = "off"),
       width = 7, height = 4, units = "in", dpi = 300)

#Create the MSI input to be used in the marine survival forecast
MSI_results <- MSI_ESU %>%
  select(ReturnYear, MSAdjMARSS)

write_csv(MSI_results, "MSI_results.csv")  
