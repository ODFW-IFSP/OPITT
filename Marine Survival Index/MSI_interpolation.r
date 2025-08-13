#Creator: Erik Suring (erik.suring@odfw.oregon.gov)
#
#Six ODFW Life Cycle Monitoring sites formed the basis of a Marine Survival Index
#used to forecast Coho Salmon marine survival as part of PFMC Amendment 13
#fisheries management. In 2017 monitoring was cut at the NF Nehalem site. In 2025
#monitoring was cut at the WF Smith and Winchester sites. This code uses the
#remaining LCM sites to interpolate the missing marine survival estimates.
#
#MARSS parameters forked from Washington Department of Fish & Wildlife, Fish Program
#https://github.com/wdfw-fp/OPI-H-Forecast-Evaluation-2023/tree/2024_forecast
#

# Load libraries
require(tidyverse)
require(MARSS)

#Coho Salmon marine survival estimate data from the ODFW LCM project
MSI_data<-read_csv("MSI.csv")
#Convert ratios to logit
MSI_data$MSAdj = log(MSI_data$MSAdj/(1-MSI_data$MSAdj))

marss_dat<-MSI_data%>%
  pivot_wider(names_from = Site,values_from = MSAdj,id_cols=BroodYear)

marss_mat<-marss_dat%>%
  dplyr::select(-c(BroodYear))%>%
  as.matrix()%>%
  t()

model=list(
  Q= "equalvarcov",#"unconstrained",
  R= diag(rep(0,nrow(marss_mat))),#"diagonal and equal","diagonal and unequal",
  U= matrix(rep(0,nrow(marss_mat)),nrow=nrow(marss_mat),1)
)
fit=MARSS(marss_mat, model=model,control=list(maxit=200,allow.degen=T))

fitted<-t(fit$states)
colnames(fitted)<-gsub("X.","mle_",colnames(fitted))

MSI_fitted<-marss_dat%>%
  bind_cols(fitted)%>%
  pivot_longer(cols=c(everything(), -BroodYear),names_to = "Site")%>%
  mutate(type=ifelse(grepl("mle_",Site),"mle","obs"))%>%
  mutate(Site=gsub("mle_","",Site))%>%
  filter(type=="mle")%>%
  group_by(BroodYear,Site)%>%
  summarise(MSAdj=mean(value))

#Convert logit back to ratios
MSI_fitted$MSAdj <- (exp(MSI_fitted$MSAdj)/(exp(MSI_fitted$MSAdj)+1))

write_csv(MSI_fitted, "MSI_fitted.csv")
