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
require(Metrics)
require(gt)

#Coho Salmon marine survival estimate data from the ODFW LCM project
MSI_data<-read_csv("MSI.csv")
#Convert ratios to logit
MSI_data$MSAdj = log(MSI_data$MSAdj/(1-MSI_data$MSAdj))

MSI_table<-MSI_data%>%
  pivot_wider(names_from = Site,values_from = MSAdj,id_cols=ReturnYear)

#Original method: average of six LCM sites
MSI_Original <- MSI_data%>%
  group_by(ReturnYear)%>%
  filter(between(ReturnYear, 1999, 2017))%>%
  summarise(MSAdjOG=mean(MSAdj))
MSI_Original$MSAdjOG <- (exp(MSI_Original$MSAdjOG)/(exp(MSI_Original$MSAdjOG)+1))
MSI_ESU <- MSI_Original

#Current method: weighted average of the remaining five sites
MSI_Current_Weight <- MSI_data%>%
  group_by(ReturnYear)%>%
  filter(Site != "YaquinaMill")%>%
  summarise(MSAdj=mean(MSAdj))%>%
  mutate(Site = "PostNehalemWeight")

MSI_Current <- bind_rows(MSI_Current_Weight, MSI_data)
MSI_Current <- MSI_Current%>%
  group_by(ReturnYear)%>%
  filter(Site != "NF Nehalem Total")%>%
  summarise(MSAdjCurrent=mean(MSAdj))
MSI_Current$MSAdjCurrent <- (exp(MSI_Current$MSAdjCurrent)/(exp(MSI_Current$MSAdjCurrent)+1))
MSI_ESU <- left_join(MSI_Current, MSI_ESU)

#Alternative One: use MARSS to interpolate missing values
marss_mat<-MSI_table%>%
  dplyr::select(-c(ReturnYear))%>%
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

MSI_fitted<-MSI_table%>%
  bind_cols(fitted)%>%
  pivot_longer(cols=c(everything(), -ReturnYear),names_to = "Site")%>%
  mutate(type=ifelse(grepl("mle_",Site),"mle","obs"))%>%
  mutate(Site=gsub("mle_","",Site))%>%
  filter(type=="mle")%>%
  group_by(ReturnYear,Site)%>%
  summarise(MSAdj=mean(value))

MSI_fitted <- MSI_fitted%>%
  group_by(ReturnYear)%>%
  summarise(MSAdjMARSS=mean(MSAdj))
MSI_fitted$MSAdjMARSS <- (exp(MSI_fitted$MSAdjMARSS)/(exp(MSI_fitted$MSAdjMARSS)+1))

MSI_ESU <- left_join(MSI_ESU, MSI_fitted)
  
#Alternative Two: Use a weighted average of the three remaining sites
MSI_Three_Weight <- MSI_table%>%
  dplyr::select(-c("NF Nehalem Total", "West Fork Smith", "Winchester Creek"))%>%
  mutate(Cascade=Cascade*(2.5/6))%>%
  mutate(SiletzMill=SiletzMill*(2.5/6))%>%
  mutate(YaquinaMill=YaquinaMill*(1/6))
  #SiletzMill needs more weight in 1998 as Cascade is NA
  MSI_Three_Weight$SiletzMill[3] <- MSI_Three_Weight$SiletzMill[3]*2

MSI_Three_Weight <- MSI_Three_Weight%>%
  mutate(MSAdjThree = rowSums(MSI_Three_Weight[ , c(2,3,4)], na.rm=TRUE))
MSI_Three <- MSI_Three_Weight%>%
  dplyr::select(-c("Cascade", "SiletzMill", "YaquinaMill"))
MSI_Three$MSAdjThree <- (exp(MSI_Three$MSAdjThree)/(exp(MSI_Three$MSAdjThree)+1))
MSI_ESU <- left_join(MSI_ESU, MSI_Three)

print(ggplot(data=MSI_ESU, aes(x=ReturnYear, y=MSAdjCurrent)) +
  xlab("Return Year") +
  ggtitle(paste("MSI Alternatives")) +
  geom_line(color="black",linewidth=1.5) + 
  geom_point(color="black", size = 2) +
  geom_line(aes(y=MSAdjMARSS), col='orange', linewidth=1.5) +
  geom_point(aes(y=MSAdjMARSS), col='orange', size = 2) +
  geom_line(aes(y=MSAdjThree), col='blue', linewidth=1.5) +
  geom_point(aes(y=MSAdjThree), col='blue', size = 2) +
  geom_line(aes(y=MSAdjOG), col='green', linewidth=1.5) +
  geom_point(aes(y=MSAdjOG), col='green', size = 2) +
  theme_bw())



#Progressive cross-validation
#Run the MARSS alternative with the 1 to 8 years of MS LCM data removed
for (i in 1:8){
  #Remove WF Smith and Winchester from the dataset
  MSI_data_cv <- MSI_data %>%
    filter(!((between(ReturnYear, (2025-i), 2024)) & (Site == "West Fork Smith" | Site == "Winchester Creek")))
  
  MSI_table_cv<-MSI_data_cv%>%
    pivot_wider(names_from = Site,values_from = MSAdj,id_cols=ReturnYear)
  
  #Alternative One: use MARSS to interpolate missing values
  marss_mat<-MSI_table_cv%>%
    dplyr::select(-c(ReturnYear))%>%
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
  
  MSI_fitted<-MSI_table_cv%>%
    bind_cols(fitted)%>%
    pivot_longer(cols=c(everything(), -ReturnYear),names_to = "Site")%>%
    mutate(type=ifelse(grepl("mle_",Site),"mle","obs"))%>%
    mutate(Site=gsub("mle_","",Site))%>%
    filter(type=="mle")%>%
    group_by(ReturnYear,Site)%>%
    summarise(MSAdj=mean(value))
  
  MSI_fitted <- MSI_fitted%>%
    group_by(ReturnYear)%>%
    summarise(MSAdjMARSScv=mean(MSAdj))
  MSI_fitted$MSAdjMARSScv <- (exp(MSI_fitted$MSAdjMARSScv)/(exp(MSI_fitted$MSAdjMARSScv)+1))
  colnames(MSI_fitted)[2] <- paste0("MSAdjMARSSm", i)
  nrow(MSI_fitted)
  MSI_fitted[(nrow(MSI_fitted)-(i-1)):nrow(MSI_fitted),]
  
  MSI_ESU <- left_join(MSI_ESU, MSI_fitted[(nrow(MSI_fitted)-(i-1)):nrow(MSI_fitted),])
  
  MSline <- sym(paste0("MSAdjMARSSm", i))
  print(ggplot(data=MSI_ESU, aes(x=ReturnYear, y=MSAdjCurrent)) +
    xlab("Return Year") +
    ggtitle(paste("OCV minus", i)) +
    geom_line(color="black",linewidth=1.5) + 
    geom_point(color="black", size = 2) +
    geom_line(aes(y=!!MSline), col='orange', linewidth=1.5) +
    geom_point(aes(y=!!MSline), col='orange', size = 2) +
    geom_line(aes(y=MSAdjThree), col='blue', linewidth=1.5) +
    geom_point(aes(y=MSAdjThree), col='blue', size = 2) +
    geom_line(aes(y=MSAdjOG), col='green', linewidth=1.5) +
    geom_point(aes(y=MSAdjOG), col='green', size = 2) +
    theme_bw())
  
}


#Create evaluation table
MSI_eval_columns <- c("Years_removed", "MARSS_MAE", "Three_site_MAE", "MARSS_MAPE", "Three_site_MAPE")
MSI_eval <- data.frame(matrix(nrow = 0, ncol = length(MSI_eval_columns)))
for (i in 1:8){
  startrow <- nrow(MSI_ESU)+1-i
  endrow <- nrow(MSI_ESU)
  MARSScol <- 5+i
  MSI_eval <- rbind(MSI_eval, c(i, mae(MSI_ESU$MSAdjCurrent[startrow:endrow], pull(MSI_ESU[startrow:endrow,MARSScol])), mae(MSI_ESU$MSAdjCurrent[startrow:endrow], MSI_ESU$MSAdjThree[startrow:endrow]), mape(MSI_ESU$MSAdjCurrent[startrow:endrow], pull(MSI_ESU[startrow:endrow,MARSScol])), mape(MSI_ESU$MSAdjCurrent[startrow:endrow], MSI_ESU$MSAdjThree[startrow:endrow])))
}
colnames(MSI_eval) <- MSI_eval_columns


print(ggplot(data=MSI_eval, aes(x=Years_removed, y=MARSS_MAE)) +
        xlab("Years Removed") +
        ggtitle("MAE") +
        geom_line(color="black",linewidth=1.5) + 
        geom_point(color="black", size = 2) +
        geom_line(aes(y=Three_site_MAE), col='orange', linewidth=1.5) +
        geom_point(aes(y=Three_site_MAE), col='orange', size = 2) +
        theme_bw())

print(ggplot(data=MSI_eval, aes(x=Years_removed, y=MARSS_MAPE)) +
        xlab("Years Removed") +
        ggtitle("MAPE") +
        geom_line(color="black",linewidth=1.5) + 
        geom_point(color="black", size = 2) +
        geom_line(aes(y=Three_site_MAPE), col='orange', linewidth=1.5) +
        geom_point(aes(y=Three_site_MAPE), col='orange', size = 2) +
        theme_bw())

#Print the evaluation table
MSI_eval_table <- gt(MSI_eval) %>%
  tab_header(
    title = "MSI alternative evaluation",
  ) %>%
  fmt_number(
    columns = c(MARSS_MAE, Three_site_MAE, MARSS_MAPE, Three_site_MAPE),
    decimals = 3
  ) %>%
  cols_label(
    Years_removed = "Years Removed",
    MARSS_MAE = "MAE MARSS",
    Three_site_MAE = "MAE Three Sites",
    MARSS_MAPE = "MAPE MARSS",
    Three_site_MAPE = "MAPE Three Sites"
  ) %>%
  tab_style(
    style = list(
      cell_text(weight = "bold")
    ),
    locations = cells_column_labels(everything())
  )
print(MSI_eval_table)
