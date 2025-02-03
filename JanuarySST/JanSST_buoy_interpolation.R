#Creator: Erik Suring (erik.suring@odfw.oregon.gov)

#OCN Covariate Charleston Water Temperature (CHAO3) January SST is sometimes unavailable due to sensor loss
#This code uses other buoy data to predict the CHAO3 January SST
#Previous years used only buoy 46050, but sensor loss in 2025 necessitated using additional buoys

library(dplyr)
library(lubridate)

#OPITT buoy dataset from Mark Sorel with mle interpolation for missing data 
buoy_sst_d <- read.csv("buoy_interp.csv")
#The OCN CHAO3 data is a monthly average - aggregate buoy data to a monthly January average
buoy_sst_d <- buoy_sst_d %>% 
  mutate(year = year(date))
buoy_sst_d <- buoy_sst_d %>% 
  mutate(month = month(date))
buoy_sst_d <- buoy_sst_d %>% 
  mutate(day = day(date))
buoy_sst_d <- buoy_sst_d %>% filter(month==1 & type=="mle" & year>2010)
buoy_sst_m <- aggregate(value~buoyid+year+month,buoy_sst_d,mean)
buoy_sst_jan <- tapply(buoy_sst_m$value, list(buoy_sst_m$year, buoy_sst_m$buoyid), mean)

#Charleston Water Temperature monthly average created by OCN forecast data download code
CWT <- read.csv("CWT.mon.csv")
#Match CWT timeseries to buoy data
CWT_jan <- CWT %>% filter (Year>2010)
CWT_jan <- CWT_jan[,1:2]
#Set up result dataframes
corTemp <- data.frame(matrix(NA, ncol = 4, nrow = 1))
colnames(corTemp) <- colnames(buoy_sst_jan)
predTemp <- data.frame(matrix(NA, ncol = 6, nrow = nrow(CWT_jan)))
colnames(predTemp) <- c("year", colnames(buoy_sst_jan), "weighted_avg")
predTemp$year <- CWT_jan$Year

#Use correlation for weights
for(i in 1:4){
  corTemp[i] <- cor(CWT_jan$Jan, buoy_sst_jan[,i], use = "complete.obs")
  }

#Make prediction based on linear relationship between years with shared data
for(j in 1:4){
  y <- CWT_jan$Jan
  x <- buoy_sst_jan[,j]
  predx <- as.data.frame(buoy_sst_jan[,j])
  model <- lm(y ~ x)
  predTemp[,j+1] <- round(predict(model, predx),2)
  plot(x,y)
  plot(model$residuals)
}

#Create weighted average
for(k in 1:nrow(predTemp)){
  predTemp[k,6] <- round(weighted.mean(predTemp[k,2:5], corTemp),2)  
}

predTemp[k,6]
plot(CWT_jan$Jan,predTemp$weighted_avg)
