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
MSI_table<-MSI_data%>%
  pivot_wider(names_from = Site,values_from = MSAdj,id_cols=ReturnYear)
returnYear <- MSI_table$ReturnYear

marss_mat<-MSI_table%>%
  dplyr::select(-c(ReturnYear))%>%
  as.matrix()%>%
  t()

#Alternative MARSS models

#Model.1 default model structure
model.1=list(
  Q= "diagonal and unequal",
  R= "diagonal and equal",
  U= "unequal"
)
fit.1=MARSS(marss_mat, model=model.1,control=list(maxit=1000,allow.degen=T))

#Model.2 
model.2=list(
  Q= "equalvarcov",
  R= "diagonal and equal",
  U= "unequal"
)
fit.2=MARSS(marss_mat, model=model.2,control=list(maxit=1000,allow.degen=T))

#Model.3 
model.3=list(
  Q= "equalvarcov",
  R= "diagonal and unequal",
  U= "unequal"
)
fit.3=MARSS(marss_mat, model=model.3,control=list(maxit=1000,allow.degen=T))

#Model.4 
model.4=list(
  Q= "equalvarcov",
  R= diag(rep(0,nrow(marss_mat))),
  U= "unequal"
)
fit.4=MARSS(marss_mat, model=model.4,control=list(maxit=1000,allow.degen=T))

#Model.5 
model.5=list(
  Q= "equalvarcov",
  R= "diagonal and unequal",
  U= matrix(rep(0,nrow(marss_mat)),nrow=nrow(marss_mat),1)
)
fit.5=MARSS(marss_mat, model=model.5,control=list(maxit=1000,allow.degen=T))

#Model.6 
model.6=list(
  Q= "equalvarcov",
  R= diag(rep(0,nrow(marss_mat))),
  U= matrix(rep(0,nrow(marss_mat)),nrow=nrow(marss_mat),1)
)
fit.6=MARSS(marss_mat, model=model.6,control=list(maxit=1000,allow.degen=T))

#Model.7 
model.7=list(
  Q= "equalvarcov",
  R= "diagonal and equal",
  U= matrix(rep(0,nrow(marss_mat)),nrow=nrow(marss_mat),1)
)
fit.7=MARSS(marss_mat, model=model.7,control=list(maxit=1000,allow.degen=T))

c(fit.1$AICc, fit.2$AICc, fit.3$AICc, fit.4$AICc, fit.5$AICc, fit.6$AICc, fit.7$AICc)

resids <- MARSSresiduals(fit.2, type = "tt1")
for (i in 1:6) {
  plot(resids$model.residuals[i, ], ylab = "model residuals", 
       xlab = "")
  abline(h = 0)
  title(rownames(marss_mat)[i])
}

for (i in 1:6) {
  plot(returnYear, fit.2$states[i, ], ylab = "log subpopulation estimate", 
       xlab = "", type = "l")
  lines(returnYear, fit.2$states[i, ] - 1.96 * fit.2$states.se[i, 
  ], type = "l", lwd = 1, lty = 2, col = "red")
  lines(returnYear, fit.2$states[i, ] + 1.96 * fit.2$states.se[i, 
  ], type = "l", lwd = 1, lty = 2, col = "red")
  title(rownames(marss_mat)[i])
}

