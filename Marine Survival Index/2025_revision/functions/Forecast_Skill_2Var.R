#Define time period for initial calibration
StartYr = 1999
EndYr = 2013
#strtRow = StartYr - (1950 + RM_Yrs - 1) + 1
strtRow = 30
#endRow = EndYr - (1950 + RM_Yrs - 1) + 1
endRow = 44
endRow0 = endRow
nSteps = forecastYr - EndYr 
FC_Years = (EndYr+1):forecastYr  #forecastYr is last year in data array
theObs = rep(NA,nSteps)              
theForecast = rep(NA,nSteps)
theError = rep(NA,nSteps)
theNS = rep(NA,nSteps)
theGCV = rep(NA,nSteps)

#lstRow = forecastYr - (1950 + RM_Yrs - 1) + 1
lstRow = 44
y = forecastData[strtRow:lstRow,2]                  #Change Y here
yobs = y
#if (logFlag) yobs = exp(yobs)  
theYears = forecastData[strtRow:lstRow,1]

#Arrays for all forecast and ensemble mean
AllFCs = array(NA,c(length(FC_Years),length(Var1))) 
ensembleAveFC = array(0,c(length(FC_Years),4)) 
colnames(ensembleAveFC) = c("Year","Obs","Forecast","Weighted_Forecast")
ensembleAveFC[,1] = FC_Years
ensembleAveFC[,2] = NA

x_save = array(NA,c(length(Var1),2))
yfc_save = rep(NA,length(Var1))
yfc_comp = array(NA,c(length(Var1),2))
TheIntercept = rep(NA,length(Var1))
TheAICs = array(NA,c(length(FC_Years),length(Var1)))
AkaikeWeights = array(NA,c(length(FC_Years),length(Var1)))

leg.txt = dataFlag 
theColor = 1   
theColors = 1                    
#dev.new()
#win.metafile(file="xxfcplot.emf",width=6,height=6)
theTitle="" #"Model forecasts (1990 to current year)"
plot(theYears,yobs, type="b",xlab="Year",ylab=theYlabel,lwd=2,ylim=theylims,cex.axis=0.8,cex.lab=0.9,yaxt="n",xaxt="n",main=theTitle) 
axisTicksX1 = seq(1970,2015,1)
axisTicksX2 = seq(1970,2015,5) 
axisTicks1 = c(seq(1,9,by=1),seq(10,90,by=10),seq(100,900,by=100),seq(1000,9000,by=1000),10000)
axisTicks2 = c(10,50,100,500)       
tickLabels = c("10","50","100","500")
axis(2, at=axisTicks1, labels=F) # draw y axis with required labels 
axis(2, at=axisTicks2, labels=tickLabels, cex.axis=0.8) # draw y axis with required labels 
axis(1, at=axisTicksX1, labels=F) # draw y axis with required labels 
axis(1, at=axisTicksX2, labels=T, cex.axis=0.8) # draw y axis with required labels             
                               
varCount = rep(0,nSteps)
               
## Data for first predictor
i = Var1[1]
TheListofVars2 = i
icol = i
var.names2 = Var.allnames[i]
AllPredictors2 = forecastDataRM[strtRow:lstRow,icol]   
               
for (iPlot in 1:length(Var1)) {
#iPlot=1
  i = Var1[iPlot]
  j = Var2[iPlot]
  icol = i 
  jcol = j 

  #if (i == 6) {
  #  #legtemp = paste("PDO.MJJ-4","+",Var.allnames[j],"+",Var.allnames[k], sep="") 
  #  legtemp = paste(Var.allnames[j],"+",Var.allnames[k], sep="") 
  #} else {
    legtemp = paste(Var.allnames[i],"+",Var.allnames[j], sep="") 
  #} 
  endRow = endRow0

  theColor = theColor + 1
  if (theColor == 7) {
    theColor = theColor + 1
  } 
  
  for (iStep in 1:(nSteps+1)) {
  #iStep=1       
    #Subset of data for model calibration
    y = forecastData[strtRow:endRow,2]                #Change Y here
    theYears = forecastData[strtRow:endRow,1]       
    x1 = forecastDataRM[strtRow:endRow,icol] 
    x2 = forecastData[strtRow:endRow,jcol]
    x = cbind(x1,x2)  
                              
    dummy <- data.frame(y = y, x = x)  
    
    #Calibrate GAM
   
    #nKnots = 4
    GAModel <- gam(y~s(x.x1,k=nKnots)+s(x.x2,k=nKnots),data=dummy)     
    modSum = summary(GAModel)
    ypred = predict(GAModel, dummy)
    
    #NS Efficiency
    meany = mean(y)
    ss1 = sum((y - meany)^2)
    ss2 = sum((y - ypred)^2)
    NS_cal = 1 - ss2/ss1
    NS_cal = round(NS_cal,2)
    #NS_cal = modSum[[14]]
    
    #Plot the fit of the regression for the time period of initial calibration
    if (iStep == 1) {
      #if (logFlag) {
      #  ypred = exp(ypred)
      #}
      #legtemp = paste(legtemp,"  R2 = ", NS_cal, sep="") 
      thelty = 1
      if (iPlot >= 8) thelty = 2 
      #lines(theYears,ypred, col=theColor, lty=thelty, lwd=1)

#     dev.new()
#     par(mfrow=c(1,3),pty = "s")
#     plot(GAModel,pages=1,se=T,res=T,pch=1) #covariate effect
 
      if (sum(TheListofVars2 == j) == 0){
        var.names2 = c(var.names2,Var.allnames[j])
        AllPredictors2 = cbind(AllPredictors2,forecastData[strtRow:lstRow,jcol])  
        TheListofVars2 = c(TheListofVars2,j)
      }     
    }
    
    if (iStep <= nSteps) {
      TheAICs[iStep,iPlot] = GAModel$aic
      #Subset of data for forecast
      y = forecastData[strtRow:(endRow+1),2]          #Change Y here
      theYears = forecastData[strtRow:(endRow+1),1]       
      x1 = forecastDataRM[strtRow:(endRow+1),icol] 
      x2 = forecastData[strtRow:(endRow+1),jcol]
      x = cbind(x1,x2)
      dummy <- data.frame(y = y, x = x) 
      
      yfc = predict(GAModel, dummy)
      theNS[iStep] = round(NS_cal,3)
      theGCV[iStep] = GAModel[[24]]
      theForecast[iStep] = yfc[length(yfc)]
      theObs[iStep] = y[length(yfc)]
      theError[iStep] = theForecast[iStep] - theObs[iStep]  
      
      AllFCs[iStep,iPlot] = theForecast[iStep]
      
      if (!is.na(theForecast[iStep])) {                                                                    
        ensembleAveFC[iStep,3] = ensembleAveFC[iStep,3] + theForecast[iStep]
        varCount[iStep] = varCount[iStep] + 1
      }
      
      
      #Get data for additive effects plots
      if (iStep == nSteps - (2011 - 2011)) {
        x_save[iPlot,] = x[length(yfc),]
        yfc_save[iPlot] = theForecast[iStep]
        tmp_comp = predict(GAModel, dummy, type="terms")    
        yfc_comp[iPlot,] = tmp_comp[length(yfc),]
        TheAttributes = attributes(tmp_comp)
        TheIntercept[iPlot] = as.double(TheAttributes[3])
        
        if (iPlot==1) {
          GAM_list = list(rep(GAModel,length(Var1)))
          GAM_list[[iPlot]] = GAModel
        } else {
          GAM_list[[iPlot]] = GAModel
        }    

      }
      
      endRow = endRow + 1
    }
    
  }  #iStep loop end
  
  #RMSE
  nStp = nSteps - 1
  temp = theError^2
  RMSE = sum(temp[1:nStp])
  RMSE = RMSE / length(RMSE)
  RMSE = round(RMSE,3)
  rm(temp) 

  #NS Efficiency
  meany = mean(theObs[1:nStp], na.rm=TRUE)  
  ss1 = sum((theObs[1:nStp] - meany)^2, na.rm=TRUE)
  ss2 = sum((theObs[1:nStp] - theForecast[1:nStp])^2, na.rm=TRUE)
  NS = 1 - ss2/ss1
  theStats[iPlot,8] = NS
  #print("HFS:")
  #print(NS)
  NS = round(NS,2)
      
  #if (logFlag) {
  #  #theForecast = exp(theForecast)
  #  ypred = exp(ypred)
  #}

  thisYrForecast[iPlot] = theForecast[length(theForecast)]
  theStats[iPlot,9] = thisYrForecast[iPlot]
  
  thelty = 1
  if (iPlot >= 8) thelty = 2 
  lines(FC_Years,exp(theForecast), col=theColor, lty=thelty, lwd=1)                                                      

  ensembleAveFC[1:(length(theForecast)-1),2] = theObs[1:nStp]
  
  #Update text for legend
  theColors = c(theColors,theColor)
  legtemp = paste(legtemp, "  HFS =", NS)
  leg.txt = c(leg.txt, legtemp) 

} # iPlot (iVar) loop end

ensembleAveFC[,3] = ensembleAveFC[,3] / varCount

for (iStep in 1:nSteps) { 
  DeltaAIC = TheAICs[iStep,] - min(TheAICs[iStep,], na.rm=TRUE)
  tempcalc = exp(-DeltaAIC/2)
  Sumtempcalc = sum(tempcalc, na.rm=TRUE)
  AkaikeWeights[iStep,] = tempcalc / Sumtempcalc
  ensembleAveFC[iStep,4] = sum(AkaikeWeights[iStep,] * AllFCs[iStep,])
}


#Draw legend
abline(v=EndYr+0.5, lwd=2)
thelty = rep(1,(length(Var1)+1))
if (length(Var1)>=8) thelty[9:(length(Var1)+1)] = 2
thelwd = rep(1,(length(Var1)+1))
thelwd[1] = 2
#xPos = xPos - 1
legend(xPos, yPos, leg.txt, col = theColors,cex=0.65, lty=thelty, lwd=thelwd, bty="o", bg="white")
dev.off()

#write.table(ensembleAveFC, file="Forecasts.txt",row.names = FALSE)

# Partial regression plots with latest forecast
print("Additive effects...")
print(yfc_comp)

for (iPlot in 1:length(Var1)) {

  if (iPlot==1) {
    #dev.new()
    filename = paste(here::here("results/AdditiveEffectPlot.pdf"),sep="")
    pdf(file=filename,width=7.5,height=10)
    #filename = paste("AdditiveEffectPlot1.emf",sep="")
    #win.metafile(file=filename,width=7.5,height=10)
    #win.graph(width = 7.5, height = 10)
    par(mfrow=c(3,2),pty = "s")
    #par(mfrow=c(1,1),pty = "s")
  }  

  i = Var1[iPlot]
  j = Var2[iPlot]

  #dev.new()
  #filename = paste("AdditiveEffectPlot_",iPlot,".emf",sep="")
  #win.metafile(file=filename,width=7.5,height=2.5)
  #par(mfrow=c(1,3),pty = "s")
  #win.metafile(file="GAMplot.emf",width=6,height=6)
  #par(mfrow=c(3,3),pty = "s")
  #par(mar=c(5,5,1,1))
  #layout(matrix(1:9,3,3))

  #get forecast by components
  Term = rep(NA,2)
  yfc = yfc_save[iPlot]
  for (igraph in 1:2) {
    Term[igraph] = yfc_comp[iPlot,igraph]
    #print(Term[igraph])
  }

  #print("The Intercept")
  #print(TheIntercept[iPlot]))
  

  #if (i == 6) {
  #  theXlabels = c(Var.allnames[i],Var.allnames[j],Var.allnames[k])
  #} else {
    theXlabels = c(Var.allnames[i],Var.allnames[j]) 
  #}
  
  theYlims=NULL
  yText = "Additive effect"
  #if (iPlot > 1) yText = ""
  #for (igraph in 2:2) {
  for (igraph in 1:2) {
    #tempText = "PDO.MJJ-4"
    if (igraph > 0) {
      tempText = theXlabels[igraph]
      theYlims = c(-0.8,0.8)
      #yText = ""
      if (igraph == 2) {
        if (j == 127) {
          tempText = paste(theXlabels[igraph],"(day)") 
        } else if (j == 36) {
          tempText = expression(paste(UWI.JAS," (m"^3,"s"^{-1},"100m"^{-1},")", sep=""))          
        } else if (j == 38) {
          tempText = expression(paste(UWI.SON," (m"^3,"s"^{-1},"100m"^{-1},")", sep=""))          
        } else if (j == 61) {
          tempText = paste(theXlabels[2],"(mm)")
        }  
      } else if (igraph == 3) {
        if (k == 36) {
          tempText = expression(paste(UWI.JAS," (m"^3,"s"^{-1},"100m"^{-1},")", sep=""))  
        } else if (k == 75) {
          tempText = expression(paste(SST.AMJ," ("^o,"C)", sep=""))  
        }  else if (k == 84) {
          tempText = expression(paste(SST.J," ("^o,"C)", sep=""))  
        }
      }
    }
    theYlims = c(-1.5,1.25)

    if (igraph == 1) {      
      theTitle=paste("Predictor effect: model ",iPlot)
      #theTitle=paste("")      
      plot(GAM_list[[iPlot]],select=igraph,residuals=TRUE,rug=FALSE,se=FALSE,pages=0,pch=1,cex.lab=1.3,cex.axis=1.2,cex=0.9,cex.main=0.9,
        xlab=tempText,ylab=yText,ylim=theYlims,yaxt="n",main=theTitle,col=4) #covariate effect  
    } else {  
      plot(GAM_list[[iPlot]],select=igraph,residuals=TRUE,rug=FALSE,se=FALSE,pages=0,pch=1,cex.lab=1.3,cex.axis=1.2,cex=0.9,
        xlab=tempText,ylab=yText,ylim=theYlims,yaxt="n",col=4) #covariate effect  
    }

    axisTicks1 = seq(-2,2,0.5)      
    axisTicks2 = seq(-2,2,1)
    axis(2, at=axisTicks1, labels=F,cex.axis=1.2)
    axis(2, at=axisTicks2, labels=T,cex.axis=1.2)
    
    x = x_save[iPlot,igraph]
    y = Term[igraph]
    Xv = c(x,x)
    Yv = c(-2,y)
    lines(Xv,Yv,col=2)
    Xh = c(-200,x)
    Yh = c(y,y)
    lines(Xh,Yh,col=2)
      
  }

}
dev.off()

  

