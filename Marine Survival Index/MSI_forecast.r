require(stats); require(graphics)
library(mgcv)

#Number of Years over which to calculate running mean
nKnots = 3
loopFlag = FALSE #if TRUE, loop through all possible combinations of three variables
allplots = FALSE
spawnFlag = TRUE





  fishfile = "OCNForecastData.csv"
  #fname = paste("../../Data/", fishfile, sep = "")
  OCNData <- read.csv(fishfile)
  theYlabel = "Abundance (thousands)"
  theTitle ="OCN Rivers 1970-2009"
  xPos = 1969 #1978
  yPos = 40 #55 #500  

InData = OCNData

#Convert ratios to logit
InData[,2] = log(InData[,2]/(1-InData[,2]))
InData[,3] = log(InData[,3]/(1-InData[,3]))
#Covert abundances to log
InData[,4] = log(InData[,4])
InData[,5] = log(InData[,5])

Nyrs = length(OCNData[,4])
earliestYear = OCNData[1,1]
latestYear = OCNData[Nyrs,1]

#My additions to use the small data subset
modelData = InData
modelDataRM = InData
Var.allnames = names(InData)

nYears = length(modelData[,1])
forecastYr = modelData[nYears,1]


  # Trim off last year which has no dependent variable (NA) and save for forecast
  forecastData = modelData
  forecastDataRM = modelDataRM  
  modelData = modelData[1:(nYears-1),]
  modelDataRM = modelDataRM[1:(nYears-1),]  
  
      Var1 = c(  3,   6,  3,  3,  6,  7,  6)
      Var2 = c(  6,   7, 14,  9, 11, 14, 10)

  thisYrForecast = rep(NA,length(Var1))
  xCorr = array(NA, c(2,2,length(Var1)))
  theStats0 = array(NA,c(length(Var1),12))
  colnames(theStats0) = c("Var1","Var2","Var3","GCV","AIC","R2","OCV","HFS","FC","EDF1","EDF2","EDF3")
  theStats = as.data.frame(theStats0)
  
  for (iPlot in 1:length(Var1)) {
  #iPlot = 1
    i = Var1[iPlot]
    j = Var2[iPlot]
    theStats[iPlot,1] = Var.allnames[i]
    theStats[iPlot,2] = Var.allnames[j]
    icol = i  
    jcol = j  
        
    y = modelData[!is.na(modelData[,2]),2]  
	yobs = y
    theYears = modelData[!is.na(modelData[,2]),1]
    x1 = modelDataRM[!is.na(modelData[,2]),icol] 
    x2 = modelData[!is.na(modelData[,2]),jcol]
     
    x = cbind(x1,x2)
                              
    xCorr[,,iPlot] = cor(x, use = "pairwise.complete.obs")
    print(xCorr[2])
   
    dummy <- data.frame(y = y, x = x) 
    
    if (iPlot == 1){
      ensembleAve = array(0,c(length(y),4))
      TheResiduals = array(0,c(length(y),length(Var1)))
    }
    
    if (!allplots) {
      if (iPlot == 1){
        yobs = (exp(y)/(exp(y)+1))
        leg.txt = dataFlag 
        theColor = 1
        theColors = 1                      
        #dev.new()
        #win.metafile(file="xxplot.emf",width=6,height=6)
        #win.graph(width = 5, height = 8)
        pdf(file='FittedGAMS.pdf', width = 5, height = 8)
        par(mfrow = c(2, 1))
        par(oma=c(0,0,1,1)) 
        par(mar=c(4.5,4,0,0))
        
        theylims = c(0.0000001,0.2)
        theTitle="" #"Fitted GAMs"
        plot(theYears,yobs, type="b",xlab="",ylab=theYlabel,lwd=2,ylim=theylims,cex.axis=0.8,cex.lab=0.9,log="y",yaxt="n",xaxt="n",main=theTitle) 
        axisTicksX1 = seq(1970,2015,1)
        axisTicksX2 = seq(1970,2015,5)        
        axisTicks1 = c(seq(1,9,by=1),seq(10,90,by=10),seq(100,900,by=100),seq(1000,9000,by=1000),10000)
        axisTicks2 = c(10,50,100,500)       
        tickLabels = c("10","50","100","500")  
        axis(2, at=axisTicks1, labels=F) # draw y axis with required labels 
        axis(2, at=axisTicks2, labels=tickLabels, cex.axis=0.8) # draw y axis with required labels  
        axis(1, at=axisTicksX1, labels=F) # draw y axis with required labels 
        axis(1, at=axisTicksX2, labels=T, cex.axis=0.8) # draw y axis with required labels 
        
        ensembleAve = array(0,c(length(y),4))    #Array for ensemble averaging
        ensembleAve[,1] = theYears
#       ensembleAve[,2] = log(yobs)
        ensembleAve[,2] = y
        
      } 
    }  
    
    if (i == 6) {
      #legtemp = paste("PDO.MJJ-4","+",Var.allnames[j], sep="")
      legtemp = paste(Var.allnames[j], sep="")
    } else {
      legtemp = paste(Var.allnames[i],"+",Var.allnames[j], sep="")    
    }
     
    source(here::here"functions/TheGAM_2Var.R")
    ensembleAve[,3] = ensembleAve[,3]+ypred
    TheResiduals[,iPlot] = yres 
    source(here::here"functions/OCV_GAM_2Var.R")
    ensembleAve[,4] = ensembleAve[,4]+y_fc
        
    #Add model predictions to array of all model predictions.  These are the GAM fits, NOT the forecasts
    if (iPlot == 1) {
        allPredictions = cbind(theYears,yobs,ypred)
    } else {
        allPredictions = cbind(allPredictions,ypred)
    }    
    
    if (!allplots) {
      leg.txt = c(leg.txt, legtemp) 
    }         
  }
  
  #Calculate ensemble averages and their R2 and OCV*  
  ensembleAve[,3:4] = ensembleAve[,3:4] / length(Var1)
  meany = mean(ensembleAve[,2])
  ss1 = sum((ensembleAve[,2] - meany)^2)
  ss2 = sum((ensembleAve[,2] - ensembleAve[,3])^2)
  R2 = 1 - ss2/ss1
  R2_cal = round(R2,2)
  ss2 = sum((ensembleAve[,2] - ensembleAve[,4])^2)
  OCV = 1 - ss2/ss1
  OCV_cal = round(OCV,2)
  print("Ensemble average R2:")
  print(R2_cal)
  print("Ensemble average OCV*:")  
  print(OCV_cal)
  ensembleAve[,2:4] = exp(ensembleAve[,2:4])/(exp(ensembleAve[,2:4])+1)
  
  if (!allplots) {
    thelty = rep(1,(length(Var1)+1))
    if (length(Var1)>=8) thelty[9:(length(Var1)+1)] = 2
    if (jackFlag)  thelty[2] = 2
    thelwd = rep(1,(length(Var1)+1))
    thelwd[1] = 2
    legend(xPos, yPos, leg.txt, col = theColors,cex=0.65, lty=thelty, lwd=thelwd)
  }

  ResMean = rowMeans(TheResiduals)

  ## Plot ACF of residuals  
  dev.new()
  acf0 = acf(ResMean, lag.max=14, main="Autocorrelation function of residuals", cex.main=1.5, cex.lab=1.5,ci=0.95, ci.type = "white", xlab="lag (year)")

  win.graph(width = 6, height = 8)
  par(mfrow = c(3, 3))
  par(oma=c(0,0,1,1)) 
  par(mar=c(4.5,4,0,0))
  for (iPlot in 1:length(Var1)) {
      acf0 = acf(TheResiduals[,iPlot], lag.max=14, main="Autocorrelation function of residuals", cex.main=1.5, cex.lab=1.5,ci=0.95, ci.type = "white", xlab="lag (year)")
  }
  #stop("m")

    #source(here::here"functions/Forecast_Skill_2Var.R")
    #Calculate HFS of ensemble averages
    meany = mean(ensembleAveFC[1:(length(ensembleAveFC[,1])-1),2])
    ss1 = sum((ensembleAveFC[1:(length(ensembleAveFC[,1])-1),2] - meany)^2)
    ss2 = sum((ensembleAveFC[1:(length(ensembleAveFC[,1])-1),2] - ensembleAveFC[1:(length(ensembleAveFC[,1])-1),3])^2)
    HFS = 1 - ss2/ss1
    HFS_cal = round(HFS,2)
    print("Ensemble average HFS:")
    print(HFS_cal)
    ensembleAveFC[,2:3] = exp(ensembleAveFC[,2:3])/(exp(ensembleAveFC[,2:3])+1)
    
    theText = paste("Forecasts for ", forecastYr,":", sep="")
    print(theText)
    #print(exp(thisYrForecast))
    print(exp(thisYrForecast)/(exp(thisYrForecast)+1))
    tempFC = thisYrForecast # log(thisYrForecast)
    meanFC = mean(tempFC,na.rm=TRUE)
    #meanFC = exp(meanFC)
    meanFC = (exp(meanFC)/(exp(meanFC)+1))
    theText = paste("Ensemble mean forecast for ", forecastYr,":", sep="")
    print(theText)
    print(meanFC)

    source(here::here"functions/Conf_Intervals_2Var.R")  

fullStats = theStats
theStats = cbind(fullStats[,1:3],(exp(fullStats[,9])/(exp(fullStats[,9])+1)),fullStats[,6],fullStats[,7])
colnames(theStats) = c("Var1", "Var2", "Var3", "Forecast", "R2", "OCV")

cat("Forecast year:",forecastYr, "\n", file = "ForecastStatsByModel.csv", sep = ",", fill = FALSE, labels = NULL, append = FALSE)
cat("Fish data file:",fishfile, "\n", file = "ForecastStatsByModel.csv", sep = ",", fill = FALSE, labels = NULL, append = TRUE)
cat("Predictor data file:",fishfile, "\n", file = "ForecastStatsByModel.csv", sep = ",", fill = FALSE, labels = NULL, append = TRUE)
write.table(theStats, file="ForecastStatsByModel.csv", row.names = FALSE, sep=",", append = TRUE)
cat("Ensemble Mean"," ", " ", meanFC, "\n", file = "ForecastStatsByModel.csv", sep = ",", fill = FALSE, labels = NULL, append = TRUE)
cat("Lower 90% P.I.", " ", " ", exp(ensembleAveFC[nSteps,4]), "\n", file = "ForecastStatsByModel.csv", sep = ",", fill = FALSE, labels = NULL, append = TRUE)
cat("Upper 90% P.I.", " ", " ", exp(ensembleAveFC[nSteps,5]), "\n", file = "ForecastStatsByModel.csv", sep = ",", fill = FALSE, labels = NULL, append = TRUE)

write.table(fullStats, file="bestmodelsStats.txt",row.names = FALSE)
write.table(allPredictions, file="Predictions.txt",row.names = FALSE)

colnames(AllPredictors2) = var.names2
dev.off()


