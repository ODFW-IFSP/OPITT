### Estimate confidence intervals on forecast by quasi-bootstrapping
### David E. Rupp
### 28 Oct 2010
 
ProcessStartTime = proc.time() 
print("Generating prediction intervals...")

nMod = length(Var1) 
#Define time period for initial calibration
StartYr = 1999
EndYr = 2013
#strtRow = StartYr - (1950 + RM_Yrs - 1) + 1
strtRow = 30
#endRow = EndYr - (1950 + RM_Yrs - 1) + 1
endRow = 43
endRow0 = endRow
nSteps = forecastYr - EndYr 
FC_Years = (EndYr+1):forecastYr  #forecastYr is last year in data array
theObs = rep(NA,nSteps)              
theForecast = array(NA,c(nMod,nSteps))
theSE = array(NA,c(nMod,nSteps))
theSE_II = array(NA,c(nMod,nSteps))
#theError = rep(NA,nSteps)

#lstRow = forecastYr - (1950 + RM_Yrs - 1) + 1
lstRow = 44
y = forecastData[strtRow:lstRow,2]      #Change Y here  
yobs = y
#if (logFlag) yobs = exp(yobs)  
theYears = forecastData[strtRow:lstRow,1]

#Array for ensemble mean
ensembleAveFC = array(NA,c(length(FC_Years),5))
colnames(ensembleAveFC) = c("Year","Obs","Forecast","lowerCL","upperCL")
ensembleAveFC[,1] = FC_Years
ensembleAveFC[,3] = 0
countMod = rep(0,nSteps)

#Arrays for randomized forecasts
B = 1000000 #Number of resamples
randFC = array(NA,c(B,nMod))
allFCs = rep(NA,c(nMod*B)) 
randCIs = array(NA,c(2,nMod,length(FC_Years)))

leg.txt = dataFlag 
theColor = 1   
theColors = 1                      
#dev.new()
#win.metafile(file="xxfcplot.emf",width=6,height=6)
filename = paste("Forecast_with_CF.pdf",sep="")
#pdf(file=filename,width=7.5,height=10)

theTitle=paste("Ensemble mean forecast with 90% prediction intervals")
plot(theYears,yobs, type="b",xlab="Year",ylab=theYlabel,lwd=2,ylim=theylims,cex.axis=0.8,cex.lab=0.9,yaxt="n",xaxt="n",main=theTitle,cex.main=0.95) 
axisTicksX1 = seq(1970,2015,1)
axisTicksX2 = seq(1970,2015,5) 
axisTicks1 = c(seq(10,90,by=10),seq(100,900,by=100),seq(1000,9000,by=1000),10000)
axisTicks2 = c(10,50,100,500)       
tickLabels = c("10","50","100","500")
axis(2, at=axisTicks1, labels=F) # draw y axis with required labels 
axis(2, at=axisTicks2, labels=tickLabels) # draw y axis with required labels  
axis(1, at=axisTicksX1, labels=F) # draw y axis with required labels 
axis(1, at=axisTicksX2, labels=T) # draw y axis with required labels            
                                  
endRow = endRow0 

for (iStep in 1:nSteps) {

  #print("iStep")
  #print(iStep)   
  
  for (iMod in 1:nMod) {

    #print(iMod)
  
    i = Var1[iMod]
    j = Var2[iMod]
    icol = i 
    jcol = j
                
    #Subset of data for model calibration
    y = forecastData[strtRow:endRow,2]      #Change Y here
    ndata = length(y)
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
    yres = y - ypred  
    var_res = sum(yres^2)/(ndata-2)
    sd_res = sqrt(var_res) 
                                        
    #Subset of data for forecast
    y = forecastData[endRow+1,2]      #Change Y here
    theYears = forecastData[endRow+1,1] 
    #print(theYears)      
    x1 = forecastDataRM[endRow+1,icol] 
    x2 = forecastData[endRow+1,jcol]
    
    x = cbind(x1,x2)
    if (is.na(sum(x))) next #skip to next model if predictor variable missing
    
    dummy <- data.frame(y = y, x = x) 
    
    # Make forecast and get standard error
    output = predict(GAModel, dummy, se=TRUE)
    yfc = as.double(output[1])
    theSE[iMod,iStep] = as.double(output[2])
    theSE_II[iMod,iStep] = sqrt(var_res + (theSE[iMod,iStep])^2)
    theObs[iStep] = y                                 
    theForecast[iMod,iStep] = yfc
    ensembleAveFC[iStep,3] = ensembleAveFC[iStep,3] + theForecast[iMod,iStep]
 
    ### Generata B randomized future outcomes
    #randFC[,iMod] = rnorm(B,mean = yfc, sd = theSE[iMod,iStep])
    randFC[,iMod] = rnorm(B,mean = yfc, sd = theSE_II[iMod,iStep])

    #Combine randomized forecasts from each model into one set
    allFCs[((iMod-1)*B+1):(iMod*B)] = randFC[,iMod]
    countMod[iStep] = countMod[iStep] + 1
    
  }  #iMod loop end                                                                              

  endRow = endRow + 1   
  alpha = 0.1
  #Lower and upper quantiles of ensemble mean forecast from combining randomized forecast from all models
  theQuantiles = quantile(allFCs, probs = c(alpha/2,(1-alpha/2)), na.rm = TRUE, names = TRUE, type = 6)  
  ensembleAveFC[iStep,4] = theQuantiles[1]
  ensembleAveFC[iStep,5] = theQuantiles[2] 

} #iStep loop end


#Ensemble mean forecast
ensembleAveFC[,3] = ensembleAveFC[,3] / countMod                                 
ensembleAveFC[1:nSteps,2] = theObs[1:nSteps]

toPrint = paste("90% prediction intervals on ensemble mean forecast for ",theYears[nSteps],":",sep="")
print("90% prediction intervals on ensemble mean forecast:")
toPrint = c(exp(ensembleAveFC[nSteps,4]),exp(ensembleAveFC[nSteps,5]))
print(toPrint)

# Plot ensemble mean forecast and confidence intervals
#Forecast
lines(ensembleAveFC[,1],exp(ensembleAveFC[,3]), col=4)
points(ensembleAveFC[,1],exp(ensembleAveFC[,3]),col=4,pch=16)
#Lower CI
lines(ensembleAveFC[,1],exp(ensembleAveFC[,4]), col=4,lty=2)
points(ensembleAveFC[,1],exp(ensembleAveFC[,4]),col=4,pch=3,cex=0.7)
#Upper CI
lines(ensembleAveFC[,1],exp(ensembleAveFC[,5]), col=4,lty=2)
points(ensembleAveFC[,1],exp(ensembleAveFC[,5]),col=4,pch=3,cex=0.7)

                                                          
#Draw legend
abline(v=EndYr+0.5, lwd=2)      
#thelty = rep(1,(length(Var1)+1))
if (length(Var1)>=8) thelty[9:(length(Var1)+1)] = 2
thelwd = rep(1,(length(Var1)+1))
thelwd[1] = 2
theColors = c(1,4,4,3,2)
thelty = c(1,1,2,3,4)
leg.txt = c("OCN Rivers","Ensemble mean forecast","90% pred. int.")
#xPos = xPos - 1
legend(xPos, yPos, leg.txt, col = theColors,cex=0.75, lty=thelty, lwd=thelwd,, bty="o", bg="white")

dev.off()

write.table(ensembleAveFC, file=here::here("results/ForecastsandCIs.txt"),row.names = FALSE)

# Various distribution plots for confidence interval estimation
#dev.new()
#par(mfrow = c(3, 2), pty = "s")
#for (iMod in 1:nMod) {
#  theTitle <- paste("Forecasts, GAM #", iMod)
#  hist(randFC[,iMod], freq=FALSE, breaks = B/10000, col="lightblue", border="blue", xlab="log Recruits",ylab="frequency",main = theTitle, plot=TRUE)
#} 
#
#dev.new()
#par(mfrow = c(1, 1), pty = "s")
#theTitle <- paste("Ensemble forecasts")
#
#hist(allFCs, freq=FALSE, breaks = B/10000, col="lightblue", border="blue", xlab="log Recruits",ylab="frequency",main = theTitle, plot=TRUE)
#
#
#dev.new()
#par(mfrow = c(1, 1), pty = "s")
#
#x = seq(0,10,0.001)
#y = dnorm(x,mean=theForecast[1,nSteps],sd=theSE_II[1,nSteps])
#theColor=1
#theXlims = c((ensembleAveFC[nSteps,3]-1.4),(ensembleAveFC[nSteps,3]+1.4))
#plot(x,y,col=theColor,type="l",xlim=theXlims,xlab="log Recruits",ylab="Density (%)")
#for (iMod in 2:nMod) {
#  theColor = theColor + 1
#  y = dnorm(x,mean=theForecast[iMod,nSteps],sd=theSE_II[iMod,nSteps])
#  lines(x,y,col=theColor)
#}
#histData = hist(allFCs, breaks = B/5000, plot=FALSE)
#x = histData$mids
#y = histData$density
#lines(x,y,col=1,type="l",lwd=2)


#source("Plot_Conf_Intervals_3Var.R")

ElapsedTime = proc.time() - ProcessStartTime
ElapsedTime = ElapsedTime / 60 # minutes
print("Elapsed time for generation of  prediction intervals (minutes)...")
print(ElapsedTime)
