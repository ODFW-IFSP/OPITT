#require(stats); require(graphics)
#library(mgcv)
library(scatterplot3d)
library(rgl)
   
#nKnots = 2
GAModel <- gam(y~s(x.x1,k=nKnots)+s(x.x2,k=nKnots),data=dummy)     
  
#print(GAModel)
modSum = summary(GAModel)
theStats[iPlot,10:11] = modSum$edf

#print(modSum)
#print("Smoothing parameters:")
#print(GAModel$sp)

theGCV = GAModel$gcv.ubre 
theAIC = GAModel$aic
theStats[iPlot,4] = theGCV
theStats[iPlot,5] = theAIC
pvalues = modSum[[8]]

ypred = predict(GAModel, dummy)
yobs = y
yres = ypred - y

#NS Efficiency
meany = mean(yobs)
ss1 = sum((y - meany)^2)
ss2 = sum((y - ypred)^2)
NS = 1 - ss2/ss1
theStats[iPlot,6] = NS
print("GAM Efficiency:")
print(NS)
NS = round(NS,2)
#Rsqrd[k] = NS
#if (logFlag) {
#  ypred = exp(ypred) 
#  yobs = exp(yobs)                       
#}

thefilename = paste("partialreg_",iPlot,".emf",sep="")

if (allplots) {
  
  #dev.new()
  #win.metafile(file=thefilename,width=6,height=3.5)
  #par(mfrow=c(1,2),pty = "s")
  
  if (iPlot == 1) {
    #dev.new()
    win.metafile(file="GAMplot.emf",width=6,height=4)
    layout(matrix(1:6,2,3))
    #layout(matrix(1:1,1,1))
    par(mar=c(5,5,1,1))
  }

  #if (i == 6) {
  #  theXlabels = c("PDO.MJJ-4",Var.allnames[j])
  #} else {
    theXlabels = c(Var.allnames[i],Var.allnames[j])  
  #}
  
  #if (j == 127) {
  #  tempText = paste(theXlabels[2],"(day)") 
  #} else if (j == 36) {
  #  tempText = expression(paste(UWI.JAS," (m"^3,"s"^{-1},"100m"^{-1},")", sep=""))  
  #} else if (j == 61) {
  #  tempText = paste(theXlabels[2],"(mm)")
  #} else if (j == 84) {
  #  tempText = expression(paste(SST.J," ("^o,"C)", sep=""))  
  #  #tmp = theXlabels[2]
  #  #tempText = substitute(x (^o*C), list(x=tmp))
  #} else if (j == 86) {
  #  tempText = paste(theXlabels[2],"(-1000 mb)")  
  #} else if (j == 128){
  #  tempText = "log Spawners"
  #} else {
    tempText = theXlabels[2]
  #}
  
  #if (iPlot == 1) theXlabels = c("PDO.MJJ-4","log Spawners")
  theYlims=NULL
  ngraphs = 2 #c(1,2)
  yText = "Additive effect"
  if (iPlot > 2) yText = ""
  for (igraph in ngraphs) {
    #if (igraph == 2) theYlims = c(-1.0,0.8) #c(-0.8,0.8)
    plot(GAModel,select=igraph,residuals=TRUE,rug=FALSE,se=FALSE,pages=0,pch=1,cex.lab=1.3,cex.axis=1,cex=1,
      xlab=tempText,ylab=yText,ylim=theYlims) #covariate effect
	}
  
  #if (NS > 0.6) {
  if (theGCV < 0) {
    dev.new()
    plot(theYears,yobs, type="b",xlab="Year") 
    theTitle = paste("GAM: ",Var.allnames[i],"-",Var.allnames[j]," R2 = ",NS, sep="")
    title(theTitle)
    lines(theYears,ypred, col="red")
  }

} else {
#  dev.new()
#  plot(theYears,yobs, type="b",xlab="Year") 
#  theTitle = paste("GAM: ",Var.allnames[i],"-",Var.allnames[j]," R2 = ",NS, sep="")
#  title(theTitle)

  theColor = theColor + 1
  if (theColor == 7) {
    theColor = theColor + 1
  }
  thelty = 1
  if (iPlot >= 8) thelty = 2
  
  if (!jackFlag) {
    lines(theYears,ypred, col=theColor, lty=thelty)
  } else {
    #theColor = 1
    #lines(theYears,ypred, col=theColor,lty=2)
    theColor = 2
    lines(theYears,ypred, col=theColor,lty=1)
  }
  theColors = c(theColors,theColor)
  legtemp = paste(legtemp, "  R2 =", NS)
  
}
