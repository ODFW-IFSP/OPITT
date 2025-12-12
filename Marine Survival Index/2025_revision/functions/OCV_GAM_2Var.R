#Calculate the ordinary cross-validation score

nPts = length(y)
y_fc = rep(NA,nPts)

for (iPt in 1:nPts) {

  yi = y[-iPt]
  x1i = x1[-iPt]
  x2i = x2[-iPt]  
  x= cbind(x1i,x2i)  
  dummy <- data.frame(y = yi, x = x)

  GAModel <- gam(y~s(x.x1i,k=nKnots)+s(x.x2i,k=nKnots),data=dummy)     

  yi = y[iPt]
  x1i = x1[iPt]
  x2i = x2[iPt]  
  x= cbind(x1i,x2i)  
  dummy <- data.frame(y = yi, x = x)
  
  y_fc[iPt] = predict(GAModel, dummy)
  
}

#OCV and NS Efficiency 
meany = mean(y)
ss1 = sum((y - meany)^2)
ss2 = sum((y - y_fc)^2)
theOCV = ss2 / nPts
NS_OCV = 1 - ss2/ss1
print("OCV")
print(theOCV)
theStats[iPlot,7] = NS_OCV
NS_OCV = round(NS_OCV,2)
print(NS_OCV)

#if (logFlag) {
#  y_fc = exp(y_fc)                      
#}
#
#theColor = theColor + 1
#if (theColor == 7) {
#  theColor = theColor + 1
#}
#
#theLineType = 1
#if (iPlot >= 8) theLineType = 2
#lines(theYears,y_fc, col=theColor, lty=theLineType)
#
#theColors = c(theColors,theColor)
#legtemp = paste(legtemp, "  OCV* =", NS_OCV)
