##----------------------------------------------------------------
## Source the scripts with functions in the "functions" folder
files <- dir("functions", full.names=TRUE)
for(i in 1:length(files)) source(files[i])

## Use the ctsmr package
library(ctsmr)
## Make a list with global parameters
prm <- list()

## latitude and longitude for the location of the house
prm$latitude <- 55.791038
prm$longitude <- 12.525545
##----------------------------------------------------------------


##----------------------------------------------------------------
## Read the X210 East room data
data <- read.csv("data.csv",stringsAsFactors=FALSE)
## Take the needed series
data <- data[,c(1,2,4,5,6,3)]
## Give names to the series
names(data) <- c("timedate","yTi","Ta","Gv","Ph","Tn")


data$timedate <- asP(data$timedate)
## timedate is the time in POSIXct, make t in hours since begining
data$t <- asHours(data$timedate-data$timedate[1])
## Make a plot of the data
## Remember that the graphical parameter 'cex' can control size of plotting text and symbols, very useful for making fine plots in reports
## use the helping function in plotTSBeg to setup plot
plotTSBeg(3,cex=0.7)
## Plot the time series
plot(data$timedate,data$Ph,type="l")
plot(data$timedate,data$Ta,type="n",ylim=c(0,45))
lines(data$timedate,data$Ta)
lines(data$timedate,data$yTi,col=2)
lines(data$timedate,data$Tn,col=3)
plot(data$timedate,data$Gv,type="l")
plotTSXAxis(data$timedate)
##----------------------------------------------------------------


##----------------------------------------------------------------
## Fit a SDE model for the heat dynamics of the East room
## Generate a new object of class ctsm
model <- ctsm()
## Add a system equation and thereby also a state
model$addSystem(dTi ~ ( 1/(Ci*Ria)*(Ta-Ti) + 1/Ci*Ph )*dt + exp(p11)*dw1)
## Set the names of the inputs
model$addInput(Ta,Gv,Ph)
## Set the observation equation: Ti is the state, yTi is the measured output
model$addObs(yTi ~ Ti)
## Set the variance of the measurement error
model$setVariance(yTi ~ exp(e11))

## Set the initial value (for the optimization) of the value of the state at the starting time point
model$setParameter(  Ti0 = c(init=15  ,lb=0     ,ub=25 ) )
## Set the initial value for the optimization
model$setParameter(  Ci = c(init=1   ,lb=1E-5  ,ub=20 ) )
model$setParameter( Ria = c(init=5   ,lb=1     ,ub=10) )
model$setParameter( p11 = c(init=1   ,lb=-30   ,ub=10 ) )
model$setParameter( e11 = c(init=-1  ,lb=-50   ,ub=10 ) )
## Run the parameter optimization
fit <- model$estimate(data)
##----------------------------------------------------------------


##----------------------------------------------------------------
## Analyze the result
## See the summary of the fit
summary(fit)
## See the extended summary, with the correlation matrix of parameter estimates
summary(fit, extended=TRUE)

## Calculate the one-step predictions of the state
pred <- predict(fit)

## Calculate the residuals and put them with the data
data$yTiHat <- pred[[1]]$output$pred$yTi
data$residuals <- data$yTi - data$yTiHat


## Cut out a period if specified, i.e. outcommend. Use the helping function per(), see functions/per.R
X <- data#[per("1983-10-11 10:00",Dat$data$timedate,"1983-10-12 18:00"),]

## Plot the auto-correlation function and cumulated periodogram in a new window
par(mfrow=c(1,3),mar=c(3,3,2,1),mgp=c(2,0.7,0),cex=0.8)
## The blue lines indicates the 95% confidence interval, meaning that if it is
##  white noise, then approximately 19 out of 20 lag correlations will be inside.
acf(X$residuals, lag.max=8*24)
title(main="One state model",line=-2,cex.main=2)
## The periodogram is the estimated energy spectrum in the signal
spec.pgram(X$residuals)
## The cumulated periodogram 
cpgram(X$residuals)

## Time series plots of the inputs and residuals
## Plot the residuals
plotTSBeg(5,cex=0.7)
gridSeq <- seq(asP("2009-01-01"),by="days",len=365)
##
plot(X$timedate,X$residuals,type="n")
abline(v=gridSeq,h=0,col="grey92")
lines(X$timedate,X$residuals)
title(main="One-step ahead residuals",line=-2)
## 
plot(X$timedate,X$Ph,type="n")
abline(v=gridSeq,h=0,col="grey92")
lines(X$timedate,X$Ph)
## 
plot(X$timedate,X$yTi,type="n")
abline(v=gridSeq,h=0,col="grey92")
lines(X$timedate,X$yTi)
lines(X$timedate,X$yTiHat,col=2)
legend("bottomright",c("Measured","Predicted"),lty=1,col=1:2,bg="grey95")
## 
plot(X$timedate,X$Ta,type="n")
abline(v=gridSeq,h=0,col="grey92")
lines(X$timedate,X$Ta)
##
plot(X$timedate,X$Gv,type="n")
abline(v=gridSeq,h=0,col="grey92")
lines(X$timedate,X$Gv)
##
plotTSXAxis(X$timedate,format="%Y-%m-%d")
##----------------------------------------------------------------


##----------------------------------------------------------------
## Calculate the loglikelihood value
paste("Loglikelihood ",-fit$f+fit$fpen)
##----------------------------------------------------------------


##----------------------------------------------------------------
## Be a bit smart and do the same in a function, see functions/Ti.R

model.Ti <- ctsm()
## Add a system equation and thereby also a state
model.Ti$addSystem(dTi ~ ( 1/(Ci*Ria)*(Ta-Ti) + 1/Ci*Ph )*dt + exp(p11)*dw1)
## Set the names of the inputs
model.Ti$addInput(Ta,Ph)
## Set the observation equation: Ti is the state, yTi is the measured output
model.Ti$addObs(yTi ~ Ti)
## Set the variance of the measurement error
model.Ti$setVariance(yTi ~ exp(e11))
## Set the initial value (for the optimization) of the value of the state at the starting time point
model.Ti$setParameter(  Ti0 = c(init=15  ,lb=0     ,ub=25 ) )
## Set the initial value for the optimization
model.Ti$setParameter(  Ci = c(init=1   ,lb=1E-5  ,ub=20 ) )
model.Ti$setParameter( Ria = c(init=20  ,lb=1     ,ub=1E4) )
model.Ti$setParameter( p11 = c(init=1   ,lb=-30   ,ub=10 ) )
model.Ti$setParameter( e11 = c(init=-1  ,lb=-50   ,ub=10 ) )
## Run the parameter optimization
fit.Ti <- model.Ti$estimate(data)

## Analyze the result
analyzeFit(fit.Ti)
##----------------------------------------------------------------


##----------------------------------------------------------------
## A two-state model implemented in functions/TiTm.R
model.TiTm <- ctsm()
## Add a system equation and thereby also a state
model.TiTm$addSystem(dTi ~ ( 1/(Ci*Rim)*(Tm-Ti) + 1/(Ci*Ria)*(Ta-Ti) + 1/Ci*Ph  )*dt + exp(p11)*dw1 )
model.TiTm$addSystem(dTm ~ ( 1/(Cm*Rim)*(Ti-Tm))*dt + exp(p22)*dw2 )
## Set the names of the inputs
model.TiTm$addInput(Ta,Ph)
## Set the observation equation: Ti is the state, yTi is the measured output
model.TiTm$addObs(yTi ~ Ti)
## Set the variance of the measurement error
model.TiTm$setVariance(yTi ~ exp(e11))
## Set the initial value (for the optimization) of the value of the state at the starting time point
model.TiTm$setParameter(  Ti = c(init=25  ,lb=0     ,ub=40) )
model.TiTm$setParameter(  Tm = c(init=25  ,lb=0     ,ub=40) )
## Set the initial value for the optimization
model.TiTm$setParameter(  Ci = c(init=1   ,lb=1E-5  ,ub=20) )
model.TiTm$setParameter(  Cm = c(init=1   ,lb=1E-4  ,ub=100))
model.TiTm$setParameter( Ria = c(init=20  ,lb=1     ,ub=1E5))
model.TiTm$setParameter( Rim = c(init=1   ,lb=1E-4  ,ub=100))
model.TiTm$setParameter( p11 = c(init=1   ,lb=-30   ,ub=10) )
model.TiTm$setParameter( p22 = c(init=1   ,lb=-30   ,ub=10) )
model.TiTm$setParameter( e11 = c(init=-1  ,lb=-50   ,ub=10) )    
## Run the parameter optimization
fit.TiTm <- model.TiTm$estimate(data)    

## Analyze the result
analyzeFit(fit.TiTm)
##----------------------------------------------------------------


##----------------------------------------------------------------
## Perform a likelihood ratio test: lambda = lik(smallerModel)/lik(largerModel) ,
## where the smallerModel is submodel of the largerModel and lambda is chi2(f)
## distributed with f=dim(smallerModel)-dim(largerModel). Page 20 in Madsen2006.
##
## Get the logLikelihood for both models from their fit
logLikSmall <- fit.Ti$loglik
logLikLarge <- fit.TiTm$loglik
## Calculate the test statistic
chisqStat <- -2 * (logLikSmall - logLikLarge)
## It this gives a p-value smaller than confidence limit, e.g. 5%, then the
## larger model is significant better than the smaller model
prmDiff <- fit.TiTm$model$NPARAM - fit.Ti$model$NPARAM
## The p-value
1 - pchisq(chisqStat, prmDiff)
##----------------------------------------------------------------
