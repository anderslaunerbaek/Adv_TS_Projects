analyzeFit <- function(fit,tPer=NA,plotit=TRUE)
  {
    ## See the summary of the estimated parameters
    ## print(summary(fit))
    print(summary(fit,extended=TRUE))
    print(paste("loglikelihood =",-fit$f+fit$fpen))

    ##----------------------------------------------------------------
    ## Calculate the one-step predictions of the state (i.e. the residuals)
    
    pred <- predict(fit)
    
    ## Calculate the residuals and put them with the data
    data$yTiHat <- pred[[1]]$output$pred$yTi
    data$residuals <- data$yTi - data$yTiHat


    ## Cut out a period if specified
    if(is.na(tPer[1])) X <- data
    else X <- data[per(tPer[1], data$timedate, tPer[2]),]
    ##----------------------------------------------------------------
    
    if(plotit)
      {
        ##----------------------------------------------------------------
        ## Plot the auto-correlation function and cumulated periodogram in a new window
        ## New plotting device, if you use RStudio you might want not to use this
        x11()
        par(mfrow=c(1,3))
        ## The blue lines indicates the 95% confidence interval, meaning that if it is
        ##  white noise, then approximately 19 out of 20 lag correlations will be inside.
        acf(data$residuals, lag.max=8*24)
        title(main=as.character(match.call())[2],line=-2,cex.main=2)
        ## The periodogram is the estimated energy spectrum in the signal
        spec.pgram(data$residuals)
        ## The cumulated periodogram 
        cpgram(data$residuals)
        ##----------------------------------------------------------------

        ##----------------------------------------------------------------
        ## Time series plots of the inputs and residuals
        ## New plotting device, if you use RStudio you might want not to use this
        x11()
        ## Plot the residuals
        plotTSBeg(5,cex=0.7)
        gridSeq <- seq(asP("2009-01-01"),by="days",len=365)
        ##
        plot(X$timedate,X$residuals,type="n")
        abline(v=gridSeq,h=0,col="grey92")
        lines(X$timedate,X$residuals)
        title(main=as.character(match.call())[2],line=-2,cex.main=2)
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
      }

    ## Return the residuals
    invisible(data$residuals)
  }
