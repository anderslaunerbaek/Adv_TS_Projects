analyzeFitEW <- function(modelAndFit,tPer=NA,plotit=TRUE)
  {
    ## Take out the model,fit, and data from the list modelAndFit
    model <- modelAndFit$model
    fit <- modelAndFit$fit
    Dat <- modelAndFit$Dat

    ## See the summary of the estimation
    print(summary(fit))
    print(summary(fit,extended=TRUE))

    ##----------------------------------------------------------------
    ## Calculate the one-step predictions of the state (i.e. the residuals)
    tmp <- model$pred(fit,Dat,opt="p",k=1)

    ## Calculate the residuals and put them with the data in a data.frame X
    Dat$data$residualsE <- Dat$data$yTiE - tmp$output$pred[,1]
    Dat$data$residualsW <- Dat$data$yTiW - tmp$output$pred[,2]
    Dat$data$yTiEHat <- tmp$output$pred[,1]
    Dat$data$yTiWHat <- tmp$output$pred[,2]

    ## Cut out a period if specified
    if(is.na(tPer[1])) X <- Dat$data
    else X <- Dat$data[per(tPer[1],Dat$data$timedate,tPer[2]),]
    ##----------------------------------------------------------------
    
    if(plotit)
      {
        ##----------------------------------------------------------------
        ## Plot the auto-correlation function and cumulated periodogram in a new window
        x11()
        par(mfrow=c(1,3))
        ## The blue lines indicates the 95% confidence interval, meaning that if it is
        ##  white noise, then approximately 19 out of 20 lag correlations will be inside.
        acf(Dat$data$residualsE, lag.max=8*24, main="East")
        title(main=as.character(match.call())[2],line=-2,cex.main=2)
        ## The periodogram is the estimated energy spectrum in the signal
        spec.pgram(Dat$data$residualsE)
        ## The cumulated periodogram 
        cpgram(Dat$data$residualsE)
        ##----------------------------------------------------------------
        ##----------------------------------------------------------------
        ## Plot the auto-correlation function and cumulated periodogram in a new window
        x11()
        par(mfrow=c(1,3))
        ## The blue lines indicates the 95% confidence interval, meaning that if it is
        ##  white noise, then approximately 19 out of 20 lag correlations will be inside.
        acf(Dat$data$residualsW, lag.max=8*24, main="West")
        title(main=as.character(match.call())[2],line=-2,cex.main=2)
        ## The periodogram is the estimated energy spectrum in the signal
        spec.pgram(Dat$data$residualsW)
        ## The cumulated periodogram 
        cpgram(Dat$data$residualsW)
        ##----------------------------------------------------------------

        ##----------------------------------------------------------------
        ## Time series plots of the inputs and residuals
        ## A new plotting window
        x11()
        ## Plot the residuals
        plotTSBeg(5)
        gridSeq <- seq(asP("2009-01-01"),by="days",len=365)
        ##
        plot(X$timedate,X$residualsE,type="n")
        abline(v=gridSeq,h=0,col="grey92")
        lines(X$timedate,X$residualsE)
        title(main=as.character(match.call())[2],line=-2,cex.main=2)
        ## 
        plot(X$timedate,X$PhE,type="n")
        abline(v=gridSeq,h=0,col="grey92")
        lines(X$timedate,X$PhE)
        ## 
        plot(X$timedate,X$yTiE,type="n")
        abline(v=gridSeq,h=0,col="grey92")
        lines(X$timedate,X$yTiE)
        ##    lines(X$timedate,X$yTiHat,col=2)
        ##    legend("bottomright",c("Measured","Predicted"),lty=1,col=1:2,bg="grey95")
        ## 
        plot(X$timedate,X$Ta,type="n")
        abline(v=gridSeq,h=0,col="grey92")
        lines(X$timedate,X$Ta)
        ##
        plot(X$timedate,X$Ps,type="n")
        abline(v=gridSeq,h=0,col="grey92")
        lines(X$timedate,X$Ps)
        ##
        plotTSXAxis(X$timedate,format="%Y-%m-%d")
        ##----------------------------------------------------------------

        ##----------------------------------------------------------------
        ## Time series plots of the inputs and residuals
        ## A new plotting window
        x11()
        ## Plot the residuals
        plotTSBeg(5)
        gridSeq <- seq(asP("2009-01-01"),by="days",len=365)
        ##
        plot(X$timedate,X$residualsW,type="n")
        abline(v=gridSeq,h=0,col="grey92")
        lines(X$timedate,X$residualsW)
        title(main=as.character(match.call())[2],line=-2,cex.main=2)
        ## 
        plot(X$timedate,X$PhW,type="n")
        abline(v=gridSeq,h=0,col="grey92")
        lines(X$timedate,X$PhW)
        ## 
        plot(X$timedate,X$yTiW,type="n")
        abline(v=gridSeq,h=0,col="grey92")
        lines(X$timedate,X$yTiW)
        ##    lines(X$timedate,X$yTiHat,col=2)
        ##    legend("bottomright",c("Measured","Predicted"),lty=1,col=1:2,bg="grey95")
        ## 
        plot(X$timedate,X$Ta,type="n")
        abline(v=gridSeq,h=0,col="grey92")
        lines(X$timedate,X$Ta)
        ##
        plot(X$timedate,X$Ps,type="n")
        abline(v=gridSeq,h=0,col="grey92")
        lines(X$timedate,X$Ps)
        ##
        plotTSXAxis(X$timedate,format="%Y-%m-%d")
        ##----------------------------------------------------------------
      }
    
    ##----------------------------------------------------------------
    ## Calculate the loglikelihood value
    print(pst("Loglikelihood ",-fit$f+fit$fpen))
    ##----------------------------------------------------------------

    invisible(Dat$data$residuals)
  }
