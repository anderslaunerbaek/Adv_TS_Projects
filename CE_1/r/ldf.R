## Make a simple ldf(x,lags)
ldf <- function(x, lags, nBoot = 30, plotIt = TRUE, plotFits = FALSE, confidence_interval = 0.95) {  
  ## Calculate Lag Dependence Functions
  ## by local 2-order polynomial regression with the loess() function,
  ## and leave one out to find the best bandwidth. Finally an approximate
  ## 95% confidence interval is calculated with simple bootstrapping.
  ## Input:
  ## x, is the time series to be analysed
  ## lags, are the values for which the LDF is calculated
  ## nBoot, is the number of bootstrapping samples to make
  ## plotIt, should the ldf be plotted
  ## plotFits, should the smoothed be plotted
  
  ## The result is kept in val
  val <- vector()
  ##
  for (i in 1:length(lags)) {
    ## Take the k
    k <- lags[i]
    ## print text
    # print(paste("Calculating ldf no.", i, "of", length(lags)))
    message(paste("Calculating ldf no.", i, "of", length(lags)))
    ## Dataframe for modelling: xk is lagged k steps
    D <- data.frame(x=x[-(1:k)],xk=x[-((length(x)-k+1):length(x))])
    ## Leave one out optimization of the bandwidth with loess
    RSSk <- leaveOneOut(D, plotFits)
    ## Calculate the ldf
    RSS <- sum((D$x - mean(D$x))^2)
    val[i] <- (RSS - RSSk) / RSS
  }      
  
  ## Very simple bootstrapping
  iidVal <- vector()
  for (i in 1:nBoot) {
    ## Print to entertain the modeller ;-)
    # print(paste("Calculating bootstrap no.", i, "of", nBoot))
    message(paste("Calculating bootstrap no.", i, "of", nBoot))
    ## Bootstrapping to make a confidence band
    xr <- sample(x, min(length(x), 100), replace = TRUE)
    ## Dataframe for modelling
    DR <- data.frame(x = xr[-1],
                     xk = xr[-length(xr)])
    RSSk <- leaveOneOut(DR)
    ## The ldf is then calculated
    RSS <- sum((DR$x - mean(DR$x))^2)
    iidVal[i] <- (RSS - RSSk) / RSS
  }
  
  ## Plot the ldf
  if (plotIt) {
    dev.new()
    plot(c(0,lags), c(1,val), type="n", ylim=c(-1,1), ylab="LDF", main="Lag Dependence Functions", xaxt="n", xlab="lag")
    axis(1,c(0,lags))
    abline(0,0,col="gray")
    lines(c(0,lags), c(1,val), type="h")
    ## Draw the approximate 95% confidence interval
    abline(h=quantile(iidVal,confidence_interval), col="blue", lty=2)
  }
  
  
  
  return(data.frame(lag = c(0,lags), 
                    ldf = c(1,val), 
                    confi = rep(quantile(iidVal,confidence_interval), (length(lags) + 1)),
                    zero = 0))
}
