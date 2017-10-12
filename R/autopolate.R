devtools::use_package("fda")
devtools::use_package("data.table")

autopolate = function(dataframe, timeCol, timeFrmt, valueCol, breaksGen="normal", segmentSize= NULL, targetRate, basisRatio=0.1, smoothingAgent=0, missingIntervalSize=NULL, plot=FALSE, basis="spline",RMSE=FALSE, rate=60){

  #todo Validate input
  if( length(segmentSize) >= length(dataframe) ) stop("segment size should be smaller than dataframe lenght")
  #todo Validate input

  dt = data.table::data.table(dataframe[,c(timeCol,valueCol)])
  dt = na.omit(dt)

  x = strptime( dt[[timeCol]], timeFrmt )
  if(anyNA(x)) stop("There is a problem in the provided time format or ambiguous time entries")
  y = dt[[valueCol]]

  fncts = list()

  if(breaksGen=='gap-aware'){
    if(!is.null(segmentSize) ){
      nbSegments = floor(length(x)/segmentSize)
      remainder = length(x) %% segmentSize

      breaks = gapsBreaks(x[1:segmentSize],basisRatio,rate)
      basis = fda::create.bspline.basis(c(x[1],x[segmentSize]),breaks = breaks, norder = 4)
      penalty = fda::fdPar(basis,1,lambda = smoothingAgent)
      s = fda::smooth.basis(x[1:segmentSize],y[1:segmentSize],penalty)
      fncts[[1]] = s

      for( i in 2:nbSegments){
        breaks = gapsBreaks(x[((i-1)*segmentSize):(segmentSize*i)],basisRatio, rate)
        basis = fda::create.bspline.basis(c(x[(i-1)*segmentSize],x[segmentSize*i]),breaks = breaks, norder = 4)
        penalty = fda::fdPar(basis,1,lambda = smoothingAgent)
        s = fda::smooth.basis(x[((i-1)*segmentSize):(segmentSize*i)],y[((i-1)*segmentSize):(segmentSize*i)],penalty)
        fncts[[i]] = s
      }
      breaks = gapsBreaks(x[(nbSegments*segmentSize):(nbSegments*segmentSize+remainder)], basisRatio, rate)
      basis = fda::create.bspline.basis(c(x[nbSegments*segmentSize],x[nbSegments*segmentSize+remainder]),breaks = breaks, norder = 4)
      penalty = fda::fdPar(basis,1,lambda = smoothingAgent)
      s = fda::smooth.basis(x[(nbSegments*segmentSize):(nbSegments*segmentSize+remainder)],y[(nbSegments*segmentSize):(nbSegments*segmentSize+remainder)],penalty)
      fncts[[nbSegments+1]]=s
    }
    else{
      breaks = gapsBreaks(x,basisRatio,rate)
      basis = fda::create.bspline.basis(c( head(x,1), tail(x,1) ), breaks = breaks, norder= 4)
      penalty = fda::fdPar(basis,1,lambda = smoothingAgent)
      s = fda::smooth.basis(x,y,penalty)
      fncts[[1]] = s
    }
  }

  if(breaksGen=='normal'){
  if( !is.null(segmentSize) ){
    nbSegments = floor(length(x)/segmentSize)
    remainder = length(x) %% segmentSize

    basis = fda::create.bspline.basis(c(x[1],x[segmentSize]),nbasis =ceiling(segmentSize*basisRatio), norder = 4)
    penalty = fda::fdPar(basis,1,lambda = smoothingAgent)
    s = fda::smooth.basis(x[1:segmentSize],y[1:segmentSize],penalty)
    fncts[[1]] = s

    for( i in 2:nbSegments){
      basis = fda::create.bspline.basis(c(x[(i-1)*segmentSize],x[segmentSize*i]),nbasis =ceiling(segmentSize*basisRatio), norder = 4)
      penalty = fda::fdPar(basis,1,lambda = smoothingAgent)
      s = fda::smooth.basis(x[((i-1)*segmentSize):(segmentSize*i)],y[((i-1)*segmentSize):(segmentSize*i)],penalty)
      fncts[[i]] = s
    }

    basis = fda::create.bspline.basis(c(x[nbSegments*segmentSize],x[nbSegments*segmentSize+remainder]),nbasis =ceiling(remainder*basisRatio), norder = 4)
    penalty = fda::fdPar(basis,1,lambda = smoothingAgent)
    s = fda::smooth.basis(x[(nbSegments*segmentSize):(nbSegments*segmentSize+remainder)],y[(nbSegments*segmentSize):(nbSegments*segmentSize+remainder)],penalty)
    fncts[[nbSegments+1]]=s
  }
  else{
    basis = fda::create.bspline.basis(c( head(x,1), tail(x,1) ), nbasis = ceiling( basisRatio*length(x) ), norder= 4)
    penalty = fda::fdPar(basis,1,lambda = smoothingAgent)
    s = fda::smooth.basis(x,y,penalty)
    fncts[[1]] = s
    }
    }


  newX = seq(head(x,1),tail(x,1),targetRate)

  newY = c()
  for(v in newX){
    found=0
    for(fnct in fncts){
      if( head(fnct$argvals,1)<=v && tail(fnct$argvals,1)>=v){
        newY = c(newY,fda::eval.fd(v,fnct$fd))
        found=1
        break
      }
    }
    if(found==0) stop( cat("Couldn't find a generated smoother for value ", toString( as.POSIXct(v, origin='1970-01-01') ) ) )
  }

  newDt = data.table::data.table(newX,newY)
  names(newDt) = c(timeCol,valueCol)

  if( !is.null(missingIntervalSize) ){
  missing = missingIntervals(x,missingIntervalSize)
  missingValues = c()
  for(v in missing){
    missingValues = c(missingValues,seq(v[1]+targetRate,v[2]-targetRate,targetRate))
    }
  newDt = newDt[which(newX %in% missingValues),c('newY'):=NA,]
  }

  if(isTRUE(plot)){
    plot(x,y,col=rgb(0,0.1,1,0.65))
    for(fnct in fncts){
      denseInterval = seq(head(fnct$argvals,1),tail(fnct$argvals,1),10)
      lines(denseInterval,fda::eval.fd(denseInterval,fnct$fd))
    }
  }

  if(isTRUE(RMSE)){
    print("##################### RMSE #####################")
    all = c()
    counter=1
    for(fnct in fncts){
      RMSE = sqrt( sum( (fnct$y - fda::eval.fd(fnct$argvals,fnct$fd))^2 ) / length(fnct$argvals) )
      msg = paste("RMSE of function ",counter, " : ", RMSE)
      print(msg)
      counter = counter+1
      all = c(all,RMSE)
    }
    msg = paste("Average RMSE : ", mean(all))
  }

  return( as.data.frame(newDt) )
}

missingIntervals = function(x,threshold){
  intervals = list()
  counter = 1
  for(i in 1:(length(x)-1)){
    difference = x[i+1]-x[i]
    units(difference)= 'secs'
    if(difference>threshold){
      intervals[[counter]] = c(x[i],x[i+1])
      counter = counter+1
    }
  }
  return(intervals)
}

variationBreaks = function(x,y){
  breaks = c()
  for(i in 2:(length(y)-1)){

    #from decreasing to increasing
    if(y[i-1]>y[i] & y[i+1]>y[i]){breaks= c(breaks,x[i])}

    #from increasing to decreasing
    if(y[i-1]<y[i] & y[i+1]<y[i]){breaks= c(breaks,x[i])}

  }
  midBreaks = c()
  for(k in 1:(length(breaks)-1)){
    midBreaks = c(midBreaks,(breaks[k]+breaks[k+1])/2.0)
  }
  return(midBreaks)
}

thresholdBreaks = function(x,y,sthreshold,lthreshold){

}

gapsBreaks = function(x, basisRatio, rate ){
  missing = missingIntervals(x,rate)
  #breaks = seq(head(x,1),tail(x,1),(tail(x,1)-head(x,1))/(as.double(tail(x,1)-head(x,1),units = 'secs')*basisRatio))
  breaks = seq(head(x,1),tail(x,1),length.out = length(x)*basisRatio )
  for(interval in missing){
    breaks = breaks[breaks<=interval[1] | breaks>=interval[2] ]
  }
  return(breaks)
}


