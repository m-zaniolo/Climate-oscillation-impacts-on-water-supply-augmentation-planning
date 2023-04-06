# basic function for wavelet analysis. Adapted from Github_WaveletComp: https://github.com/cran/WaveletComp

analyze.wavelet <- function(my.data, my.series = 1, loess.span = 0.75, 
                            dt = 1, dj = 1/20, 
                            lowerPeriod = 2*dt, upperPeriod = floor(nrow(my.data)/3)*dt, 
                            make.pval = TRUE, method = "white.noise", params = NULL,
                            n.sim = 100, 
                            date.format = NULL, date.tz = NULL, 
                            verbose = TRUE, nor) {         
  
  if(verbose == T){
    out <- function(...){ cat(...) }
  }
  else{
    out <- function(...) { }
  }  
  
  ###################################################################################################
  ## The following function smoothes the series in a data frame.
  ## Input: a data frame with dates as row names
  ## Output: a data frame with the same row and column names, with smoothed series
  ###################################################################################################
  
  loess.data.frame = function(x, loess.span)  {
    x.smoothed = x
    for (i in 1:ncol(x))  {
      day.index = 1:nrow(x)
      my.loess.x = loess(x[, i] ~ day.index, span = loess.span)
      # smoothed series = fitted values:
      x.loess = as.numeric(predict(my.loess.x, data.frame(x = 1:nrow(x))))
      x.smoothed[, i] = x.loess
    }
    return(x.smoothed)
  }
  
  ###################################################################################################
  ## Select the time series to be analyzed
  ###################################################################################################
  
  if (is.numeric(my.series)) { 
    my.series = names(my.data)[my.series] 
  }
  
  if (length(my.series) != 1) { stop('Please select (only) one series for analysis!\n') }   
  if (is.element('date', my.series)) { stop('Please review your selection of series!\n') }  
  
  ind = which( names(my.data) == my.series )
  x = data.frame(my.data[,ind])   
  colnames(x) = my.series
  rownames(x) = rownames(my.data)
  
  ###################################################################################################
  ## Some initial tests
  ###################################################################################################
  
  if ( !is.numeric(x[[my.series]]) ) { stop('Some values in your time series do not seem to be interpretable as numbers.\n') }
  
  if ( sum(is.na(x[[my.series]]))>0 ) { stop('Some values in your time series seem to be missing.\n') }
  
  if ( sd(x[[my.series]]) == 0 ) { stop('Your time series seems to be constant, there is no need to search for periodicity.\n') }
  if ( lowerPeriod > upperPeriod ) { stop('Please choose lowerPeriod smaller than or (at most) equal to upperPeriod.\n') }
  
  ###################################################################################################
  ## Smooth the data (if requested)
  ###################################################################################################
  
  if (loess.span != 0) {
    out("Smoothing the time series...\n")
    x.trend = loess.data.frame(x, loess.span)
    x = x-x.trend
    x = cbind(x, x.trend)
    colnames(x) = c(my.series, paste(my.series,'.trend',sep=''))
  }
  
  ###################################################################################################
  ## Add date column if available
  ###################################################################################################  
  
  if (is.element('date',names(my.data))) {x = cbind(date = my.data$date, x)}  
  
  ###################################################################################################
  ## Start the analysis of wavelets
  ###################################################################################################
  
  out("Starting wavelet transformation...\n")
  if (make.pval == T) { out("... and simulations... \n") }
  my.wt = wt(x=x[[my.series]], start = 1, 
             dt = dt, dj = dj, 
             lowerPeriod = lowerPeriod, upperPeriod = upperPeriod,
             make.pval = make.pval,
             method = method,
             params = params,
             n.sim = n.sim, save.sim = F, nor = nor)
  
  ##################################################################################################
  ## Compute the power ridge
  ##################################################################################################
  
  Ridge = ridge(my.wt$Power)
  
  ##################################################################################################  
  ## Prepare the output  
  ##################################################################################################
  
  
  output <- list(series = x, loess.span = loess.span, dt = dt, dj = dj,
                 Wave = my.wt$Wave, Phase = my.wt$Phase, Ampl = my.wt$Ampl,
                 Power = my.wt$Power, Power.avg = my.wt$Power.avg,
                 Power.pval = my.wt$Power.pval, Power.avg.pval = my.wt$Power.avg.pval,  
                 Ridge = Ridge,     
                 Period = my.wt$Period, Scale = my.wt$Scale,                      
                 nc = my.wt$nc, nr = my.wt$nr,      
                 coi.1 = my.wt$coi.1, coi.2 = my.wt$coi.2,
                 axis.1 = my.wt$axis.1, axis.2 = my.wt$axis.2,
                 date.format = date.format, date.tz = date.tz
  )
  
  
  class(output) = "analyze.wavelet"
  
  out("Class attributes are accessible through following names:\n")
  out(names(output), "\n")
  
  return(invisible(output))
  
}

########################################################################################################################
# whether standarized timeseries or not

WaveletTransform <-function(x, dt = 1, dj = 1/20, 
                            lowerPeriod = 2*dt, upperPeriod = floor(length(x)*dt/3), nor){
  
  ###############################################################################
  ## Provide parameters (which could be useful for other transforms as well)
  ###############################################################################
  
  # Original length and length of zero padding:
  series.length = length(x)
  pot2 = trunc(log2(series.length) + 0.5)
  pad.length = 2^(pot2+1)-series.length 
  
  # Define central angular frequency omega0 and fourier factor:
  omega0 = 6
  #   fourier.factor   = (4*pi)/(omega0 + sqrt(2+omega0^2))
  fourier.factor = (2*pi)/omega0
  
  # Compute scales and periods:
  min.scale = lowerPeriod/fourier.factor             # Convert lowerPeriod to minimum scale 
  max.scale = upperPeriod/fourier.factor             # Convert upperPeriod to maximum scale 
  J = as.integer( log2(max.scale/min.scale) / dj)    # Index of maximum scale -1
  
  scales = min.scale * 2^((0:J)*dj)        # sequence of scales 
  scales.length = length(scales)           # J + 1
  periods = fourier.factor*scales          # sequence of periods
  
  # Computation of the angular frequencies
  N = series.length+pad.length
  omega.k = 1:floor(N/2)
  omega.k = omega.k * (2*pi)/(N*dt)                    # k <= N/2
  omega.k = c(0, omega.k, -omega.k[ floor((N-1)/2):1 ])
  
  ###############################################################################
  ## Define the Morlet wavelet transform function
  ###############################################################################
  
  morlet.wavelet.transform = function(x, nor) {
    
    # Standardize x and pad with zeros
    
    #############################################
    # Standardize x and pad with zeros
    
    
    if (nor == T){
      x = (x-mean(x))/sd(x) 
    }
    
    
    xpad = c(x, rep(0,pad.length)) 
    
    # Compute Fast Fourier Transform of xpad
    fft.xpad = fft(xpad)
    
    # Compute wavelet transform of x 
    # Prepare a complex matrix which accomodates the wavelet transform
    wave = matrix(0, nrow=scales.length, ncol=N)                 
    wave = wave + 1i*wave                   
    
    # Computation for each scale...
    # ... simultaneously for all time instances
    for (ind.scale in (1:scales.length)) {
      
      my.scale = scales[ind.scale]
      
      norm.factor = pi^(1/4) * sqrt(2*my.scale/dt)
      expnt       = -( (my.scale * omega.k - omega0)^2 / 2 ) * (omega.k > 0)
      daughter    = norm.factor * exp(expnt)
      daughter    = daughter * (omega.k > 0)
      
      wave[ind.scale,] = fft( fft.xpad * daughter, inverse=TRUE) / N
    }
    
    # Cut out the wavelet transform
    wave = wave[,1:series.length]
    return(wave)
  }        
  
  ###############################################################################
  ## Compute the wavelet transform, power, phases, amplitudes
  ###############################################################################
  
  Wave = morlet.wavelet.transform(x, nor) 
  
  # Compute wavelet power
  Power = Mod(Wave)^2 / matrix(rep(scales, series.length), nrow=scales.length)
  
  # Phase  
  Phase = Arg(Wave)
  
  # Amplitude
  Ampl  = Mod(Wave) / matrix(rep(sqrt(scales), series.length), nrow=scales.length)
  
  ###############################################################################
  ## Prepare the output
  ###############################################################################
  
  output = list(Wave = Wave, 
                Phase = Phase, Ampl = Ampl,
                Period = periods, Scale = scales,
                Power = Power, 
                nc = series.length, nr = scales.length)
  
  return(invisible(output))
}


###############################################################################################
wt <- function(x, start = 1, dt = 1, dj = 1/20, 
               lowerPeriod = 2*dt, upperPeriod = floor(length(x)*dt/3),
               make.pval = TRUE, method = "white.noise", params = NULL, 
               n.sim = 100, save.sim = FALSE, nor) {
  
  
  ###############################################################################
  ## Call function WaveletTransform
  ## Retrieve the wavelet transform, power, phases, amplitudes
  ###############################################################################
  
  # wavelet transform
  WT = WaveletTransform(x, dt = dt, dj = dj, 
                        lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, nor = nor)
  
  Wave  = WT$Wave 
  Phase = WT$Phase
  Ampl  = WT$Ampl
  
  Power = WT$Power
  Power.avg = rowMeans(Power)
  
  Period = WT$Period
  Scale  = WT$Scale
  nr  = WT$nr
  nc  = WT$nc
  
  rm(WT)
  
  ###############################################################################
  ## Compute p values for significance check
  ###############################################################################
  
  Power.pval = NULL
  Power.avg.pval = NULL
  series.sim = NULL
  
  if (make.pval == T) {
    
    Power.pval = matrix(0, nrow = nr, ncol = nc)
    Power.avg.pval = rep(0, nr)
    
    if (save.sim == T) { series.sim = matrix(NA, nrow=nc, ncol=n.sim) }
    
    pb = txtProgressBar(min = 0, max = n.sim, style = 3) # create a progress bar
    for(ind.sim in 1:n.sim){
      
      x.sim = SurrogateData(x, method = method)
      
      if (save.sim == T) { series.sim[,ind.sim] = x.sim }
      
      WT.sim = WaveletTransform(x.sim, dt = dt, dj = dj, 
                                lowerPeriod = lowerPeriod, upperPeriod = upperPeriod, nor)
      
      Power.sim = WT.sim$Power
      Power.avg.sim = rowMeans(Power.sim)  
      
      rm(WT.sim)
      
      Power.pval[Power.sim >= Power] = Power.pval[Power.sim >= Power] + 1
      Power.avg.pval[Power.avg.sim >= Power.avg] = Power.avg.pval[Power.avg.sim >= Power.avg] + 1
      setTxtProgressBar(pb, ind.sim) # set progress bar
    }
    close(pb) # close progress bar
    
    # p-values
    
    Power.pval = Power.pval / n.sim 
    Power.avg.pval = Power.avg.pval / n.sim  
    
  }  
  
  ###############################################################################
  ## Compute the cone of influence COI
  ###############################################################################
  
  coi = COI(start = start, dt = dt, nc = nc, nr = nr, Period = Period)
  
  ###############################################################################
  ## Prepare the output
  ###############################################################################
  
  output = list(Wave = Wave, Phase = Phase, Ampl = Ampl,
                Power = Power, Power.avg = Power.avg,
                Power.pval = Power.pval, Power.avg.pval = Power.avg.pval, 
                Period = Period, Scale = Scale,     
                coi.1 = coi$x, coi.2 = coi$y,
                nc = nc, nr = nr,    
                axis.1 = coi$axis.1, axis.2 = coi$axis.2,
                series.sim = series.sim)
  
  return(invisible(output))
}

####################################################################################################
SurrogateData <- function(x, method = "white.noise", 
                          params = list(AR = list(p = 1),
                                        ARIMA = list(p = 1, q = 1, include.mean = TRUE, sd.fac = 1, trim = FALSE, trim.prop = 0.01)
                                        #                                       ,
                                        #                                       meboot = list(trim = 0.1, force.clt = F, expand.sd = T, fiv = 5)
                          ) 
){
  
  if(method == "white.noise")  x.sur <- rnorm(length(x)) 
  if(method == "shuffle")      x.sur <- sample(x, length(x)) 
  if(method == "Fourier.rand") x.sur <- FourierRand(x) 
  
  if(method == "AR")           { 
    
    x.sur <- AR(x, params = params) 
    
  } 
  
  #   if(method == "meboot")       { 
  #   
  #      trim      = params$meboot$trim
  #      force.clt = params$meboot$force.clt
  #      expand.sd = params$meboot$expand.sd
  #      fiv       = params$meboot$fiv
  #      
  #      x.sur <- meboot(x, reps=2, trim = trim, force.clt = force.clt, expand.sd = expand.sd, fiv = fiv)$ensemble[,1]
  #      
  #   }
  
  if(method == "ARIMA")         {
    
    x.sur <- ARIMA(x, params = params)
    
  }
  
  return(invisible(x.sur))
}

###############################################################################
ridge <- function(wavelet.spectrum, band = 5, scale.factor = 0.1){
  
  min.level = scale.factor * max(wavelet.spectrum)
  
  ridge.column = function(column.vec, band=band){
    
    nrows = length(column.vec)
    
    ind = seq(1,nrows)
    band.max.vec = column.vec
    
    for (i in (1:band)) {
      
      lower.ind  = ind - i
      lower.ind[lower.ind<1] = 1
      upper.ind  = ind + i
      upper.ind[upper.ind>nrows] = nrows                 
      
      band.max.vec = pmax(band.max.vec, column.vec[lower.ind], column.vec[upper.ind])
      
    }
    
    
    my.ridge.column = rep(0,nrows)
    my.ridge.column[pmax(band.max.vec) == column.vec] = 1
    
    return(my.ridge.column)
    
  }
  
  Ridge = apply(wavelet.spectrum, 2, ridge.column, band=band)
  
  Ridge = Ridge * (wavelet.spectrum>min.level)
  
  return(invisible(Ridge))
}

############################################################################################
reconstruct <- function(WT, my.series = 1, lvl = 0, 
                        only.coi = FALSE, 
                        only.sig = TRUE, siglvl = 0.05, 
                        only.ridge = FALSE, 
                        sel.period = NULL, sel.lower = NULL, sel.upper = NULL,  
                        rescale = TRUE,
                        plot.waves = FALSE, plot.rec = TRUE, 
                        lty = 1, lwd = 1, col = 1:2, ylim = NULL,
                        show.legend = TRUE, 
                        legend.coords = "topleft", legend.horiz = FALSE, legend.text = NULL,
                        label.time.axis = TRUE, 
                        show.date = FALSE, date.format = NULL, date.tz = NULL,
                        timelab = NULL, timetck = 0.02, timetcl = 0.5,
                        spec.time.axis = list(at = NULL, labels = TRUE, 
                                              las = 1, hadj = NA, padj = NA),
                        main.waves = NULL, main.rec = NULL, main = NULL, 
                        lwd.axis = 1, 
                        verbose = TRUE) {
  
  if(verbose == T){
    out <- function(...){ cat(...) }
  }
  else{
    out <- function(...) { }
  } 
  
  #     lwd.axis = 0.25
  
  series.data = WT$series
  
  ####################################
  ## Identify the scenario
  ####################################
  
  if (class(WT) == 'analyze.wavelet') {
    
    out("Your input object class is 'analyze.wavelet'...\n") 
    
    my.series = ifelse(names(series.data)[1] == 'date', names(series.data)[2], names(series.data)[1]) 
    
    Wave  = WT$Wave
    Power = WT$Power
    Power.pval = WT$Power.pval
    Ridge = WT$Ridge     
    
  } 
  if (class(WT) == 'analyze.coherency') {   
    
    out("Your input object class is 'analyze.coherency'...\n")   
    
    if (is.numeric(my.series)) { 
      if (!is.element(my.series,c(1,2))) { stop("Please choose either series number 1 or 2!") }
      my.series = ifelse(names(series.data)[1] == 'date', names(series.data)[my.series+1], names(series.data)[my.series])  
    }
    
    ind = which( names(series.data) == my.series ) 
    which.series.num = ifelse(names(series.data)[1] == 'date', ind-1, ind)
    if (!is.element(which.series.num, c(1,2))) { stop("Your series name is not available, please check!") }
    
    if (which.series.num == 1) {
      Wave  = WT$Wave.x
      Power = WT$Power.x
      Power.pval = WT$Power.x.pval
      Ridge = WT$Ridge.x
    }
    if (which.series.num == 2) {
      Wave  = WT$Wave.y
      Power = WT$Power.y
      Power.pval = WT$Power.y.pval
      Ridge = WT$Ridge.y
    }      
    
  }   
  
  out(paste("Your time series '", my.series, "' will be reconstructed...", sep=''), '\n')   
  
  
  ####################################
  ## Prepare reconstruction components
  ####################################
  
  out("Starting the reconstruction process...\n")
  
  nc = WT$nc    
  nr = WT$nr
  
  dt = WT$dt
  dj = WT$dj
  
  Scale = WT$Scale
  Period = WT$Period   
  
  loess.span = WT$loess.span
  
  
  rec.waves = matrix(0, nrow=nr, ncol=nc)
  for (s.ind in seq_len(nr)) {
    rec.waves[s.ind,] = (Re(Wave[s.ind,])/sqrt(Scale[s.ind]))*dj*sqrt(dt)/(pi^(-1/4)*0.776)
  }
  
  # select minimum level?
  rec.waves = rec.waves * (Power >= lvl)
  comment.lvl = paste('minimum power level: ',lvl, sep='')
  
  # select ridge?
  if (only.ridge == T) {      
    rec.waves = rec.waves * Ridge  
    rec.waves[Ridge == 0] = NA
  }  
  comment.ridge = paste('only ridge: ', only.ridge, sep='')
  
  
  # use significant parts only?
  if (only.sig == T) {
    if (!is.null(Power.pval)) { 
      rec.waves = rec.waves * (Power.pval < siglvl)
      rec.waves[Power.pval >= siglvl] = NA
    }      
  }  
  if (only.sig == F) { siglvl = NA }
  comment.sig = paste('significance level: ',siglvl,sep='')
  
  # only cone of influence for reconstruction? 
  if (only.coi == T) {
    for (i in (1:nc)) {
      for (s.ind in seq_len(nr)) {
        if (Scale[s.ind] > 2^WT$coi.2[i]) {rec.waves[s.ind, i] = NA}
      }
    } 
  }  
  comment.coi = paste('only coi: ', only.coi, sep='') 
  
  # so far
  rnum.used = which(rowSums(rec.waves, na.rm=T)!=0)
  comment.periods = 'period: all relevant'
  
  
  # sel.period available?
  if (length(sel.period) != 0) {
    
    sel.rnum = numeric()
    
    # nearest available period:
    for (i in (1:length(sel.period))) {
      sel.rnum = union(sel.rnum, which(abs(Period-sel.period[i]) == min(abs(Period-sel.period[i])))) 
    }  
    
    rec.waves = rec.waves[sel.rnum,]
    comment.periods = paste('period: ',paste(as.character(round(Period[sel.rnum], 1)), collapse=', '), sep='')
    
    if (length(sel.rnum) == 1) {
      rec.waves = t(rec.waves)
    } 
    
    rnum.used = intersect(rnum.used, sel.rnum)
  } 
  
  # in case, sel.period is not available, refer to sel.upper, sel.lower
  if (length(sel.period) == 0 & ((length(sel.lower) != 0) | (length(sel.upper) != 0))) {   
    
    # in case, sel.lower (sel.upper) is not available, use minimum (maximum) period 
    if (length(sel.lower) == 0) {sel.lower = min(Period)}
    if (length(sel.upper) == 0) {sel.upper = max(Period)}             
    # in case, sel.lower > sel.upper, interchange
    if (sel.lower > sel.upper) { 
      sel.lower.h = sel.lower
      sel.lower = sel.upper
      sel.upper = sel.lower.h
    }
    sel.rnum = which(((Period >= sel.lower) & (Period <= sel.upper))) 
    rec.waves = rec.waves[sel.rnum,]
    
    # selected band / range of periods   
    sel.period.band  = Period[sel.rnum]  
    sel.period.range = as.character(round(sel.period.band, 1))
    if (length(sel.rnum) > 1) {
      sel.period.range = paste(round(range(sel.period.band), 1), collapse=' - ')
    }   
    if (length(sel.rnum) == 1) {rec.waves = t(rec.waves)}
    comment.periods = paste('period: ', sel.period.range, sep='') 
    
    rnum.used = intersect(rnum.used, sel.rnum)
    
  }
  
  
  
  ####################################
  ## Compute reconstructed series
  ####################################   
  
  # reconstructed time series
  x.r  = colSums(rec.waves, na.rm=T)
  
  # retrieve original time series
  x    = series.data[[my.series]]
  
  # rescale the reconstructed time series?
  if (rescale == T) {
    x.r  = (x.r-mean(x.r))*sd(x)/sd(x.r) + mean(x) 
  }   
  
  
  ####################################
  ## Plottings check: time axis default or not default?
  ####################################
  
  # date parameters
  
  if (is.null(date.format)) { date.format = WT$date.format }
  if (is.null(date.tz)) { date.tz = ifelse(is.null(WT$date.tz),"",WT$date.tz) }
  
  # individual time axis parameters
  
  if (!is.list(spec.time.axis))       spec.time.axis = list()
  if (is.null(spec.time.axis$at))     spec.time.axis$at = NULL
  if (is.null(spec.time.axis$labels)) spec.time.axis$labels = T
  if (is.null(spec.time.axis$las))    spec.time.axis$las = 1
  if (is.null(spec.time.axis$hadj))   spec.time.axis$hadj = NA
  if (is.null(spec.time.axis$padj))   spec.time.axis$padj = NA
  
  # initialize warning indicators for the case of individual time axis specification
  # warning: reset to time axis default?
  time.axis.warning = F 
  # warning: NAs among time axis tick marks?
  time.axis.warning.na = F
  # warning: calendar dates not chronological and/or format not standard unambiguous 
  chronology.warning = F
  
  if ( (!is.null(spec.time.axis$at)) & (label.time.axis==F) )  warning("\nPlease set label.time.axis = TRUE to make time axis specification effective.", immediate. = TRUE)
  
  ##############################
  
  # label time axis ? 
  if (label.time.axis == T) {
    
    # check if there is a user-defined time axis specification
    time.axis.default = ( is.null(spec.time.axis$at) )
    
    # check chronology and format of calendar date in case show.date = T
    if (show.date==T) { 
      
      if (is.element('date',names(series.data))) { my.date = series.data$date } else { my.date = rownames(series.data) }
      
      if (is.null(date.format)) { 
        
        chronology.warning = inherits(try(as.Date(my.date, tz=date.tz), silent=T),'try-error')
        
        if (!chronology.warning) {
          my.date = as.Date(my.date, tz=date.tz) 
          chronology.warning = ifelse( sum(is.na(my.date))> 0, TRUE, sum(diff(my.date, tz=date.tz)<0)>0 )
        }
        
      }
      if (!is.null(date.format)) {
        
        chronology.warning = inherits(try(as.POSIXct(my.date, format=date.format, tz=date.tz), silent=T),'try-error')
        
        if (!chronology.warning) {
          my.date = as.POSIXct(my.date, format=date.format, tz=date.tz) 
          chronology.warning = ifelse( sum(is.na(my.date))> 0, TRUE, sum(diff(my.date, tz=date.tz)<0)>0 )
        }
      }
      if (chronology.warning) { 
        show.date = F
        time.axis.default = T
        timelab = 'index'
      }
      
    } # end if (show.date==T)      
    
    # first: check conditions of time axis specification
    # time.axis.warning = TRUE leads to default reset
    
    if ( !time.axis.default ) { 
      
      # 1. According to show.date: Are tick marks (spec.time.axis$at) interpretable as numbers or dates?
      # show.date = F requires numeric values
      # show.date = T requires dates
      # mismatch results in default reset
      
      # case show.date==FALSE
      # Are tick marks (spec.time.axis$at) interpretable as numbers?
      
      # is.numeric returns TRUE iff there is no character value but at least one non-NA value which is interpretable as number
      # (further logical entries T and F are interpreted as 1 and 0 which could be accepted here, non-positive values as well,
      # since values are matched with the true range of values in the axis command)
      
      if (show.date==F) { time.axis.warning = (!is.numeric(spec.time.axis$at)) }
      
      # case show.date==TRUE
      # Are tick marks (spec.time.axis$at) interpretable as dates using given date format and time zone specifications?
      
      # checking the first value given in spec.time.axis$at:
      # inherits is TRUE (and thus gives a time axis warning) iff the first value cannot be formatted, regardless of subsequent values
      # date formatting may produce NAs, however, therefore checking is necessary whether there are NAs ONLY
      # (NA tick marks are omitted automatically by the axis command as long as there are non-NAs)
      
      if (show.date==T) { 
        
        if (is.null(date.format)) { 
          time.axis.warning = inherits(try(as.Date(spec.time.axis$at, tz=date.tz), silent=T), 'try-error') 
          if (!time.axis.warning) { time.axis.warning = (sum(!is.na(as.Date(spec.time.axis$at, tz=date.tz)))==0) }
        }
        
        if (!is.null(date.format)) { 
          time.axis.warning = inherits(try(as.POSIXct(spec.time.axis$at, format=date.format, tz=date.tz), silent=T),'try-error') 
          if (!time.axis.warning) { time.axis.warning =  (sum(!is.na(as.POSIXct(spec.time.axis$at, format=date.format, tz=date.tz)))==0) }
        }            
        
      } # end if (show.date==T)
      
      if (!time.axis.warning) { 
        # 2. Are labels (spec.time.axis$labels) appropriate? Either logical and single-valued or non-logical and of the same length as tick marks?
        # is.logical returns TRUE iff there is not any character or numeric value, but there could be NA, even if it is the only value.
        
        # if logical:  single-valued and non-NA? 
        if ( is.logical(spec.time.axis$labels) ) { 
          time.axis.warning = ( ifelse(length(spec.time.axis$labels)!=1,TRUE, is.na(spec.time.axis$labels)) )
        }
        
        # if non-logical: do vectors of tick marks and labels have equal length?
        if (!is.logical(spec.time.axis$labels)) {
          time.axis.warning = ( length(spec.time.axis$labels) != length(spec.time.axis$at) ) 
        }
      } # end if (!time.axis.warning)
      
    } # end if ( !time.axis.default )
    
    # default reset in case of a warning
    time.axis.default = (time.axis.default | time.axis.warning)
    
    if ( (is.null(timelab) & time.axis.default) | (!is.null(timelab) & time.axis.warning) ) {timelab=ifelse(show.date, 'calendar date','index')}
    
  } # end if (label.time.axis == T)
  
  ####################################
  ## Plottings
  ####################################
  
  ## plot of reconstruction waves
  
  if (plot.waves == T) {
    
    out("Reconstruction waves are being plotted...\n")
    
    if (is.null(main) == F) {main.waves = main}
    
    range.rec.waves = range(rec.waves, na.rm =T)
    #         matplot(WT$axis.1, t(rec.waves), type = 'l', ylim = range.rec.waves,
    matplot(1:nc, t(rec.waves), type = 'l', ylim = range.rec.waves,
            main = main.waves, sub = paste(comment.lvl, ', ', comment.sig,', ',comment.coi,', ',comment.ridge,', ', comment.periods,  sep=''),
            xaxs = 'i', xaxt = 'n', 
            xlab = '', ylab = '')
    
    # label time axis ?   
    if (label.time.axis == T) {
      
      #####################
      
      if (show.date == F) {
        
        if (time.axis.default) {
          
          A.1 = axis(1, lwd = lwd.axis, labels=NA, tck = timetck, tcl = timetcl)
          mtext(A.1, side = 1, at = A.1, line = par()$mgp[2]-0.5, font = par()$font.axis, cex=par()$cex.axis)
          
        } # end if (time.axis.default)
        
        if (!time.axis.default) {
          
          time.tick = spec.time.axis$at
          time.tick.label = spec.time.axis$labels
          
          # NAs among tick marks?
          which.na = which(is.na(time.tick))
          # NA warning
          if ( length(which.na)>0 ) { time.axis.warning.na = T } 
          
          # in case of NAs among tick marks:
          # NA tick marks (and corresponding labels) are omitted automatically by the axis command
          
          axis(1, lwd = lwd.axis, at=time.tick, labels=time.tick.label, tck = timetck, tcl = timetcl,
               las = spec.time.axis$las, hadj = spec.time.axis$hadj, padj=spec.time.axis$padj,
               mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
          
        } # end if (!time.axis.default)
        
        mtext(timelab, side = 1, line = par()$mgp[1]-1, font = par()$font.lab, cex=par()$cex.lab)
        
      } # end if (show.date == F)
      
      #####################
      
      if (show.date == T) {  
        
        par(new=TRUE)
        
        if (time.axis.default) {
          
          # empty plot, but calendar axis (w/ plotting)
          plot(my.date, seq(range.rec.waves[1], range.rec.waves[2],length.out=nc), 
               ylim = range.rec.waves,
               type="n", 
               xaxs = "i", 
               yaxt='n', 
               xlab="", ylab="",
               lwd = lwd.axis, tck=timetck, tcl=timetcl, 
               mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
          
        } # end if (time.axis.default)
        
        if (!time.axis.default) {
          
          # empty plot, but calendar axis (w/o plotting)
          plot(my.date, seq(range.rec.waves[1], range.rec.waves[2],length.out=nc), 
               ylim = range.rec.waves,
               type="n", 
               xaxs = "i", 
               xaxt='n', yaxt='n', 
               xlab="", ylab="")
          
          # user-defined calendar axis
          if (is.null(date.format)) { time.tick = as.Date(spec.time.axis$at, tz=date.tz) }
          if (!is.null(date.format)) { time.tick = as.POSIXct(spec.time.axis$at, format=date.format, tz=date.tz) }
          
          time.tick.label = spec.time.axis$labels
          
          # NAs among tick marks?
          which.na = which(is.na(time.tick))
          # NA warning
          if ( length(which.na)>0 ) { time.axis.warning.na = T } 
          
          # in case of NAs among tick marks:
          # NA tick marks (and corresponding labels) are omitted automatically by the axis command
          
          if ( is.logical(time.tick.label) ) {
            # then time.tick.label has length 1 by criterion and IS TRUE OR FALSE
            if (time.tick.label == T) {time.tick.label = time.tick}
          }   
          
          axis(1, lwd = lwd.axis, at=time.tick, labels=time.tick.label, tck = timetck, tcl = timetcl,
               las = spec.time.axis$las, hadj = spec.time.axis$hadj, padj=spec.time.axis$padj,
               mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
          
        } # end if (!time.axis.default)             
        
        mtext(timelab, side = 1, line = par()$mgp[1]-1, font = par()$font.lab, cex=par()$cex.lab)  
        
      }  # end if (show.date == T)   
      
    } # end if (label.time.axis == T)                  
    
    
  } #end if (plot.waves == T)
  
  x.rec = cbind(series.data, x.r=x.r)
  my.rec.series = paste(my.series,'.r',sep='')
  colnames(x.rec) = c(colnames(series.data),my.rec.series)
  rownames(x.rec) = rownames(series.data)
  
  if (plot.waves & plot.rec) { 
    par(ask = T)
  }
  
  
  
  ## plot of reconstructed series
  
  if (plot.rec == T) {
    
    out("Original (detrended) and reconstructed series are being plotted...\n")
    
    if (is.null(main) == F) { main.rec = main }
    
    range.rec = range(x.rec[,c(my.series,my.rec.series)], na.rm=T)
    if (is.null(ylim)) {ylim = range.rec} 
    #        matplot(WT$axis.1, x.rec[,c(my.series,my.rec.series)], type = 'l', ylim = ylim, lty = lty, lwd = lwd, col = col,
    matplot(1:nc, x.rec[,c(my.series,my.rec.series)], type = 'l', ylim = ylim, lty = lty, lwd = lwd, col = col,
            main = main.rec, sub = paste(comment.lvl, ', ', comment.sig,', ',comment.coi,', ',comment.ridge,', ', comment.periods,  sep=''),
            xaxs = 'i', xaxt = 'n', 
            xlab = '', ylab = '')
    
    par(ask = F)
    
    if (show.legend == T) {      
      
      if (is.null(legend.text)) { 
        legend.text = c(paste('original', ifelse(loess.span!=0, paste(' (detrended, span: ', loess.span, ')', sep=''),''), sep=''), 'reconstructed')
      }
      
      legend(legend.coords, horiz=legend.horiz, legend = legend.text, lty = lty, lwd = lwd, col = col, text.font = par()$font.lab, cex=par()$cex.lab)        
    }
    
    
    # label time axis ?   
    if (label.time.axis == T) {
      
      #####################
      
      if (show.date == F) {
        
        if (time.axis.default) {
          
          A.1 = axis(1, lwd = lwd.axis, labels=NA, tck = timetck, tcl = timetcl)
          mtext(A.1, side = 1, at = A.1, line = par()$mgp[2]-0.5, font = par()$font.axis, cex=par()$cex.axis)
          
        } # end if (time.axis.default)
        
        if (!time.axis.default) {
          
          time.tick = spec.time.axis$at
          time.tick.label = spec.time.axis$labels
          
          # NAs among tick marks?
          which.na = which(is.na(time.tick))
          # NA warning
          if ( length(which.na)>0 ) { time.axis.warning.na = T } 
          
          # in case of NAs among tick marks:
          # NA tick marks (and corresponding labels) are omitted automatically by the axis command
          
          axis(1, lwd = lwd.axis, at=time.tick, labels=time.tick.label, tck = timetck, tcl = timetcl,
               las = spec.time.axis$las, hadj = spec.time.axis$hadj, padj=spec.time.axis$padj,
               mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
          
        } # end if (!time.axis.default)
        
        mtext(timelab, side = 1, line = par()$mgp[1]-1, font = par()$font.lab, cex=par()$cex.lab)
        
      } # end if (show.date == F)
      
      #####################
      
      if (show.date == T) {  
        
        par(new=TRUE)
        
        if (time.axis.default) {
          
          # empty plot, but calendar axis (w/ plotting)
          plot(my.date, seq(ylim[1], ylim[2], length.out=nc), 
               ylim = ylim,
               type="n", 
               xaxs = "i", 
               yaxt='n', 
               xlab="", ylab="",
               lwd = lwd.axis, tck=timetck, tcl=timetcl, 
               mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)                         
          
        } # end if (time.axis.default)
        
        if (!time.axis.default) {
          
          # empty plot, but calendar axis (w/o plotting)
          plot(my.date, seq(ylim[1], ylim[2], length.out=nc), 
               ylim = ylim,
               type="n", 
               xaxs = "i", 
               xaxt='n', yaxt='n', 
               xlab="", ylab="")        
          
          # user-defined calendar axis
          if (is.null(date.format)) { time.tick = as.Date(spec.time.axis$at, tz=date.tz) }
          if (!is.null(date.format)) { time.tick = as.POSIXct(spec.time.axis$at, format=date.format, tz=date.tz) }
          
          time.tick.label = spec.time.axis$labels
          
          # NAs among tick marks?
          which.na = which(is.na(time.tick))
          # NA warning
          if ( length(which.na)>0 ) { time.axis.warning.na = T } 
          
          # in case of NAs among tick marks:
          # NA tick marks (and corresponding labels) are omitted automatically by the axis command
          
          if ( is.logical(time.tick.label) ) {
            # then time.tick.label has length 1 by criterion and IS TRUE OR FALSE
            if (time.tick.label == T) {time.tick.label = time.tick}
          }   
          
          axis(1, lwd = lwd.axis, at=time.tick, labels=time.tick.label, tck = timetck, tcl = timetcl,
               las = spec.time.axis$las, hadj = spec.time.axis$hadj, padj=spec.time.axis$padj,
               mgp = par()$mgp - c(0,0.5,0), font = par()$font.axis, cex.axis=par()$cex.axis)
          
        } # end if (!time.axis.default)             
        
        mtext(timelab, side = 1, line = par()$mgp[1]-1, font = par()$font.lab, cex=par()$cex.lab)  
        
      }  # end if (show.date == T)  
      
    } # end if (label.time.axis == T)                            
    
  } # if (plot.rec == T)
  
  if (chronology.warning) { 
    warning("\nPlease check your calendar dates, format and time zone: dates may not be in an unambiguous format or chronological. The default numerical axis was used instead.") 
  }  
  
  if (time.axis.warning == T) {
    warning("\nPlease check your time axis specifications. Default settings were used.")
  }
  if (time.axis.warning.na == T) {
    warning("\nNAs were produced with your time axis specifications.")
  }
  
  ####################################
  ## Output
  ####################################
  
  output <- list(series = x.rec,
                 rec.waves = rec.waves,
                 loess.span = loess.span,
                 lvl = lvl,
                 only.coi = only.coi,
                 only.sig = only.sig, siglvl = siglvl, 
                 only.ridge = only.ridge,                
                 rnum.used = rnum.used,   
                 rescale = rescale,
                 dt = dt, dj = dj,
                 Period = Period, Scale = Scale,
                 nc = nc, nr = nr, 
                 axis.1 = WT$axis.1,
                 axis.2 = WT$axis.2,
                 date.format = date.format, date.tz = date.tz
  )
  
  class(output) = "reconstruct"
  
  out("Class attributes are accessible through following names:\n")
  out(names(output), "\n")             
  
  return(invisible(output))
  
  
}

########################################################################
periodic.series = function(start.period = 100, end.period = start.period, phase = 0, length = 600, make.plot = FALSE){
  
  if (start.period < end.period) { m = start.period/end.period }
  if (start.period >= end.period) { m = end.period/start.period }
  N = length/(1-0.5*(1-m))
  x = 0:(N-1)
  if (start.period < end.period) { period = end.period }
  if (start.period >= end.period) { period = start.period }
  p = rep(period, length = N)
  y = sin((x + phase) * 2 * pi/p)
  if (start.period < end.period) { y = y[N:1] }
  # "thinning out":
  i = 1:N
  cum.keep.share = 1 - 0.5*i*(1-m)/N
  keep = floor(cum.keep.share*i)
  diff.keep = c(1, diff(keep))
  if (sum(diff.keep) < length) {diff.keep[match(0, diff.keep)] = 1}
  if (sum(diff.keep) > length) {diff.keep[match(1, diff.keep)] = 0}
  if (start.period < end.period) {
    # modify diff.keep:
    first.1.rev = match(1,diff.keep[N:1]) # first one in reversed diff.keep
    if (first.1.rev > 1)  {
      # cut off 0s at the end, add at beginning:
      diff.keep = c(rep(0, first.1.rev - 1), diff.keep[1:(N-(first.1.rev - 1))])
    }
  }
  if (start.period >= end.period) {
    # modify diff.keep:
    first.1 = match(1,diff.keep) # first one in diff.keep
    if (first.1 > 1)  {
      # cut off 0s at the beginning, add at the end:
      diff.keep = c(diff.keep[first.1:N], rep(0, first.1 - 1))
    }
  }
  z = y[diff.keep == 1]
  if (start.period < end.period) {z = rev(z)}
  if (make.plot) {ts.plot(z, xlab = 'time', ylab = 'periodic series')}
  return(z)
}


### EXAMPLE:

# x1 = 0.8*periodic.series(start.period = 100, end.period = 95, phase = 0, length = 1000)
# x2 =     periodic.series(start.period = 100, end.period = 100, phase = 0, length = 1000)
# x3 = 1.2*periodic.series(start.period = 100, end.period = 105, phase = 0, length = 1000)
# 
# ts.plot(x2, ylim = c(-2, +2), xlab = 'time', ylab = 'series with variable period')
# lines(x1, col = 'blue')
# lines(x3, col = 'red')
# legend('topleft', legend = c('speeding up (end period = 95)', 'period = 100', 'slowing down (end period = 105)'), lty = 1, col = c('blue', 'black', 'red'))