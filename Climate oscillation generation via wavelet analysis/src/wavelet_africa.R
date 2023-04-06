############################################################
#####################
# compute relative power contributed by short, med, long period range, for all sub basins across Africa: paper fig2

Africa_relative_power = function(basin_prec,file_name){
  
  basin_num = c(1:6, 8:20, 22:28, 30:83, 85:99)
  basin_index = paste(rep("basin_", 95), basin_num, sep = "")
  
  n = length(basin_index)
  dflist = list()
  
  # loop over all the basins prec, compute relative power for each sub basins
  
  for (i in 1 : n){
    power_df = relative_power(basin_prep, basin_index[i], ratio = F, q_mat_n)
    dflist[[i]] = power_df
  }
  
  power_total = do.call(rbind, dflist)
  write.csv(power_total, file = paste0(file_name, "relative_power_africa.csv"))
  
}

# function: for a given sub_basin, compute relative power contributed by short, med, long period range
# input: basin_prec ; basin_name; 
# output: power_df (dataframe) : three columns (period, basin_name, value)


relative_power = function(basin_prep, basin_name){
  
  
  # extract prep
  prep_ts = basin_prep[, basin_name]
  
  # wavelet transform
  wt_ts = waveletTransN(prep_ts, 0.8, nor = T) # normalize data: relative power
  
  # compute area under curve
  period_p = compute_area(wt_ts)

  total_p = sum(period_p)
  period_p = period_p / total_p   # percent of power contributed from each period range
  #power_df = data.frame("period" = c("noise","short", "medium", "long"), "name" = rep(basin_name, 4),  value = period_p)
  power_df = data.frame("period" = c("short", "medium", "long"), "name" = rep(basin_name, 3),  value = period_p)
    
  return (power_df)
  
}

# helper function 
# compute_area(wt) : compute area under curve for a power spectrum series
# input :wt: dataframe. wt$Period, wt$Power.avg
# output: period_p: vector of pwoer contributed from (<10, 10-30, >30)

compute_area = function(wt){
  
  ts_p1 = trapz(wt$Period[wt$Period <= 10], wt$Power.avg[ wt$Period <= 10])
  ts_p2 = trapz(wt$Period[wt$Period > 10 & wt$Period <= 30], wt$Power.avg[wt$Period > 10 & wt$Period <= 30])
  ts_p3 = trapz(wt$Period[wt$Period > 30 & wt$Period <= 60], wt$Power.avg[wt$Period > 30 & wt$Period <= 60])

  period_p = c(ts_p1, ts_p2, ts_p3)
  return(period_p)
  
}



#############################################################
# function: add noise spectrum on average power
# plot_withN(x,  nsim, loess)
# input: x: time series ;   nsim: # of simulations ; loess: detrend factor
# return:simulated quantile noise power matrix: q_mat

plot_withN = function(x,  nsim, loess){
  
  # 1. wavelet transform on time series
  
  # detrend time series
  x.detrend = detrend(x, loess)$detrend
  
  
  # wavelet transform time series
  
  wt.x =  WaveletTransform(x.detrend, dt =1, dj = 1/100, lowerPeriod = 2, upperPeriod = floor(length(x.detrend)*1/3), nor = T)
  
  periods = wt.x$Period
  Power.ave = rowMeans(wt.x$Power)
  
  
  # 2. significance: simulate 100 white noise
  
  n = length(x.detrend)
  nr = wt.x$nr
  
  sim_matrix = matrix(0, nr, nsim) # nr: periods     # matrix for noise average power
  Power.avg.pval = rep(0, nr)
  
  set.seed(6)
  for (i in 1:nsim){
    # simulate white noise
    
    x.sur = rnorm(n)
    
    # wavelet transform on white noise
    wt.sur =  WaveletTransform(x.sur, dt =1, dj = 1/100, lowerPeriod = 2, upperPeriod = floor(length(x.sur)*1/3), nor = T)
    
    # average power of simulated white noise
    power.avg.sur = rowMeans(wt.sur$Power)
    sim_matrix[,i] = power.avg.sur
    
    # pvalue vector
    
    Power.avg.pval[power.avg.sur >= Power.ave] = Power.avg.pval[power.avg.sur >= Power.ave] + 1
  }
  
  # compute pvalues
  Power.avg.pval = Power.avg.pval / nsim 
  
  # 3. given sim_matrix(100 cols, each col is noise power spectrum)
  #compute percentile for each row
  
  q_mat = matrix(0, nr,5)
  prob = c(0.25, 0.5, 0.75, 0.9)
  for ( i in 1:4){
    q_mat[,i] = apply(sim_matrix, 1, quantile, probs = prob[i], na.rm = TRUE)
  }
  q_mat[,5] = periods
  
  # 4. plot
  make_PlotWithNoise(wt.x, q_mat, Power.ave, Power.avg.pval)
  
  return(q_mat)
  
}


# helper n: given wt.x, noise percentile matrix q_mat, avg.pvalue
#make_PLotWithNoise(WT, q_mat, Power.ave, Power.avg.pval)
# large periods at the top
make_PlotWithNoise = function(WT, q_mat,  Power.ave, Power.avg.pval){
  
  Power.avg = Power.ave
  minimum.level = min(Power.avg)
  maximum.level = max(Power.avg)
  maximum.level = max(q_mat[,4])
  xlim = range(c(Power.avg,minimum.level,maximum.level))
  min.p = min(log2(WT$Period))
  max.p = max(log2(WT$Period))
  ylim = range(c(log2(WT$Period), min.p, max.p))
  
  par(cex.axis=1.8, cex.lab=1.5)
  plot(Power.avg, log2(WT$Period), xlim = xlim, 
       ylim = ylim,
       lwd = 2.5, col = 1, type = "l", 
       axes = FALSE, 
       ylab = "", xlab = '', 
       yaxs = 'i', 
       #main = 'Average wavelet power spectrum for timeseries')
       main = '')
  title(ylab="Period(years)",xlab = 'Global Wavelet Power(GWP)', mgp=c(2.5,1.5,0),cex.lab=1.8)
  
  # add sig points
  siglvl = 0.1
  P.dat = data.frame(Pvalue = Power.avg.pval, Log.period = log2(WT$Period), Average = Power.avg)
  
  with(P.dat[P.dat$Pvalue < 0.1,], points(Average, Log.period, pch = 20, col = "red", cex = 1))
  
  # add sig periods
  sig_P = P.dat[P.dat$Pvalue < 0.1,]
  sig_peaks_period = sig_P$Log.period[find_peak(sig_P$Average)]
  sig_peaks_power = sig_P$Average[find_peak(sig_P$Average)]
  if (length(sig_peaks_period)>0){
    text(sig_peaks_power, sig_peaks_period,round(2^sig_peaks_period), cex=1.5, pos=4, col="red")
  }
  
  
  # add noise spectrum
  #library(RColorBrewer)
  noise_color = c( "paleturquoise","steelblue1", "royalblue1", "royalblue3")  #"lightskyblue",
  for (i in 1:4){
    if (i == 2 || i==4){
      lines(q_mat[,i],log2(WT$Period), col = noise_color[i], lwd = 2, lty=5)     #brewer.pal(8, "Blues")[2*i]
    }
  }
  
  # show legend
  #legend("topright", legend=c("1% significance", "precipitation spectrum", "50% noise spectrum", "90% noise spectrum"), pch=c(20, NA,NA,NA),lty = c(NA, 1,2,2), col=c("red","black",noise_color[2],noise_color[4]), horiz=F, box.lwd=1, text.font = par()$font.lab, cex=par()$cex.lab, pt.cex=1)
  
  # outside box
  box(lwd = 1)
  
  # set x labels
  A.1 = axis(1, lwd = 1, labels=NA, tck = 0.02, tcl = 0.5)
  mtext(A.1, side = 1, at = A.1, line = par()$mgp[2]-0.5, font = par()$font.axis, cex=par()$cex.axis)
  #mtext('average wavelet power', side = 1, line = par()$mgp[1]-1, font = par()$font.lab, cex=par()$cex.lab) 
  
  
  # label period (y) axis
  period.tick <- unique(trunc(log2(WT$Period)))
  period.tick[period.tick<log2(WT$Period[1])] = NA
  
  period.tick = na.omit(period.tick)
  period.tick.label <- 2^(period.tick) 
  
  axis(2, lwd = 1, at = period.tick, labels = NA, tck=0.02, tcl=0.5)
  axis(4, lwd = 1, at = period.tick, labels = NA, tck=0.02, tcl=0.5)
  mtext(period.tick.label, side = 2, at = period.tick, las = 1, line = par()$mgp[2]-0.5, font = par()$font.axis, cex=par(cex.axis = 1.8)$cex.axis)
  
  #mtext("period (years)", side = 2, line = par()$mgp[1]-0.5, font = par()$font.lab, cex=par()$cex.lab)
  
}

# helper n: detrend time series
# detrend(x, loess)
detrend = function(x, loess){
  
  index = 1:length(x)
  my.loess.x = loess(x ~ index, span = loess)
  # smoothed series = fitted values:
  x.loess = as.numeric(predict(my.loess.x, data.frame(index)))
  
  x.detrend = x - x.loess
  detrend_df = data.frame('trend_line' = x.loess, 'detrend' = x.detrend )
  
  return (detrend_df)
  
}

# helper n : given timeseries, find local peak. Return local peak index
# find_peak function
find_peak = function (x, thresh = 0) {
  
  pks <- which(diff(sign(diff(x, na.pad = FALSE)), na.pad = FALSE) < 0) + 1
  if (!missing(thresh)) {
    pks[x[pks - 1] - x[pks] > thresh]
  }
  else pks
}
