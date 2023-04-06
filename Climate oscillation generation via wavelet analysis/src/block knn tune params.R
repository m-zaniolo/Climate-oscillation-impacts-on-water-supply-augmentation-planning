##tune params for block knn algorithm
# plot
# validation test for knn generator (1)simulated power spectrum vs real  (2) boxplot for <10, 10-30 and 30-60 power 

tune_params = function(file_folder,basin_prec, basin_name,test_NS, test_NN, params, spectrum,Boxplots, Component,peaks, Bsize, Nsize){
  
  
  # wavelet analysis
  
  prep_ts = basin_prec[,basin_name]
  wt.tsN = waveletTransN(prep_ts, 0.8, nor = T) # normalize data: relative power
  
  # step1 generate simualted time series:  call function annual_generator
  

  anndf_N = annual_generator(wt.tsN, test_NS, test_NN, peaks, Bsize, Nsize, nor = T) 
  
  annSyn_df_N = anndf_N[[1]]  # annual normalize simulated timeseries
  
  annS_df_N = anndf_N[[3]] # annual normalize short osi
  
  annM_df_N = anndf_N[[4]] # annual normalize med osi
  
  annL_df_N = anndf_N[[5]] # annual normalize long osi
  
  anndf = annual_generator(wt.ts, 251, 100, peaks, Bsize, Nsize, nor = F)  # each timeseries, 251 data points. 100 timeseries in total
  anndf_N = annual_generator(wt.tsN, 251, 100, peaks, Bsize, Nsize, nor = T) 
  
  annSyn_df = anndf[[1]]  # annual simulated timeseries
  annSyn_df_N = anndf_N[[1]]  # annual normalize simulated timeseries
  
  
  annS_df = anndf[[3]] # annual short osi
  annS_df_N = anndf_N[[3]] # annual normalize short osi
  
  annM_df = anndf[[4]] # annual med osi
  annM_df_N = anndf_N[[4]] # annual normalize med osi
  
  annL_df = anndf[[5]] # annual long osi
  annL_df_N = anndf_N[[5]] # annual normalize long osi
  
  
  
  # step 2: plot simulated power vs real power spectrum
  
  if (spectrum == T){
    
    png(paste0(file_folder, basin_name, params,"_test power.png"), width = 1200, height = 500)
    par(mar = c(1,1,1,1))
    par(mfrow = c(1,2))
    test_spectrum(wt.tsN, annSyn_df_N,nor = T, xrange = c(0, 1.5), "relative")
    dev.off()
    
  }
  
  # step3: simulated power boxplot vs real power boxplot
  
  if (Boxplots == T){
    png(paste0(file_folder, basin_name, params,"_test boxplot.png"), width =550, height = 360)
    boxplot = test_boxplot(wt.tsN, annSyn_df_N, nor = T, "s")
    boxplot
    dev.off()
  }
  
  
  
}






