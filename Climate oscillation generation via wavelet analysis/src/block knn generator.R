# function for block knn generator and tune params
# this used for generating annual prep series / annual oscillation series and monthly timeseries for rainfall-runoff model

source("block knn helper.R")
source("WaveletTransform_github.R")
source("wt_github.R")
source("analyzeWavelet_github.R")
source("data processing.R")
source("SurrogateData_github.R")
source("ridge_github.R")
source("COI_github.R")
source("reconstruct_github.R")



##########################
annual_generator = function(basin_name, file_name,Short, Med, Long,S1,Resi, NS, NN, peaks, Bsize, Nsize, nor = T){
  
  prep_ts = basin_prep[,basin_name]
  wt.tsN = waveletTransN(prep_ts, 0.8, nor = T)   # normalize and detrend data before wavelet transform
  
  
  #call  function annual_generator
  
  anndf_N = annual_generator(wt.tsN, NS, NN, peaks, Bsize, Nsize, nor = T)    # in normalized space
  annSyn_df_N = anndf_N[[1]]  # annual simulated timeseries
  annS1_df_N = anndf_N[[2]] 
  annS_df_N = anndf_N[[3]] # annual short osi
  annM_df_N = anndf_N[[4]] # annual med osi
  annL_df_N = anndf_N[[5]] # annual long osi
  annResi_df_N = anndf_N[[6]]
  
  #save annual osi file (NS * NN. Each column vector is one simulated timeseries)
  
  if (Short == T){
    write.csv(annS_df_N, paste0(file_name, "osi_", basin_code, "_s.csv"))
  }
  if (Med == T){
    write.csv(annM_df_N, paste0(file_name, "osi_", basin_code, "_m.csv"))
  } 
  if (Long == T){
    write.csv(annL_df_N, paste0(file_name, "osi_", basin_code, "_l.csv"))
  }
  if (S1 == T){
    write.csv(annS1_df_N, paste0(file_name, "osi_", basin_code, "_s1.csv"))
  }
  if (Resi == T){
    write.csv(annResi_df_N, paste0(file_name, "osi_", basin_code, "_resi.csv"))
  }
  # test simulated power
  png(paste0(file_name, basin_name, "_test boxplot.png"), width = 800, height = 350)
  test_boxplot(wt.tsN, annSyn_df_N, nor = T)
  dev.off()
  
  return(annSyn_df_N)
}

#####################################
# output monthly timeseries in csv
month_generator = function(basin, annSyn_df_N, file_name){
  
  raw_prec = read_CMIP6(read_cmip = T)
  pre_his = raw_prec[[1]]
  pre_fut = raw_prec[[2]]
  prec_annual = compute_yearly(pre_his,pre_fut,basin)  

  wt_basin = waveletTransN(prec_annual, 0.8, nor = T)

  #loop over NN, generate NN rows of monthly simulations
  # for each row, one simulation of NS-year monthly (NS*12) data
  
  monthSyn_df = data.frame(matrix(0, nrow = NN, ncol = 100 *12 ))
  
  set.seed(6)
  for (i in 1:NN){
    
    annPsyn = annSyn_df_N[,i]
    monthPsyn = monthly_knn(wt_basin, pre_his, pre_fut, basin, annPsyn, Trend = T)
    monthP_vec = as.vector(t(as.matrix(monthPsyn))) # reshape into a vector
    
    monthSyn_df[i, ] = monthP_vec
    write.csv(monthSyn_df, paste0(file_name, "rainfall_", basin_code, ".csv"))
    
  }
  
}

#############################
# functions for tune params

# a. simulated power vs real power spectrum
test_spectrum <- function(wt, annSyn, nor = F, xrange, title){
  
  ave_p = wt$Power.avg
  period_t = wt$Period
  
  #png(paste0(file_folder, basin_name, paste0("/", params,"_test power.png")), width = 600, height = 400)
  
  plot(ave_p, log2(period_t), type = 'l', col = 'black', lwd = 2, main = title,  xlab = 'average_power', ylab = 'period', xlim = xrange)
  
  x.tick <- unique(trunc(log2(period_t)))
  x.tick.label = 2^x.tick
  axis(1, at = x.tick, labels = x.tick.label)
  
  
  for (i in 20:29){
    wt.sim = waveletTransN(annSyn[,i], 0, nor = F)
    lines(wt.sim$Power.avg, log2(period_t), col = brewer.pal(10,"Set3")[i-19], lwd = 2)
    
  }
  
}


# b. simulated power boxplot vs real power boxplot

test_boxplot <- function(wt, annSyn, nor = F, basins){
  
  # 100 simulated timeseries's power matrix  100*3
  spower_mat = matrix(0, nrow = 100, ncol = 3)
  
  for (i in 1:100){
    
    wt.sim = waveletTransN(annSyn[,i],0, nor = nor)
    spower_mat[i, 1] = trapz(wt.sim$Period[wt.sim$Period <= 10], wt.sim$Power.avg[wt.sim$Period <= 10] )
    spower_mat[i, 2] = trapz(wt.sim$Period[wt.sim$Period > 10 & wt.sim$Period <= 30], wt.sim$Power.avg[wt.sim$Period > 10 & wt.sim$Period <= 30])
    spower_mat[i, 3] = trapz(wt.sim$Period[wt.sim$Period > 30 & wt.sim$Period <= 60], wt.sim$Power.avg[wt.sim$Period > 30 & wt.sim$Period <= 60])
    
  }
  
  # real timeseries' power 
  
  ts_p1 = trapz(wt$Period[wt$Period <= 10],wt$Power.avg[wt$Period <= 10])
  #ts_p1 = trapz(wt$Period[wt$Period > 4 & wt$Period <= 10], wt$Power.avg[wt$Period > 4 & wt$Period <= 10])
  
  ts_p2 = trapz(wt$Period[wt$Period > 10 & wt$Period <= 30], wt$Power.avg[wt$Period > 10 & wt$Period <= 30])
  ts_p3 = trapz(wt$Period[wt$Period > 30 & wt$Period <= 60], wt$Power.avg[wt$Period > 30 & wt$Period <= 60])
  
  names = c("<= 10 years", "10-30 years", "30-60 years")
  real_p = data.frame('block' =  names, "power" = c(ts_p1, ts_p2, ts_p3))
  
  # boxplot
  # prepare data for box plot
  spower_df = data.frame(spower_mat)
  id = c(1:100)
  spower_df[,4] = id
  colnames(spower_df) = c(names, "id")
  # melt dataframe
  melt_spower =  reshape2:: melt(spower_df, id = "id")
  
  
  boxplot <- ggplot(data = NULL) + 
    geom_boxplot(data = melt_spower, aes(x = factor(variable), y = value), show.legend = FALSE)+
    geom_point(data = real_p, aes(x = factor(block), y = power, col ="red"), size = 3,show.legend = FALSE)+
    scale_color_manual(values = "red",label = "historical power")+
    coord_cartesian(ylim=c(0,3))+
    labs(title= basins ,x="", y = "Relative Power") +
    theme(plot.title = element_text(hjust = 0.5)) +
    my_theme +
    theme(legend.position="bottom", legend.box = "horizontal")
  
  
  return (boxplot)
}







