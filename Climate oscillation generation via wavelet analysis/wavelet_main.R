#  main script for 
#1. wavelet block knn synthetic (monthly) precipitation generator 
#2. wavelet power spectrum analysis for sub-basin across Africa

##################################################
# load package

library(dplyr)
library(pracma)
library(ggplot2)
library(reshape)
library(ggpubr)
library(RColorBrewer)
library(ggpubr)
library(RColorBrewer)
#library(tidyr)    might mask extract function is raster
library(zoo)
library(data.table)
library(gridExtra)
library(ncdf4)
library(Rcpp)
library(lubridate)
library(raster)

source("/src/data processing.R")
source('/src/block knn generator')
source('/src/block knn helper')
source('/src/block knn tune params.R')
source("/src/wavelet_africa.R")
source("/src/wavelet_basic.R")






#####################################################################################
# 1. wavelet block knn synthetic (monthly) precipitation generator


### step0: read and process precipitation timeseries data at sub-basin scale


sub_basin_his = read.csv('/data/sub_basin_his.csv')
sub_basin_fut = read.csv('/data/sub_basin_fut.csv')
basin_prec = data_process(sub_basin_his, sub_basin_fut )

### step1: generate synthetic annual time series based on oscillatory of chosen basin

file_name = "/result/"
basin_name = "basin_88"  # for paper analysis, basin 83 = sl,  77 = sm, 88 = s
basin_code = "s"           # choose from "sl", "sm", "s", , "l", "m"

#plot average power spectrum with noise spectrum add on 
if (plot_power == T) {
  prep_x = basin_prec[,basin_name]
  
  png(paste0(file_name, basin_name, "_average_power_withNoise.png"),width = 700, height = 400)
  plot_withN(prep_x, 1000, 0.8)   # 1000 simulated white noise spectrum
  dev.off()
}

# whether contain short, med, long oscillation components for the selected basin
Short = T
Med = F
Long = F
S1 = F
Resi = F

# NN =  how many simulations? NS = how many years in one synthetic timeseries?
NN = 100
NS = 100

# set params for block knn generator
# the params are tuned for three selected basins. Users could tune params for any selected basin following the code the following part
Bsize = c(10, 10, 20, 20, 20) #block size in knn

Nsize = 20  # neighbors numbers in knn 
# period range for short, medium and long oscillation component

if (basin_code == "sl" | basin_code == "l"){ 
  peaks = c(2,4,  4,9, 13,15,  30,60 )   # s1, s, m, l
} else if (basin_code == "sm" | basin_code == "m"){
  peaks =  c(2,4,  4,8, 10,16,  30,60 )  # s1, s, m, l
} else if (basin_code == "s"){
  peaks = c(0,0,  2,  7, 18, 22, 0,0)  # s1, s, m, l
}
#############################################################

# annual prec time series generator

annSyn_df_N = annual_generator(basin_name, file_name,Short, Med, Long,S1,Resi, NS, NN, peaks, Bsize, Nsize, nor = T)


### Step3: Monthly prec time series generator : monthly pattern reserved for selected basin. Our case is Mombasa, Kenya region.
# save monthly synthetic timeseries as csv file
kenya <- extent(35, 40, -5, 0)
prep_kenya = compute_yearly(pre_his,pre_fut,ken)$year_prep
wt_kenya = waveletTransN(prep_kenya, 0.8, nor = T)
month_generator(kneya, annSyn_df_N, file_name)



########### tune parms for block knn generator

# plot
# validation test for knn generator (1)simulated power spectrum vs real  (2) boxplot for <10, 10-30 and 30-60 power 

file_folder = "/plot/"
basin_name = "basin_88"
params_num = "trail_1"

# three plots
spectrum = T
Boxplots = T
Component = T
# set params to be tunes

peaks = c(0,0,  2,  7, 18, 22, 0,0) # s1, s, m, l

Bsize = c(10, 10, 20, 20, 20) #block size in knn

Nsize = 20  # neighbors numbers in knn 

test_NN = 100
test_NS = 250

# generate validation plots
tune_params (file_folder,basin_prec, basin_name,test_NS, test_NN, params_num, spectrum,Boxplots, Component,peaks, Bsize, Nsize)

############################################################################
# 2.wavelet power spectrum analysis for Sun-basin across Africa

# save relative power contributed by short, med, long period range for all sub basins across Africa: paper figure2
Africa_relative_power(basin_prec,file_name)


