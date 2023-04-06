# process sub_basin file
# read his and fut file
data_process = function(sub_basin_his,sub_basin_fut ){

  # add year column for future periods
  sub_basin_fut = sub_basin_fut[1:1032,]
  Year = rep(2015:2100, each = 12)
  sub_basin_fut = cbind(Year, sub_basin_fut)
  # combine his and future
  sub_basin  = rbind(sub_basin_his, sub_basin_fut)
  
  # aggregate to annual prep
  basin_prep = aggregate(sub_basin[,3:101], by = list(sub_basin$Year), FUN = sum)
  
  colnames(basin_prep)[1] = 'Year'
  #unit conversion
  basin_prep[, 2:100] = basin_prep[, 2:100]*86400*365/12
  
  return (basin_prep)
}

####################################################
# process CMIP6 climate model grid scale output
# extract and process prec timeseries for selected region
#CMIP6 climate data process and basin data extract

## region_annual

region_monthly = function(basin){
  raw_prec = read_CMIP6(read_cmip = T)
  pre_his = raw_prec[[1]]
  pre_fut = raw_prec[[2]]
  
  prec_annual = compute_yearly(pre_his,pre_fut,basin)
  
  return(prec_annual)
  
}

## read CMIP6 data

read_CMIP6 = function(read_cmip = T){
  library(ncdf4)
  library(Rcpp)
  library(lubridate)
  library(dplyr)
  library(raster)
  # CMIP6 basin data extract
  file_loc = '/data/'
  
  file_name = c('pr_Africa_FIO-ESM-2-0_historical_r1i1p1f1_gn_18500116-20141216_v20191209.nc','pr_Africa_FIO-ESM-2-0_ssp245_r1i1p1f1_gn_20150116-21001216_v20191226.nc')
  
  paste(file_loc, file_name[1], sep="")  # historical Africa CMIP6 prec
  paste(file_loc, file_name[2], sep="")   # future Africa CMIP6 Prec
  
  # use raster to read data
  pre_his <- brick(paste(file_loc,file_name[1],sep=''), varname="pr")
  
  pre_fut <- brick(paste(file_loc,file_name[2],sep=''), varname="pr")
  
  detach(ncdf4)
  detach(Rcpp)
  detach(lubridate)
  detach(dplyr)
  library(raster)
  
  return(list(pre_his, pre_fut))
}



# function: compute_yearly()
# input: pre_his, pre_fut,basin
# compute annual prep  and do the unit conversion: *86400*365/12
# output: yearly_prep(three columns: years,year_prep_raw,year_prep)

compute_yearly = function(pre_his,pre_fut,basin){
  
  pre_df = basin_pre(pre_his,pre_fut,basin)
  # add years and months columns
  years <- rep(1850:2100,each =12)
  months <- rep(c(1:12),times=251)
  pre_df$years = years
  pre_df$months = months
  
  # aggregate to annual data
  yearly_prep <- pre_df %>%
    group_by(years) %>%
    summarise(year_prep_raw = sum(pre_raw))
  
  # unit conversion
  yearly_prep$year_prep = yearly_prep$year_prep_raw*86400*365/12
  
  
  return(yearly_prep)
}


# # function: basin_pre()
# input: pre_his, pre_fut,basin(an extent object)
# do the extract, spatial average and combine hist and future
# output: pre_df . dataframe. one column indicate pre_raw for each month, all years
basin_pre = function(pre_his,pre_fut,basin){
  # extract
  
  pre_b_his = data.frame(extract(pre_his,basin))
  pre_b_fut = data.frame(extract(pre_fut,basin))
  
  # spatial average
  ave_his <- data.frame(colMeans(pre_b_his))
  colnames(ave_his) <- 'pre_raw'
  
  ave_fut <- data.frame(colMeans(pre_b_fut))
  colnames(ave_fut) <- 'pre_raw'
  
  # bind historical and future dataframe together
  pre_df = rbind(ave_his,ave_fut)
  
  return (pre_df)
  
}
