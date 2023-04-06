# find_knn(f_vector, ts , K)
# K: K neighbors (default:K = as.integer(sqrt(N)) )
# input: feature vector of size B, a time series ts
# two outout: selected neighbor (vector of size B' associated last time index TT)

# then the subsequent values starting from TT+1.... will be the simulated values and next feature vector

find_knn = function (f_vector, ts, K){
  
  B = length(f_vector) 
  N = length(ts)
  
  vec_m = all_vector(B,ts)
  vec_df = as.data.frame(vec_m)
  
  # step 2: compute distance vector  
  # step 1 : build all length B vectors of time
  
  distance = apply(vec_df,1,ec_dis)
  
  #K = as.integer(sqrt(N)) # or should be sqrt(total blocks = N-B+1) # try increase number of k
  
  ordered_index = order(distance)[2:(K+1)]  # which block has the smallest, 2nd smallest... ?
  # omit 1st smallest, because that would be the block itself
  
  # step 4:  weighed sample a vector from ordered_index
  
  sample_index = weighted_sample(K,ordered_index)
  
  # step 5: the selected neighbor and its last index T
  
  TT = sample_index + B -1
  
  # add a guard: TT must < (N-B) 
  while (TT > (N-B)){
    sample_index = weighted_sample(K,ordered_index)
    TT = sample_index + B -1
  }
  
  sample_neighbor = ts[sample_index:T]
  
  return (TT)
  
}


# helper function: build all vectors of size B of a time serie
# all_vector(B,ts)
# output: a matrix. each row is a length B vector

all_vector = function (B,ts){
  
  N = length(ts)
  num_block = N - B + 1 # number of length B vectors (blocks)
  
  vec_m = matrix (0,nrow = num_block, ncol = B) # init a matrix
  
  for (i in 1:num_block){
    
    vec_m[i,] = ts[i:(i+B-1)]
    
  }
  
  return (vec_m)
  
  
}

# helper function: compute euclidean distance of a vector with feature vector
# ec_dis(vec, feature vector)
# output: is a value : distance

ec_dis = function(vec){
  
  distance_v = sqrt(sum((vec-f_vector)^2))
  
  return (distance_v)
}

# helper function: weighted sample a vector from ordered_index
# weighted_sample(k,ordered_index)
# output: sample index ( a value). the index is xxth vector in vec_df

weighted_sample = function(k,ordered_index){
  
  # build kernel
  w = 1:k
  w = 1/w
  w = w/sum(w)  # the probability vector of each index being selected  
  
  # try this:  w = 1/k   make each neighbor same prob to be selected
  
  #w = rep(1/k, times = k)
  w = cumsum(w) # build cdf of the above probability
  
  # generate uniform random number in [0,1]
  uni = runif(1,0,1)
  
  # map that random number to which index
  
  uni_cdf = c(uni,w)
  map_index = rank(uni_cdf)[1]  # the rank of the uni random number is the corresponding index
  
  sample_index = ordered_index[map_index]
  
  return (sample_index)
  
}


# function to implement knn to build synthetic timeseries
# knn_simulation(ts,NS, lower.t, upper.t)
# ts: raw timeseries 
# NS: length of desired synthetic series
# lower.t, upper.t: lower and upper period of input component series
# K:number of neighbors

# output: synthetic time series vector

knn_simulation = function(ts, NS, lower.t, upper.t, K, Bsize){
  
  
  
  # step 0. init
  N = length(ts)
  #B = as.integer((lower.t + upper.t) / (2*4))
  B = Bsize
  
  
  BS = as.integer(NS/B)  # how many blocks to simulate
  
  # step 1: random select a start Ts
  Ts = sample(B:(N-B), 1)
  # the subsequent values become the first block in synthetic series and the feature vector
  f_vector <<- ts[(Ts+1):(Ts+B)]    # set f_vector as global variable (because need to be used in other function )
  ts_syn = f_vector
  
  # step 2: loop to add more blocks into ts_syn
  for ( i in 1:BS){
    
    Ts = find_knn(f_vector,ts, K)
    f_vector <<- ts[(Ts+1) : (Ts+B)]
    
    ts_syn = append(ts_syn, f_vector)
    
    i = i +1 
    
  }
  
  # step 3: let the length of ts_syn  = NS
  
  ts_syn = ts_syn[1:NS]
  
  return (ts_syn)
  
  
}


# function  simulate multiple synthetic time series

# knn_generator (ts,NS, lower.t, upper.t, NN, K)
# NN = number of synthetic time series
# NS = length of synthetic time series

# output: dataframe, each column is a time series (NS * NN)

# note the lower.t and upper.t = 4 is the min value could be select .B =1

knn_generator = function(ts, NS, lower.t, upper.t, NN, K, Bsize){
  
  syn_df = matrix(0,nrow = NS, ncol = NN)
  
  set.seed(6)
  for ( i in 1:NN){
    #set.seed(i)
    
    syn_df[,i] = knn_simulation(ts, NS,lower.t,upper.t, K, Bsize)
    
    
  }
  
  syn_df = as.data.frame(syn_df)
  return (syn_df)
  
}

# function simulate annual time series

# annual_generator(wt, NS, NN)
# NS = length of synthetic time series
# NN = number of synthetic time series
# peaks (vector) = lower and upper of peaks
# Bsize: block size for block knn
# K: number of nerighbors


# output: annSyn_df. dataframe. each column is a synthetic annual time series. NS * NN
# output: individual short, mid, long ossi dataframe

# synthetic annual series are at normalized space


# could also return a list of annSyn_df and annOsc_df (only add up component simulations)
# customized periods you want here



annual_generator = function(wt, NS, NN, peaks, Bsize, K, nor){
  
  # peaks variables
  s1.low = peaks[1]
  s1.upp = peaks[2]
  s.low = peaks[3]
  s.upp = peaks[4]
  m.low = peaks[5]
  m.upp = peaks[6]
  l.low = peaks[7]
  l.upp = peaks[8]
  
  
  # step 1 : get raw_ts of each period component [could customized] (omit edge data)
  raw_ts_s1 = ts_builder(s1.low, s1.upp, wt, nor = nor )
  
  raw_ts_s = ts_builder(s.low, s.upp,wt, nor = nor)
  
  raw_ts_m = ts_builder(m.low, m.upp, wt, nor = nor)
  
  raw_ts_l = ts_builder(l.low, l.upp, wt, nor = nor) 
  
  ## residuals series
  # whether consider normalized or non-normalized series
  ts_detrend_nor = ts_builder(1, 1, wt, noise = TRUE,nor = nor)  # full length time series 
  ts_s1 = ts_sub(wt, s1.low, s1.upp, rescale = FALSE)$series$year_prep.r 
  ts_s = ts_sub(wt, s.low, s.upp, rescale = FALSE)$series$year_prep.r   
  ts_m = ts_sub(wt, m.low, m.upp, rescale = FALSE)$series$year_prep.r
  ts_l = ts_sub(wt, l.low, l.upp, rescale = FALSE)$series$year_prep.r
  
  raw_resi = ts_detrend_nor - ts_s1 - ts_s - ts_m - ts_l
  
  # step 2: call knn_generator to build simulation matrix for each component
  
  
  # whether slice raw_ts and want to make adjustment
  
  ts_df_s1 = knn_generator(raw_ts_s1, NS, lower.t, upper.t, NN, K, Bsize[1] ) 
  ts_df_s = knn_generator(raw_ts_s[20:200], NS, lower.t, upper.t, NN, K, Bsize[2])
  ts_df_m = knn_generator(raw_ts_m[100:200], NS, lower.t, upper.t, NN, K, Bsize[3])
  ts_df_l = knn_generator( raw_ts_l[0:150], NS, lower.t, upper.t, NN, K, Bsize[4])   
  ts_df0 = knn_generator(raw_resi, NS, lower.t, upper.t, NN, K, Bsize[5])
  
  # smooth medium osi
  sm_ts_df_m = smooth_df(ts_df_m, 0.3, NS, NN)
  # smooth long osi
  sm_ts_df_l = smooth_df(ts_df_l, 0.5, NS, NN) 
  
  # step 3: add these component dataframe together
  
  annOsc_df = ts_df_s  #+ sm_ts_df_l 
  
  annSyn_df = ts_df_s1 + ts_df_s  + sm_ts_df_m + sm_ts_df_l  +ts_df0
  #annSyn_df = sm_ts_df_m + ts_df0
  #annSyn_df = sm_ts_df_l + ts_df0
  
  return (list(annSyn_df, ts_df_s1, ts_df_s , sm_ts_df_m, sm_ts_df_l, ts_df0))
  
  
}


# helper function :ts_sub(wt,lower,upper,rescale)
# get sub component series

ts_sub = function(wt,lower,upper,rescale){
  
  
  
  sub_series = reconstruct(wt,only.sig = FALSE,sel.lower = lower,sel.upper = upper,plot.rec = FALSE,rescale = rescale)
  
  
  return (sub_series)
}


# helper function : ts_builder(wt, lower.t, upper.t, wt,noise(T or F))

# return the raw_ts time series of specfic period ready for use for annual simulation
## if noise == TRUE, return normalized detrend time series, which will be used for generate noise series

# raw_ts: omit the block at both edge of time period


ts_builder = function(lower.t, upper.t, wt, noise = FALSE, nor){
  
  
  
  if (noise == FALSE){
    
    sub_ts = ts_sub(wt, lower.t, upper.t, rescale = FALSE)
    omit = as.integer((lower.t + upper.t) / 2)
    N = length(sub_ts$series$year_prep.r)
    raw_ts = sub_ts$series$year_prep.r[(omit + 1): (N - omit)]
    
    return (raw_ts)
    
  } else{
    
    # return original normalized detrend series
    
    # return original detrend series (not normalized)
    ts_detrend = wt$series$year_prep
    
    ts_detrend_nor = (ts_detrend- mean(ts_detrend)) / sd(ts_detrend)
    
    if (nor == T){
      return (ts_detrend_nor)
      
    } else{
      
      return(ts_detrend)
    }
    
    
  }
  
  
}


# helper function smooth_df(osi_df, spar)
# given a simualted osi df, return a smoothed version df
# sp: smooth parameter. higher, more smoother

smooth_df = function(osi_df, sp, NS, NN){
  
  smo = matrix(0, NS, NN)
  
  for (i in 1 : NN){
    smoothingSpline = smooth.spline(1:NS, osi_df[,i], spar= sp)
    smo[,i] = smoothingSpline$y
  }
  
  smo_df = data.frame(smo)
  return(smo_df)
}


#### disaggregate to monthly

# function  monthly_knn(wt, pre_his, pre_fut, basin,annPsyn, Trend (T or F))
# input: annPsyn: a vector, len = NS. simulated annual time series (normalized space)
# basin: ken 
#wt: ken's wt

# output:monthPsyn:  dataframe (NS * 12). each value is monthly simulated values (detrend)

# if annOsc, then output monthly osscillation


monthly_knn = function(wt, pre_his, pre_fut, basin, annPsyn, Trend = FALSE, annOsc = 0){
  
  # step 1. monthly proportion matrix of Kenya
  P = month_p(pre_his, pre_fut, basin)
  
  # step 2. original (Kenya)detrend annual P
  
  #wt = waveletco(pre_his,pre_fut,basin,loess=0.8,only_his=FALSE) # compute outside to save time
  
  annP = wt$series$year_prep  
  
  # step 3. rescale annPsyn to reflect original ken's mean and variability
  
  annPsyn = (annPsyn - mean(annPsyn))/ sd(annPsyn) * sd(annP) + mean(annP) 
  
  
  # if want add trend back 
  if (Trend == TRUE){
    
    
    # add first 100 year mean
    annPMean = mean(wt$series$year_prep.trend[1:100])
    annPsyn = annPsyn + annPMean
    annP = annP + annPMean
  }
  
  
  # step 4. resample using KNN
  
  N = length(annP)
  K = as.integer(sqrt(N))
  NS = length(annPsyn)
  # init
  monthPsyn = data.frame(matrix(0,nrow = NS, ncol =12 ))
  
  for (i in 1: NS){
    
    # find k nearest neighbors in historical annP: according to distance between annPsyn[i] and annP
    knn_index = order(abs(annPsyn[i] - annP))[1:K]
    
    # weighted resample an index from knn_index
    sample_index =  weighted_sample(K, knn_index)
    
    # fill monthPsyn[i,] based on monthly proportion matrix at year = sample_index
    
    # whether simulate timeseries or oscillation
    if (length(annOsc) == 1){
      monthPsyn[i,] = annPsyn[i] * P[sample_index,]   # simulate time series
      
    } else{
      
      monthPsyn[i,] = annOsc[i] * P[sample_index,]    # simulate oscillation
    }
    
    
  }
  
  
  return (monthPsyn)
}






# helper function generate monthly proportion matrix(N*12)
# month_p(pre_his,pre_fut,basin)

# output: proportion matrix of monthly prep P
# P: dataframe. N* 12. 

month_p = function(pre_his, pre_fut, basin){
  
  pre_df = basin_pre(pre_his,pre_fut,basin)
  
  # unit conversion from mm/s to mm/month
  pre_df$pre = pre_df$pre_raw * 86400*365/12
  # reshape to matrix (N*12)
  
  monthPMatrix = matrix(pre_df$pre, ncol = 12, byrow = TRUE)
  
  # calculate rowsum (annual value)   # a vector
  annSum = rowSums(monthPMatrix)
  
  # compute proportion matrix of monthQ
  P = data.frame(monthPMatrix / annSum)
  
  return (P)
  
}


