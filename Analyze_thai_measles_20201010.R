### uses function Helper_Measles_30Jan2017_R to pull in some random functions that are used
library(tidyverse)

### simulation code for: Low amplitude seasonality in transmission enhances measles persistence in the vaccine era
### Authors: Amy Wesolowski (awesolowski@jhu_edu), Alex Becker

### interpolate monthly data to estimate biweekly data 
interp_monthly_data<-function(values){
  new_values<-rep(NA, (length(values)*2)+(length(values)/12)*2);
  count = 1
  add_extra<-c(seq(2,(length(values)+2),by=12), seq(7,(length(values)+2),by=12))
  for(bb in 1:(length(values)-1)){
    if(is_element(bb, add_extra)){
      new_values[c(count,count+1,count+2)] = seq(values[bb], values[bb+1], length = 4)[1:3]
      count = count+3
    }
    else{
      new_values[c(count,count+1)] = seq(values[bb], values[bb+1], length = 3)[1:2]
      count = count+2}}
  new_values[length(new_values)-1] = mean(new_values[length(new_values)-2],values[length(values)])
  new_values[length(new_values)] = values[length(values)]
  return(new_values)}

### basic measles simulation - not stochastic

forward_sim_measles<-function(Tmax, smean, imean, alpha, ps, births, beta_values, season_index){
  storeS<-rep(NA, Tmax); storeI<-rep(NA, Tmax)
  storeS[1] = smean*ps[1]
  storeI[1] = imean*ps[1]
  for(jj in 2:Tmax){
    if(jj<(5*26)){
      storeI[jj] = pmax(beta_values[season_index[jj]]*(storeI[(jj-1)]^alpha)*(storeS[(jj-1)])/ps[(jj-1)],2)
    }
    else{
      storeI[jj] = beta_values[season_index[jj]]*(storeI[(jj-1)]^alpha)*(storeS[(jj-1)])/ps[(jj-1)]
    }
    
    storeS[jj] = pmax(storeS[(jj-1)] - storeI[jj] + births[jj],2)}
  return(list(storeS = storeS, storeI = storeI))}


### basic measles simulation - stochastic
forward_sim_measles_stochastic<-function(Tmax, smean, imean, alpha, ps, births, beta_values, season_index){
  storeS<-rep(NA, Tmax); storeI<-rep(NA, Tmax)
  storeS[1] = smean*ps[1]
  storeI[1] = imean*ps[1]
  for(jj in 2:Tmax){
    if(jj<(5*26)){
      lambda_est = pmax(beta_values[season_index[jj]]*storeS[(jj-1)]*storeI[(jj-1)]^alpha/ps[(jj-1)],2)
    }
    else{
      lambda_est = beta_values[season_index[jj]]*storeS[(jj-1)]*storeI[(jj-1)]^alpha/ps[(jj-1)]
    }
    if(storeI[(jj-1)] == 0){storeI[jj] = 0}
    else{storeI[jj] = rnbinom(1, mu = lambda_est, size = storeI[(jj-1)])}
    storeS[jj] = pmax(storeS[(jj-1)] - storeI[jj] + births[jj],2)}
  return(list(storeS = storeS, storeI = storeI))}

### run simulations to create grid plots looking at different amplitudes 

rm(list=ls())
require(grid)
require(deSolve)
require(ggplot2)
require(doRNG)
require(foreach)
require(doParallel)

args<-commandArgs(TRUE)

## this tells your computer if you're in slurm or not
## this is useful so you can change WDs and inference as you need to based on your computer
inSLURM <- c(Sys_getenv('SLURM_JOB_ID') != "")

if(!inSLURM){
  args <- '20'
  setwd('~/Dropbox/GridsForAmy')
}

R0 <- as_numeric(args)

burnin <- 50
years <- 70

seir_step <- function(x, params) {
  
  delta_t <- 1/365
  time <- 1:(years/delta_t)
  S <- E <- I <- R <- N <- rep(NA,length(time))
  
  S[1] = x[1]
  E[1] = x[2]
  I[1] = x[3]
  R[1] = x[4]
  N[1] <- S[1] + E[1]+ I[1] + R[1]
  
  v = params['v']
  mu = params["mu"]
  R0 = params["R0"]
  beta1 = params["beta1"]
  sigma = params["sigma"]
  gamma = params["gamma"]
  
  beta0 <- R0 * ( sigma + mu ) / sigma * ( gamma + mu )
  
  seas <- beta0 * (1 + beta1 * cos(2 * pi * time/365))
  
  for(it in 2:length(S)){
    
    births <- v*N[it-1]
    
    dw <- runif(1,0.95,1.05)
    
    foi <- seas[it]*S[it-1]*I[it-1]/N[it-1] * dw
    inf <- sigma*E[it-1]
    rec <- gamma*I[it-1]
    
    S[it] = S[it-1] + delta_t*(births - foi - mu*S[it-1])
    E[it] = E[it-1] + delta_t*(foi - inf - mu*E[it-1])
    I[it] = I[it-1] + delta_t*(inf - rec - mu*I[it-1])
    R[it] = R[it-1] + delta_t*(rec -  mu*R[it-1])
    
    N[it] = S[it] + E[it] + I[it] + R[it]
    
  }
  
  return(data_frame(cbind(time,S,I,R,N)))
}


if(inSLURM){
  
  require(doMPI)
  
  cl <-startMPIcluster()
  
  registerDoMPI(cl)
  ## do a ton of samples if in slurm
  CBR_options <- seq(5,50,0.5)
  amp_options <- seq(0.01,0.5,0.002)
  
  
}else{
  
  require(doParallel)
  cl <- makeCluster(4,'PSOCK')
  registerDoParallel(cl)
  ## if not
  CBR_options <- seq(5,50,2)
  amp_options <- seq(0.01,0.5,0.05)
  
}

d <- expand_grid(CBR = CBR_options, amp = amp_options,R0=R0)
d$freq <- rep(NA,nrow(d))
d$PR <- NA

w1 <- getDoParWorkers()

bifur_grid <- foreach(it=1:nrow(d),
                      .combine=rbind,
                      .errorhandling="remove",
                      /options_multicore=list(set_seed=TRUE)
) %dorng% {
  
  CBR <- d[it,1]
  amp <- d[it,2]
  
  v <- CBR/1000
  
  N <- 3e6
  
  mu <- 1/50
  
  paras = c(v=v, mu = mu, N = N,
            R0 = R0, beta1 = amp,
            sigma = 365/8, gamma = 365/5)
  
  xstart = c(S = 1/paras['R0'], E = 5e-5, I = 5e-5, R = 1-1/paras['R0'])*N
  
  out <- seir_step(xstart,paras)
  
  out <- out[(burnin*365):nrow(out),]
  
  scaled_I <- sqrt(out$I/out$N)
  
  spec <- spectrum(scaled_I,plot=F)
  options <- 1/spec$freq/ 365
  power <- spec$spec
  
  
  freq <- options[which_max(power)]
  
  bi_ind_ll <-which_min(abs(options - 1.75))
  bi_ind_ul <-which_min(abs(options - 2.25))
  
  an_ind_ll <-which_min(abs(options - 0.75))
  an_ind_ul <-which_min(abs(options - 1.25))
  
  yr2 <- power[bi_ind_ll:bi_ind_ul]
  yr1 <- power[an_ind_ll:an_ind_ul]
  
  pr <- sum(yr2)/sum(yr1)
  
  data_frame(cbind(CBR,amp,R0,freq,pr,mu,N))
  
}

bifur_grid$roundfreq <- round(bifur_grid$freq)



if(inSLURM){
  save(bifur_grid,file=paste0('stoch_amp_CBR_grid_R0=',R0,'_RData'))
}


#### The next set of functions are used to simulate SIA campaigns using transmission values with different amplitudes_ 

### basic forward simulation, stochastic with an SIA
forward_sim_measles_stochastic_sia<-function(Tmax, smean, imean, alpha, ps, births, mcv_1, sia, beta_values, season_index){
  births<-(1-mcv_1)*births
  storeS<-rep(NA, Tmax)## each row is a time step
  storeI<-rep(NA, Tmax)
  storeS[1] = smean*ps[1]
  storeI[1] = imean*ps[1]
  for(jj in 2:Tmax){
    lambda_est = beta_values[season_index[jj]]*storeS[(jj-1)]*storeI[(jj-1)]^alpha/ps[(jj-1)]
    storeI[jj] = pmax(rnbinom(1, mu = lambda_est, size = storeI[(jj-1)]),1e-3)
    storeS[jj] = pmax(sia[jj]*(storeS[(jj-1)] - storeI[jj] + births[jj]),1e-3)
  }
  return(list(storeS = storeS, storeI = storeI))}

## helper function to correct negative values
correct_neg<-function(values){
  return(ifelse(values<1,1,values))}

## create beta values that are amplified
amp_betas<-function(beta_values){
  amp_beta_matrix<-matrix(,50,26) ## for 50 locations, 26 transmission values
  mean_values<-rep(mean(beta_values, na_rm=T), 26)
  dev_from_mean<-beta_values-mean_values
  amp_beta_matrix[1,] = beta_values
  amp_beta_matrix[2,] = mean_values
  x1<-rep(seq(0.001,5,length=26), 2) ## range of amplitude values
  for(ii in 2:26){
    amp_beta_matrix[ii,] = correct_neg(mean_values+dev_from_mean/x1[ii])
  }
  for(jj in 27:50){
    amp_beta_matrix[jj,] = correct_neg(mean_values+dev_from_mean*x1[jj])
  }
  return(amp_beta_matrix)
}

## run simulations without an SIA
run_sims_no_sia<-function(smean, imean, test_pops, test_births, vacc_cov_long, beta_est){
  n_runs = 25
  infected_matrix<-matrix(,n_runs,45*26)
  for(ii in 1:n_runs){
    results<-forward_sim_measles_stochastic_sia(Tmax = Tmax, smean = smean, imean = imean, alpha = 0.975, ps = test_pops, births = test_births, mcv_1 = vacc_cov_long, sia = rep(1,Tmax), beta_values = beta_est, season_index = season_index)  
    infect_values<-results$storeI[((10*26)+1):Tmax]-1
    infect_values[which(infect_values < 0)] = 0
    infected_matrix[ii,] = infect_values
  }
  return(infected_matrix)
}

### run simulations with an SIA
run_sims<-function(smean, imean, test_pops, test_births, vacc_cov_long, beta_est, sia_value){
  
  sia_none<-rep(1,Tmax)
  sia<-rep(1,Tmax); sia[seq(((10+4)*26),Tmax, by = (4*26))] = sia_value
  
  n_runs = 25
  infected_sia_matrix<-matrix(,n_runs,45*26)
  infected_no_sia_matrix<-matrix(,n_runs,45*26)
  for(ii in 1:n_runs){
    results_sia<-forward_sim_measles_stochastic_sia(Tmax = Tmax, smean = smean, imean = imean, alpha = 0.975, ps = test_pops, births = test_births, mcv_1 = vacc_cov_long, sia = sia, beta_values = beta_est, season_index = season_index)  
    infect_values_sia<-results_sia$storeI[((10*26)+1):Tmax]-1
    infect_values_sia[which(infect_values_sia<0)] = 0
    suscep_values_sia<-results_sia$storeS[((10*26)+1):Tmax]
    
    sia_none<-rep(1,Tmax)
    results_no_sia<-forward_sim_measles_stochastic_sia(Tmax = Tmax, smean = smean, imean = imean, alpha = 0_975, ps = test_pops, births = test_births, mcv_1 = vacc_cov_long, sia = sia_none, beta_values = beta_est, season_index = season_index)  
    infect_values_no_sia<-results_no_sia$storeI[((10*26)+1):Tmax]-1
    infect_values_no_sia[which(infect_values_no_sia<0)] = 0
    suscep_values_no_sia<-results_no_sia$storeS[((10*26)+1):Tmax]
    infected_sia_matrix[ii,] = infect_values_sia
    infected_no_sia_matrix[ii,] = infect_values_no_sia
  }
  return(list(sia = infected_sia_matrix, no_sia = infected_no_sia_matrix))
}

### run simulations with multiple SIAs
run_sims_vary_numb_sias<-function(smean, imean, test_pops, test_births, vacc_cov_long, beta_est, sia_value, numb_sias){
  full_numb<-seq((10+2)*26,Tmax, by=2*26)
  sia<-rep(1,Tmax); sia[full_numb[1:numb_sias]] = sia_value
  n_runs = 25
  infected_sia_matrix<-matrix(,n_runs,45*26)
  for(ii in 1:n_runs){
    results_sia<-forward_sim_measles_stochastic_sia(Tmax = Tmax, smean = smean, imean = imean, alpha = 0.975, ps = test_pops, births = test_births, mcv_1 = vacc_cov_long, sia = sia, beta_values = beta_est, season_index = season_index)  
    infect_values_sia<-results_sia$storeI[((10*26)+1):Tmax]-1
    infect_values_sia[which(infect_values_sia<0)] = 0
    suscep_values_sia<-results_sia$storeS[((10*26)+1):Tmax]
    infected_sia_matrix[ii,] = infect_values_sia
  }
  return(list(sia = infected_sia_matrix))
}

### count the number of fade outs
count_numb_fades<-function(infect_values){
  diff_infect<-diff(infect_values)
  multi_zeros<-which(diff_infect == 0 & infect_values[1:length(infect_values)-1] == 0)
  numb_multi_zeros<-length(multi_zeros)
  return(numb_multi_zeros)
}


### loop to write files -- look at odds ratio of fade outs 
sia_cov_rates<-c(0.7,0.8,0.9,0.95) ## different SIA coverage values
sia_cov_rates_neg<-1-sia_cov_rates
## start at 1 
for(jj in 1:nrow(beta_results_2)){
  print(jj)
  beta_est = beta_results_2[jj,]
  amp_sample_beta<-amp_betas(beta_est)
  cv_values<-apply(amp_sample_beta, 1, function(x) round(sd(x)/mean(x),3))
  test_pops<-c(rep(thai_pop_matrix_45[jj,1], extra_years*26), rep(thai_pop_matrix_45[jj,], each = 26))
  test_births<-c(rep(birth_values_45[1]*thai_pop_matrix_45[jj,1]/1000/26, extra_years*26), rep(birth_values_45*thai_pop_matrix_45[jj,]/1000/26, each = 26))
  smean = start_values_prop_2[jj,1]#0_04526316
  imean = start_values_prop_2[jj,2]#0_00448331
  total_value = 45*26
  loc_matrix<-matrix(,50,length(sia_cov_rates)*3)
  colnames(loc_matrix)<-c(sapply(sia_cov_rates, function(x) c(paste('lower', 100*x, sep = '-'), paste('mean', 100*x, sep = '-'), paste('upper', 100*x, sep = '-'))))
  rownames(loc_matrix)<-cv_values
  for(ii in 1:length(sia_cov_rates)){
    loc_matrix_index<-seq(ii*(3-2),ii*3)
    for(bb in 1:50){
      sim_values<-run_sims(smean, imean, test_pops, test_births, vacc_cov_long, amp_sample_beta[bb,], sia_value = sia_cov_rates_neg[ii])
      reg_amp_sia_fade<-apply(sim_values$sia, 1, function(x) count_numb_fades(x))
      reg_amp_no_sia_fade<-apply(sim_values$no_sia, 1, function(x) count_numb_fades(x))
      ratio_values<-(reg_amp_sia_fade/total_value)/(reg_amp_no_sia_fade/total_value)
      ratio_values[is_infinite(ratio_values)]<-NA
      summary_values<-t_test(ratio_values, conf_level = 0.5)
      loc_matrix[bb,loc_matrix_index] = c(summary_values$conf_int[1], summary_values$estimate, summary_values$conf_int[2])
    }
  }
}
