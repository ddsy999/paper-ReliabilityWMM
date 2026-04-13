library(clue)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(patchwork)
library(purrr)

source("wmm_functions.R")


############################
# Algorithm 
############################

wmm_EM2 <- function(df,
                   theta,
                   maxGEMiter = 1e+3,
                   tol = 1e-6,verbose=FALSE
){
  
  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=1)
  
  N=nrow(df)
  K=ncol(gamma)
  
  time_vec  = df$time
  event_vec = df$event
  
  trace <- vector("list", length(maxGEMiter))
  
  bw = 0
  for (it in 1:maxGEMiter) {
    ### M-step ###
    new_pi = colSums(gamma)/N
    
    new_beta1 = barrier_beta1(beta[1],event_vec,time_vec,gamma,bw=bw)
    # new_beta3 = barrier_beta3(beta[3],event_vec,time_vec,gamma,bw=bw)
    new_beta2 = barrier_beta1(beta[2],event_vec,time_vec,gamma,bw=bw)
    # new_beta = c(new_beta1,new_beta2,new_beta3)
    new_beta = c(new_beta1,new_beta2)
    
    new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
    
    ### organize ###
    parameter_diff = sqrt(sum((beta-new_beta)^2))
    beta=new_beta; pi=new_pi; lambda=new_lambda;
    
    dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
    # dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
    ### Stopping rule ###
    if(parameter_diff<tol || it==maxGEMiter){
      if(verbose) cat("EM ","[",it,"]"," beta :",beta ,"\n")
      break
    }
    ### E-step ###
    gamma = weibull_estep_annealed(df,pi,lambda,beta,r=1)
    
    trace[[it]] <- list(pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1 )#, dQbeta3=dQbeta3)
  }
  
  return(list(trace = trace, lambda = lambda, beta = beta))
}


wmm_DAEM2 <- function(df,
                     theta,
                     maxGEMiter = 1e+3,
                     nsteps = 50,
                     r_init = 0.1,
                     r_end = 1,
                     tol=1e-6,verbose=FALSE
){
  method = "DAEM"
  r_grid <- exp(seq(log(r_init), log(r_end), length.out = nsteps))
  
  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=r_init)
  
  N=nrow(df)
  K=ncol(gamma)
  
  time_vec  = df$time
  event_vec = df$event
  
  trace <- vector("list", nsteps)
  
  bw = 0
  for( hyperIter in 1:nsteps){
    r =r_grid[hyperIter]
    
    for( gemIter in 1:maxGEMiter){
      ### M-step ###
      new_pi = colSums(gamma)/N
      
      new_beta1 = barrier_beta1(beta[1],event_vec,time_vec,gamma,bw=bw)
      # new_beta3 = barrier_beta3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      new_beta = c(new_beta1,new_beta2)
      # new_beta = c(new_beta1,new_beta2,new_beta3)
      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
      
      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;
      
      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){
        
        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        # dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }
    
    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1) #dQbeta3=dQbeta3)
    
  }
  
  return(list(trace = trace, pi = pi, lambda = lambda, beta = beta))
}

wmm_BM2 <- function(df,
                   theta,
                   maxGEMiter = 1e+3,
                   nsteps = 50,
                   bw_init = 1e-1,
                   bw_end = 1e-5,
                   tol=1e-6,verbose=FALSE
){
  method="BM"         
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = nsteps))
  
  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=1)
  
  N=nrow(df)
  K=ncol(gamma)
  
  time_vec  = df$time
  event_vec = df$event
  
  trace <- vector("list", nsteps)
  
  r=1
  
  for( hyperIter in 1:nsteps){
    
    bw=bw_grid[hyperIter]
    
    for( gemIter in 1:maxGEMiter){
      ### M-step ###
      new_pi = colSums(gamma)/N
      
      new_beta1 = barrier_safe_wrapper1(beta[1],event_vec,time_vec,gamma,bw=bw)
      # new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      # new_beta = c(new_beta1,new_beta2,new_beta3)
      new_beta = c(new_beta1,new_beta2)
      
      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
      
      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;
      
      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){
        
        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        # dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"\n") #,"diifbeta3:",dQbeta3,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }
    
    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1) #, dQbeta3=dQbeta3)
    
  }
  return(list(trace = trace,pi = pi, lambda = lambda, beta = beta))
}

wmm_DHEM2 = function(df,
                    theta,
                    maxGEMiter=1e+3,
                    nsteps=1e+2,
                    r_init=0.1,
                    r_end=1,
                    bw_init=1e-1,
                    bw_end=1e-5,
                    tol=1e-6,verbose=FALSE
){
  method = "DHEM"
  r_grid  <- exp(seq(log(r_init), log(r_end), length.out = nsteps))
  bw_grid <- exp(seq(log(bw_init), log(bw_end), length.out = nsteps))
  
  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=r_init)
  
  N=nrow(df)
  K=ncol(gamma)
  
  time_vec  = df$time
  event_vec = df$event
  
  trace <- vector("list", nsteps)
  
  for( hyperIter in 1:nsteps){
    r =r_grid[hyperIter]
    bw=bw_grid[hyperIter]
    
    for( gemIter in 1:maxGEMiter){
      ### M-step ###
      new_pi = colSums(gamma)/N
      
      new_beta1 = barrier_safe_wrapper1(beta[1],event_vec,time_vec,gamma,bw=bw)
      # new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      # new_beta = c(new_beta1,new_beta2,new_beta3)
      new_beta = c(new_beta1,new_beta2)
      
      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
      
      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2+(pi-new_pi)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;
      
      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){
        
        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        # dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }
    
    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1 )#, dQbeta3=dQbeta3)
    
  }
  return(list(trace = trace,pi = pi, lambda = lambda, beta = beta))
}

wmm_DHEM_adaptive2 <- function(df,
                              theta,
                              maxGEMiter=1e+3,
                              nsteps=1e+2,
                              r_init=0.1,
                              r_end=1,
                              bw_init=1e-1,
                              eta=0.1,
                              tol=1e-6,verbose=FALSE
){
  method = "adapDHEM"
  # r schedule (r -> 1)
  r_grid <- exp(seq(log(r_init), 0, length.out = nsteps))
  
  pi = theta$pi
  lambda = theta$lambda
  beta = theta$beta
  gamma = weibull_estep_annealed(df,theta$pi,theta$lambda,theta$beta,r=r_init)
  
  N=nrow(df)
  K=ncol(gamma)
  
  time_vec  = df$time
  event_vec = df$event
  
  trace <- vector("list", nsteps)
  
  bw = bw_init
  for( hyperIter in 1:nsteps){
    r <- r_grid[hyperIter]
    bw = bw*0.5
    for( gemIter in 1:maxGEMiter){
      ### M-step ###
      new_pi = colSums(gamma)/N
      
      new_beta1 = barrier_safe_wrapper1(beta[1],event_vec,time_vec,gamma,bw=bw)
      # new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      # new_beta = c(new_beta1,new_beta2,new_beta3)
      new_beta = c(new_beta1,new_beta2)
      
      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
      
      theta0 = list(pi=pi,beta=beta,lambda=lambda)
      theta1 = list(pi=new_pi,beta=new_beta,lambda=new_lambda)
      
      # Checking ACC
      acc1 = FALSE;acc2 = FALSE;
      deltaDKL = wmm_delta_DKL(df,theta0,theta1,r)
      DKL      = wmm_DKL(df,theta0,theta1)
      deltaB   = wmm_bar_diff2(theta1,theta0)
      #cat(deltaDKL-bw*wmm_bar_diff(theta1,theta0),"\n")

      if(is.na(deltaDKL-bw*deltaB)){
        # bw = 0.1*bw
        break
      }
        
      if(deltaDKL-bw*deltaB<0){
        # Acc 1st test      
        if(!is.finite(deltaDKL)||!is.finite(DKL)||DKL<0) break
        if(deltaDKL<eta*DKL) {
          #cat(deltaDKL,eta*DKL,"\n")
          break}else{acc1=TRUE}
        # Acc 2nd test
        if(bw*abs(deltaB)>eta*DKL){
          bw = min(bw,eta*DKL/abs(deltaB))
          next
        }else{acc2 = TRUE}
      }else{
        acc1=TRUE;acc2=TRUE;
      }
      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2))
      #print(parameter_diff)
      beta=new_beta; pi=new_pi; lambda=new_lambda;
      dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
      # dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
      
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1," para diff: ", parameter_diff,"\n")
        break
      }
    }
    
    if(acc1&&acc2){
      trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1 )#, dQbeta3=dQbeta3)
      last_acc = list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1)#, dQbeta3=dQbeta3)
    }
  }
  
  return(list(trace = trace,pi = last_acc$pi, lambda = last_acc$lambda, beta = last_acc$beta,bw=last_acc$bw,dQbeta1 = last_acc$dQbeta1, dQbeta3=last_acc$dQbeta3))
}



#### ---------------------------------------------------------------------------
# Data Load
#### ---------------------------------------------------------------------------

# df <- read.table("DATA\\Aarest_data.txt", header = TRUE)
# df <- read.table("DATA\\FRT_censord.txt", header = TRUE)
df <- read.table("DATA\\LFP.txt", header = TRUE)
# df <- read.table("DATA\\SerumReversal.txt", header = TRUE)

#### ---------------------------------------------------------------------------
# Hyper-parameter (set here first)
#### ---------------------------------------------------------------------------
maxGEMiter = 1e+6
nsteps     <- 100
r_init     <- 0.1
r_end      <- 1
bw_init    <- 1e-1
bw_end     <- 1e-8   # 필요시 조정
eta        <- 0.1    # adaptive DHEM 전용
errtol     = 1e-10

#### ---------------------------------------------------------------------------
# Init-parameter
#### ---------------------------------------------------------------------------
K <- 2
pi_init   <- rep(1 / K, K)
beta_init <- c(0.2, 1)   
lambda_init <- c(1,1e-5) #wmm_lambda_init(df$time,df$event, beta_init,ratio1=0.3,ratio3=0.7)
theta_init = list(beta=beta_init,pi=pi_init,lambda=lambda_init)

#### ---------------------------------------------------------------------------
# Train (one run each)
#### ---------------------------------------------------------------------------
verbose=TRUE
# 1) EM (standard EM: r=1, bw=0)  -- wmm_EM은 내부에서 wmm_em_at_r_bw 사용
fit_EM <- wmm_EM2(df,theta=theta_init,maxGEMiter = maxGEMiter,tol=errtol,verbose=verbose)
print("end EM")
# 2) DAEM (bw=0, r schedule)
# fit_DAEM <- wmm_DAEM2(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps = nsteps,r_init = r_init,r_end=r_end,tol=errtol,verbose=verbose)
# print("end DAEM")
# 3) Barrier method (r=1, bw schedule)
# fit_BM <- wmm_BM2(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps=nsteps,bw_init = bw_init,bw_end = bw_end,tol = errtol,verbose=verbose)
# print("end BM")
# 4) DHEM (r schedule + bw schedule)
fit_DHEM <- wmm_DHEM2(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps=nsteps,r_init = r_init,r_end=r_end,bw_init = bw_init,bw_end = bw_end,tol = errtol,verbose=verbose)
print("end DHEM")
# 5) Adaptive DHEM (r schedule + adaptive bw control)
fit_adapDHEM <- wmm_DHEM_adaptive2(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps=nsteps,r_init = r_init,r_end=r_end,bw_init = bw_init,eta=eta,tol = errtol,verbose=verbose)
print("end adapDHEM")
# fit_adapBM = wmm_BM_adaptive(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps=nsteps,bw_init = bw_init,bw_end = bw_end,tol = errtol,verbose=verbose)
# print("end adapBM")



#### ---------------------------------------------------------------------------
# Quick check (final params)
#### ---------------------------------------------------------------------------

# list(pi = fit_DAEM$pi, lambda = fit_DAEM$lambda, beta = fit_DAEM$beta)
# list(pi = fit_BM$pi,   lambda = fit_BM$lambda,   beta = fit_BM$beta)
# list(pi = fit_DHEM$pi, lambda = fit_DHEM$lambda, beta = fit_DHEM$beta)
# list(pi = fit_adapDHEM$pi, lambda = fit_adapDHEM$lambda, beta = fit_adapDHEM$beta, bw = fit_adapDHEM$bw)
# list(pi = fit_adapBM$pi, lambda = fit_adapBM$lambda, beta = fit_adapBM$beta, bw = fit_adapBM$bw)

#### ---------------------------------------------------------------------------
# ggplot
#### ---------------------------------------------------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# trace(list) -> data.frame
# Columns: nstep, r, bw, pi1..piK, lambda1..lambdaK, beta1, beta3, dQbeta1, dQbeta3
# - pi / lambda의 list-column은 만들지 않고, 전부 숫자 컬럼으로 펼칩니다.
trace_to_df2 <- function(fit, nsteps, r_init, K = 2) {
  tr <- fit$trace
  if (is.null(tr) || length(tr) == 0) stop("fit$trace가 비어있습니다.")
  
  # r_grid (nsteps 기준으로 고정)
  r_grid <- exp(seq(log(r_init), 0, length.out = nsteps))
  
  rows <- vector("list", length(tr))
  idx <- 0L
  
  for (t in seq_along(tr)) {
    x <- tr[[t]]
    if (is.null(x)) next  # NULL step은 행을 만들지 않음
    
    idx <- idx + 1L
    
    out <- data.frame(
      nstep   = t,
      r       = if (t >= 1 && t <= nsteps) r_grid[t] else NA_real_,
      bw      = if (!is.null(x$bw)) x$bw else NA_real_,
      beta1   = if (!is.null(x$beta) && length(x$beta) >= 1) x$beta[1] else NA_real_,
      dQbeta1 = if (!is.null(x$dQbeta1)) x$dQbeta1 else NA_real_
    )
    
    # pi (pi1..piK)
    for (k in 1:K) out[[paste0("pi", k)]] <- NA_real_
    if (!is.null(x$pi)) {
      pi <- x$pi
      if (length(pi) < K) pi <- c(pi, rep(NA_real_, K - length(pi)))
      for (k in 1:K) out[[paste0("pi", k)]] <- pi[k]
    }
    
    # lambda (lambda1..lambdaK)
    for (k in 1:K) out[[paste0("lambda", k)]] <- NA_real_
    if (!is.null(x$lambda)) {
      lam <- x$lambda
      if (length(lam) < K) lam <- c(lam, rep(NA_real_, K - length(lam)))
      for (k in 1:K) out[[paste0("lambda", k)]] <- lam[k]
    }
    
    rows[[idx]] <- out
  }
  
  if (idx == 0L) stop("업데이트(= non-NULL trace)가 없습니다.")
  
  df <- do.call(rbind, rows[seq_len(idx)])
  df <- df[order(df$nstep), , drop = FALSE]
  rownames(df) <- NULL
  df
}


trace_to_df_EM2 <- function(trace) {
  n <- length(trace)
  
  out <- data.frame(
    nstep   = seq_len(n),
    beta1   = NA_real_,
    beta2   = NA_real_,
    dQbeta1 = NA_real_,
    pi1     = NA_real_,
    pi2     = NA_real_,
    lambda1 = NA_real_,
    lambda2 = NA_real_
  )
  
  for (i in seq_len(n)) {
    ti <- trace[[i]]
    
    out$beta1[i]   <- ti$beta[1]
    out$beta2[i]   <- ti$beta[2]
    out$dQbeta1[i] <- ti$dQbeta1
    
    out$pi1[i] <- ti$pi[1]
    out$pi2[i] <- ti$pi[2]
    
    out$lambda1[i] <- ti$lambda[1]
    out$lambda2[i] <- ti$lambda[2]
  }
  
  out
}



df_EM <- trace_to_df_EM2(fit_EM$trace)
# df_DAEM <- trace_to_df2(fit_DAEM,nsteps = nsteps,r_init=r_init, K = 2)
# df_BM <- trace_to_df2(fit_BM,nsteps = nsteps,r_init=r_init, K = 2)
df_DHEM <- trace_to_df2(fit_DHEM,nsteps = nsteps,r_init=r_init, K = 2)
df_adapDHEM <- trace_to_df2(fit_adapDHEM,nsteps = nsteps,r_init=r_init, K = 2)


find_stationary_const <- function(df, method = "L1") {
  
  # norm 계산
  if (method == "L1") {
    score <- abs(df$dQbeta1)
  } else if (method == "L2") {
    score <- sqrt(df$dQbeta1^2 )
  } else {
    stop("method must be 'L1' or 'L2'")
  }
  
  # 최소값 행 반환
  stationary_row <- df[which.min(score), ]
  
  return(stationary_row)
}
find_stationary_adap_const<- function(df, tol = 1e-4) {
  
  score <- abs(df$dQbeta1) 
  
  candidates <- df[score <= tol, ]
  
  if (nrow(candidates) == 0) {
    return(NULL)
  }
  
  tail(candidates, 1)
}


make_stationary_table_const <- function(df_DHEM, df_adapDHEM) {
  
  dhe <- find_stationary_const(df_DHEM)
  adap <- find_stationary_adap_const(df_adapDHEM)
  
  dhe$Algorithm  <- "DHEM"
  adap$Algorithm <- "Adaptive DHEM"
  
  # 열 순서 정리 (Algorithm 맨 앞으로)
  dhe  <- dhe[, c("Algorithm", names(dhe)[names(dhe) != "Algorithm"])]
  adap <- adap[, c("Algorithm", names(adap)[names(adap) != "Algorithm"])]
  
  rbind(dhe, adap)
}

stationary_table <- make_stationary_table_const(df_DHEM, df_adapDHEM)
stationary_table
df_EM %>% tail(1)






plot_dQbeta_trace <- function(df, title = NULL,
                              axis_text_y_size  = 14,
                              axis_title_y_size = 20,
                              axis_text_x_size  = 14,
                              axis_title_x_size = 14,
                              title_size        = 20
) {
  
  base_theme <- theme_bw() +
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size,angle =0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )
  
  p1 <- ggplot(df, aes(x = r, y = dQbeta1)) +
    geom_line() +
    coord_cartesian(xlim = c(r_init, 1))+
    labs(x = "Annealing parameter", y = expression(nabla~Q~beta[1])) +
    base_theme
  
  # p3 <- ggplot(df, aes(x = r, y = dQbeta3)) +
  #   geom_line() +
  #   coord_cartesian(xlim = c(r_init, 1))+
  #   labs(x = "Annealing parameter", y = expression(nabla~Q~beta[3])) +
  #   base_theme
  (p1 ) + plot_annotation(title = title) &
    theme(plot.title = element_text(size = title_size))
  # (p1 / p3) + plot_annotation(title = title) &
  #   theme(plot.title = element_text(size = title_size))
}

# 사용 예
# plot_dQbeta_trace(df_DAEM, title = expression(DAEM * ": " * nabla~beta * " trace"))
# plot_dQbeta_trace(df_DHEM, title = expression(DHEM * ": " * nabla~beta* " trace"))
# plot_dQbeta_trace(df_adapDHEM, title = expression(adap*" "*DAEM * ": " * nabla~beta * " trace"))
# plot_dQbeta_trace(df_BM, title = expression(Barrier*" "*method * ": " * beta * " trace"))

(plot_beta_trace(df_DAEM, title = expression(DAEM * ": " * beta * " trace"))|plot_dQbeta_trace(df_DAEM, title = expression(DAEM * ": " * nabla~beta * " trace")))+
  plot_annotation(title= "DAEM") &  theme(plot.title = element_text(size = 20, hjust = 0.5))

(plot_beta_trace(df_BM, title = expression(Barrier*" "*method * ": " * beta * " trace"))|plot_dQbeta_trace(df_BM, title = expression(DHEM * ": " * nabla~beta * " trace")))+
  plot_annotation(title= "Barrier method")&  theme(plot.title = element_text(size = 20, hjust = 0.5))

(plot_beta_trace(df_DHEM, title = expression(DHEM * ": " * beta * " trace"))|plot_dQbeta_trace(df_DHEM, title = expression(DHEM * ": " * nabla~beta * " trace")))+
  plot_annotation(title= "DHEM")&  theme(plot.title = element_text(size = 20, hjust = 0.5))

(plot_beta_trace(df_adapDHEM, title = expression(adapDHEM * ": " * beta * " trace"))|plot_dQbeta_trace(df_adapDHEM, title = expression(adapDHEM* ": " * nabla~beta * " trace")))+
  plot_annotation(title= "Adaptive DHEM")&  theme(plot.title = element_text(size = 20, hjust = 0.5))


