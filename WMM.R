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

wmm_EM <- function(df,
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
    new_beta3 = barrier_beta3(beta[3],event_vec,time_vec,gamma,bw=bw)
    new_beta2 = barrier_beta2(beta[2],event_vec,time_vec,gamma,bw=bw)
    new_beta = c(new_beta1,new_beta2,new_beta3)
    
    new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
  
    ### organize ###
    parameter_diff = sqrt(sum((beta-new_beta)^2))
    beta=new_beta; pi=new_pi; lambda=new_lambda;

    dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
    dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
    
    trace[[it]] <- list(pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)
    
    ### Stopping rule ###
    if(parameter_diff<tol || it==maxGEMiter){
      if(verbose) cat("EM ","[",it,"]"," beta :",beta ,"\n")
      break
    }
    ### E-step ###
    gamma = weibull_estep_annealed(df,pi,lambda,beta,r=1)

    
  }

  return(list(trace = trace, lambda = lambda, beta = beta))
}


wmm_DAEM <- function(df,
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
      new_beta3 = barrier_beta3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      new_beta = c(new_beta1,new_beta2,new_beta3)
      
      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
    
      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;

      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){
        
        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }

    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)

  }

  return(list(trace = trace, pi = pi, lambda = lambda, beta = beta))
}

wmm_BM <- function(df,
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
      new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      new_beta = c(new_beta1,new_beta2,new_beta3)
      
      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
    
      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;

      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){
        
        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }

    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)

  }
  return(list(trace = trace,pi = pi, lambda = lambda, beta = beta))
}

wmm_DHEM = function(df,
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
      new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
      new_beta2 = 1
      new_beta = c(new_beta1,new_beta2,new_beta3)
      
      new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
    
      ### organize ###
      parameter_diff = sqrt(sum((beta-new_beta)^2))
      beta=new_beta; pi=new_pi; lambda=new_lambda;

      ### Stopping rule ###
      if(parameter_diff<tol || gemIter==maxGEMiter){
        
        dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
        dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
        if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3,"\n")
        break
      }
      ### E-step ###
      gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    }

    trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)

  }
  return(list(trace = trace,pi = pi, lambda = lambda, beta = beta))
}

wmm_DHEM_adaptive <- function(df,
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
    new_beta3 = barrier_safe_wrapper3(beta[3],event_vec,time_vec,gamma,bw=bw)
    new_beta2 = 1
    new_beta = c(new_beta1,new_beta2,new_beta3)
    
    new_lambda = sapply(1:K , function(i)  sum(gamma[,i]*event_vec)/sum(gamma[,i]*(time_vec^new_beta[i])))
  
    theta0 = list(pi=pi,beta=beta,lambda=lambda)
    theta1 = list(pi=new_pi,beta=new_beta,lambda=new_lambda)

    # Checking ACC
    acc1 = FALSE;acc2 = FALSE;
    deltaDKL = wmm_delta_DKL(df,theta0,theta1,r)
    DKL      = wmm_DKL(df,theta0,theta1)
    deltaB   = wmm_bar_diff(theta1,theta0)
    
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
    # parameter_diff = abs(abs((beta-new_beta)))
    #print(parameter_diff)
    beta=new_beta; pi=new_pi; lambda=new_lambda;
    dQbeta1 = diffB_onlyB(beta[1],event_vec,time_vec,gamma,1)
    dQbeta3 = diffB_onlyB(beta[3],event_vec,time_vec,gamma,3)
     
    ### E-step ###
    gamma = weibull_estep_annealed(df,pi,lambda,beta,r=r)
    ### Stopping rule ###
    if(parameter_diff<tol || gemIter==maxGEMiter){
      if(verbose) cat(method,"[Hpyer iter: ",hyperIter,"]","[GEM iter: ",gemIter,"]"," beta :",beta ," r:",r , " bw:",bw,"diifbeta1:",dQbeta1,"diifbeta3:",dQbeta3," para diff: ", parameter_diff,"\n")
      break
    }
    }

    if(acc1&&acc2){
      trace[[hyperIter]] <- list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)
      last_acc = list(bw=bw,r = r, pi = pi, lambda = lambda, beta = beta, dQbeta1 = dQbeta1, dQbeta3=dQbeta3)
    }
  }

  return(list(trace = trace,pi = last_acc$pi, lambda = last_acc$lambda, beta = last_acc$beta,bw=last_acc$bw,dQbeta1 = last_acc$dQbeta1, dQbeta3=last_acc$dQbeta3))
}


#### ---------------------------------------------------------------------------
# Data Load
#### ---------------------------------------------------------------------------

# df <- read.table("DATA\\Aarest_data.txt", header = TRUE)
df <- read.table("DATA\\FRT_censord.txt", header = TRUE)
# df <- read.table("DATA\\LFP.txt", header = TRUE)
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
K <- 3
pi_init   <- rep(1 / K, K)
beta_init <- c(0.2, 1,5)   
lambda_init <- wmm_lambda_init(df$time,df$event, beta_init,ratio1=0.3,ratio3=0.7)
theta_init = list(beta=beta_init,pi=pi_init,lambda=lambda_init)

#### ---------------------------------------------------------------------------
# Train (one run each)
#### ---------------------------------------------------------------------------
verbose=TRUE
# 1) EM (standard EM: r=1, bw=0)  -- wmm_EM은 내부에서 wmm_em_at_r_bw 사용
fit_EM <- wmm_EM(df,theta=theta_init,maxGEMiter = maxGEMiter,tol=errtol,verbose=verbose)
print("end EM")
# 2) DAEM (bw=0, r schedule)
# fit_DAEM <- wmm_DAEM(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps = nsteps,r_init = r_init,r_end=r_end,tol=errtol,verbose=verbose)
# print("end DAEM")
# 3) Barrier method (r=1, bw schedule)
# fit_BM <- wmm_BM(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps=nsteps,bw_init = bw_init,bw_end = bw_end,tol = errtol,verbose=verbose)
# print("end BM")
# 4) DHEM (r schedule + bw schedule)
fit_DHEM <- wmm_DHEM(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps=nsteps,r_init = r_init,r_end=r_end,bw_init = bw_init,bw_end = bw_end,tol = errtol,verbose=verbose)
print("end DHEM")
# 5) Adaptive DHEM (r schedule + adaptive bw control)
fit_adapDHEM <- wmm_DHEM_adaptive(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps=nsteps,r_init = r_init,r_end=r_end,bw_init = bw_init,eta=eta,tol = errtol,verbose=verbose)
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

trace_to_df<- function(fit, nsteps, r_init, K = 3) {
  tr <- fit$trace
  if (is.null(tr) || length(tr) == 0) stop("fit$trace가 비어있습니다.")

  # r_grid (nsteps 기준으로 고정)
  r_grid <- exp(seq(log(r_init), 0, length.out = nsteps))

  rows <- vector("list", length(tr))
  idx <- 0L

  for (t in seq_along(tr)) {
    x <- tr[[t]]
    if (is.null(x)) next  # 핵심: NULL step은 행을 만들지 않음

    idx <- idx + 1L

    out <- data.frame(
      nstep   = t,
      r       = if (t >= 1 && t <= nsteps) r_grid[t] else NA_real_,
      bw      = if (!is.null(x$bw)) x$bw else NA_real_,
      beta1   = if (!is.null(x$beta) && length(x$beta) >= 1) x$beta[1] else NA_real_,
      beta3   = if (!is.null(x$beta) && length(x$beta) >= 3) x$beta[3] else NA_real_,
      dQbeta1 = if (!is.null(x$dQbeta1)) x$dQbeta1 else NA_real_,
      dQbeta3 = if (!is.null(x$dQbeta3)) x$dQbeta3 else NA_real_
    )

    # pi
    for (k in 1:K) out[[paste0("pi", k)]] <- NA_real_
    if (!is.null(x$pi)) {
      pi <- x$pi
      if (length(pi) < K) pi <- c(pi, rep(NA_real_, K - length(pi)))
      for (k in 1:K) out[[paste0("pi", k)]] <- pi[k]
    }

    # lambda
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

trace_to_df_EM <- function(trace) {
  n <- length(trace)
  
  out <- data.frame(
    nstep   = seq_len(n),
    beta1   = NA_real_,
    beta2   = NA_real_,
    beta3   = NA_real_,
    dQbeta1 = NA_real_,
    dQbeta3 = NA_real_,
    pi1     = NA_real_,
    pi2     = NA_real_,
    pi3     = NA_real_,
    lambda1 = NA_real_,
    lambda2 = NA_real_,
    lambda3 = NA_real_
  )
  
  for (i in seq_len(n)) {
    ti <- trace[[i]]
    
    out$beta1[i]   <- ti$beta[1]
    out$beta2[i]   <- ti$beta[2]
    out$beta3[i]   <- ti$beta[3]
    out$dQbeta1[i] <- ti$dQbeta1
    out$dQbeta3[i] <- ti$dQbeta3
    
    out$pi1[i] <- ti$pi[1]
    out$pi2[i] <- ti$pi[2]
    out$pi3[i] <- ti$pi[3]
    
    out$lambda1[i] <- ti$lambda[1]
    out$lambda2[i] <- ti$lambda[2]
    out$lambda3[i] <- ti$lambda[3]
  }
  
  out
}

df_EM <- trace_to_df_EM(fit_EM$trace)
# df_DAEM <- trace_to_df(fit_DAEM,nsteps = nsteps,r_init=r_init, K = 3)
# df_BM <- trace_to_df(fit_BM,nsteps = nsteps,r_init=r_init, K = 3)
df_DHEM <- trace_to_df(fit_DHEM,nsteps = nsteps,r_init=r_init, K = 3)
df_adapDHEM <- trace_to_df(fit_adapDHEM,nsteps = nsteps,r_init=r_init, K = 3)




find_stationary <- function(df, method = "L1") {
  
  # 필수 컬럼 체크
  required_cols <- c("dQbeta1", "dQbeta3")
  if (!all(required_cols %in% names(df))) {
    stop("data.frame must contain dQbeta1 and dQbeta3 columns")
  }
  
  # norm 계산
  if (method == "L1") {
    score <- abs(df$dQbeta1) + abs(df$dQbeta3)
  } else if (method == "L2") {
    score <- sqrt(df$dQbeta1^2 + df$dQbeta3^2)
  } else {
    stop("method must be 'L1' or 'L2'")
  }
  
  # 최소값 행 반환
  stationary_row <- df[which.min(score), ]
  
  return(stationary_row)
}
find_stationary_adap<- function(df, tol = 1e-4) {
  
  score <- abs(df$dQbeta1) + abs(df$dQbeta3)
  
  candidates <- df[score <= tol, ]
  
  if (nrow(candidates) == 0) {
    return(NULL)
  }
  
  tail(candidates, 1)
}

expectedCompLog_Mstep <- function(df, pi, lambda, beta, r) {
  
  t     <- as.numeric(df$time)
  event <- as.numeric(df$event)
  n <- length(t)
  K <- length(pi)
  
  # annealed E-step (gamma 생성에만 r 사용)
  gamma <- weibull_estep_annealed(df, pi, lambda, beta, r)
  
  logt <- log(t)
  
  logLik_mat <- sapply(1:K, function(k) {
    event * (log(lambda[k]) + log(beta[k]) + (beta[k] - 1) * logt) -
      lambda[k] * (t ^ beta[k])
  })
  
  logpi_mat <- matrix(log(pi), nrow = n, ncol = K, byrow = TRUE)
  
  # 표준 Q (r 곱하지 않음)
  sum(gamma * (logpi_mat + logLik_mat))
}

df_DHEM$expectedLogLik <- mapply(
  FUN = function(r, beta1, beta3, pi1, pi2, pi3,
                 lambda1, lambda2, lambda3) {
    
    pi     <- c(pi1, pi2, pi3)
    lambda <- c(lambda1, lambda2, lambda3)
    beta   <- c(beta1, 1, beta3)   # beta2 = 1
    
    expectedCompLog_Mstep(
      df = df,
      pi = pi,
      lambda = lambda,
      beta = beta,
      r = r
    )
  },
  r       = df_DHEM$r,
  beta1   = df_DHEM$beta1,
  beta3   = df_DHEM$beta3,
  pi1     = df_DHEM$pi1,
  pi2     = df_DHEM$pi2,
  pi3     = df_DHEM$pi3,
  lambda1 = df_DHEM$lambda1,
  lambda2 = df_DHEM$lambda2,
  lambda3 = df_DHEM$lambda3
)

df_adapDHEM$expectedLogLik <- mapply(
  FUN = function(r, beta1, beta3, pi1, pi2, pi3,
                 lambda1, lambda2, lambda3) {
    
    pi     <- c(pi1, pi2, pi3)
    lambda <- c(lambda1, lambda2, lambda3)
    beta   <- c(beta1, 1, beta3)   # beta2 = 1
    
    expectedCompLog_Mstep(
      df = df,
      pi = pi,
      lambda = lambda,
      beta = beta,
      r = r
    )
  },
  r       = df_adapDHEM$r,
  beta1   = df_adapDHEM$beta1,
  beta3   = df_adapDHEM$beta3,
  pi1     = df_adapDHEM$pi1,
  pi2     = df_adapDHEM$pi2,
  pi3     = df_adapDHEM$pi3,
  lambda1 = df_adapDHEM$lambda1,
  lambda2 = df_adapDHEM$lambda2,
  lambda3 = df_adapDHEM$lambda3
)

df_EM$expectedLogLik <- mapply(
  FUN = function(r, beta1, beta3, pi1, pi2, pi3,
                 lambda1, lambda2, lambda3) {
    
    pi     <- c(pi1, pi2, pi3)
    lambda <- c(lambda1, lambda2, lambda3)
    beta   <- c(beta1, 1, beta3)   # beta2 = 1
    
    expectedCompLog_Mstep(
      df = df,
      pi = pi,
      lambda = lambda,
      beta = beta,
      r = 1
    )
  },
  beta1   = df_EM$beta1,
  beta3   = df_EM$beta3,
  pi1     = df_EM$pi1,
  pi2     = df_EM$pi2,
  pi3     = df_EM$pi3,
  lambda1 = df_EM$lambda1,
  lambda2 = df_EM$lambda2,
  lambda3 = df_EM$lambda3
)


cov_beta = function(df, title = NULL,
                    vline_at=NULL,
         axis_text_y_size  = 12,
         axis_title_y_size = 10,
         axis_text_x_size  = 10,
         axis_title_x_size = 10,
         title_size        = 30
          ){
  if (is.null(vline_at)) {
    vline_at <- find_stationary(df,method = "L2")$r
  }
  
              base_theme <- theme_bw() + 
                theme(
                  axis.text.y  = element_text(size = axis_text_y_size),
                  axis.title.y = element_text(size = axis_title_y_size,angle = 0),
                  axis.text.x  = element_text(size = axis_text_x_size),
                  axis.title.x = element_text(size = axis_title_x_size)
                )
              p1 = df %>% ggplot(aes(x=r,y=beta1))+geom_point()+geom_line()+
                geom_vline(xintercept = vline_at, linetype = "dashed", color = "red") +
                labs(x = "Annealing parameter", y ="")+ggtitle( expression(Trace~of~beta[1]))+
                coord_cartesian(xlim = c(0,1),ylim = c(0,1))+
                base_theme
              p3 = df %>% ggplot(aes(x=r,y=beta3))+geom_point()+geom_line()+
                geom_vline(xintercept = vline_at, linetype = "dashed", color = "red") +
                labs(x = "Annealing parameter", y ="")+ggtitle( expression(Trace~of~beta[3]))+
                coord_cartesian(xlim = c(0,1))+
                base_theme
              p1+p3
              
}

cov_diffbeta = function(df, title = NULL,
                        vline_at=NULL,
                        axis_text_y_size  = 12,
                        axis_title_y_size = 10,
                        axis_text_x_size  = 10,
                        axis_title_x_size = 10,
                        title_size        = 30
){
  if (is.null(vline_at)) {
    vline_at <-find_stationary(df,method = "L2")$r
  }
  base_theme <- theme_bw() + 
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size,angle = 0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )
  p1 = df %>% ggplot(aes(x=r,y=dQbeta1))+geom_point()+geom_line()+
    geom_vline(xintercept = vline_at, linetype = "dashed", color = "red") +
    labs(x = "Annealing parameter", y ="")+ggtitle( expression(Trace~of~nabla~beta[1]))+
    coord_cartesian(xlim = c(0,1))+
    base_theme
  p3 = df %>% ggplot(aes(x=r,y=dQbeta3))+geom_point()+geom_line()+
    geom_vline(xintercept = vline_at, linetype = "dashed", color = "red") +
    labs(x = "Annealing parameter", y ="")+ggtitle( expression(Trace~of~nabla~beta[3]))+
    coord_cartesian(xlim = c(0,1))+
    base_theme
  p1+p3
  
}

cov_pi <- function(df,
                   vline_at = NULL,
                   title = "Convergence of π") {
  
  if (is.null(vline_at)) vline_at <- max(df$r)
  
  dlong <- df |>
    dplyr::select(r, pi1, pi2, pi3) |>
    tidyr::pivot_longer(cols = c(pi1, pi2, pi3),
                        names_to = "component",
                        values_to = "pi")
  
  p <- ggplot(dlong, aes(x = r, y = pi, color = component)) +
    geom_line(linewidth = 1) +
    geom_vline(xintercept = vline_at, linetype = "dashed", color = "red") +
    labs(x = "Annealing parameter", y ="")+ggtitle( expression(Trace~of~pi))+
    coord_cartesian(xlim = c(0, 1)) +
    theme_bw() 
  
  p
}

# cov_pi(df_DHEM,vline_at=find_stationary(df_DHEM)$r )
plot_expectedLogLik <- function(df, vline_at = NULL) {
  
  if (is.null(vline_at)) {
    vline_at <- find_stationary(df)$r
  }
  
  df %>%
    ggplot(aes(x = r, y = expectedLogLik)) +
    geom_point() +
    geom_vline(xintercept = vline_at,
               linetype = "dashed",
               color = "red") +
    labs(x = "Annealing parameter", y ="")+ggtitle( expression(Trace~of~Q))+
    theme_bw()
}



cov_beta(df_DHEM,vline_at=find_stationary(df_DHEM)$r)/cov_diffbeta(df_DHEM,vline_at=find_stationary(df_DHEM)$r)+
  plot_annotation(
    title = expression(DHEM ~ "– DeviceG Data")
  )


cov_beta(df_adapDHEM,vline_at=find_stationary_adap(df_adapDHEM)$r)/cov_diffbeta(df_adapDHEM,vline_at=find_stationary_adap(df_adapDHEM)$r)+
  plot_annotation(
    title = expression(adapDHEM ~ "– DeviceG Data" ))



make_stationary_table <- function(df_DHEM, df_adapDHEM) {
  
  dhe <- find_stationary(df_DHEM)
  adap <- find_stationary_adap(df_adapDHEM)
  
  dhe$Algorithm  <- "DHEM"
  adap$Algorithm <- "Adaptive DHEM"
  
  # 열 순서 정리 (Algorithm 맨 앞으로)
  dhe  <- dhe[, c("Algorithm", names(dhe)[names(dhe) != "Algorithm"])]
  adap <- adap[, c("Algorithm", names(adap)[names(adap) != "Algorithm"])]
  
  rbind(dhe, adap)
}

Q_weibull <- function(df, pi, lambda, beta, r = 1) {
  
  t     <- as.numeric(df$time)
  event <- as.numeric(df$event)
  
  n <- length(t)
  K <- length(pi)
  
  # latent responsibilities (annealed)
  gamma <- weibull_estep_annealed(
    df = df,
    pi = pi,
    lambda = lambda,
    beta = beta,
    r = r
  )
  
  logt <- log(t)
  
  # log f_k(t_i | theta)
  logf <- sapply(1:K, function(k) {
    log(lambda[k]) +
      log(beta[k]) +
      (beta[k] - 1) * logt -
      lambda[k] * (t ^ beta[k])
  })
  
  logpi_mat <- matrix(log(pi), nrow = n, ncol = K, byrow = TRUE)
  
  # Q(theta | theta0)
  Q_val <- sum(gamma * (logpi_mat + logf))
  
  return(Q_val)
}


stationary_table <- make_stationary_table(df_DHEM, df_adapDHEM)
stationary_table
df_EM %>% tail(1)




##########

AkRk <- function(t, betak, lambdak, lambda2, pik = 1, pi2 = 1) {
  
  Ak <- (pik * betak * lambdak) / (pi2 * lambda2)
  
  Rk <- t^(betak - 1) * exp(-lambdak * t^betak + lambda2 * t)
  
  AkRk <- Ak * Rk
  
  return(AkRk)
}


# Aarest 

t <- seq(0.1, 86, length.out = 100)

y <- AkRk(
  t,
  betak   = 0.566,
  lambdak = 0.261,
  lambda2 = 0.0250,
  pik = 0.241, pi2 = 0.505
)

df1 <- data.frame(t = t, AkRk = y)

cross_time1 <- uniroot(
  function(t) AkRk(
    t,
    betak   = 0.566,
    lambdak = 0.261,
    lambda2 = 0.0250,
    pik = 0.241,
    pi2 = 0.505
  ) - 1,
  interval = c(0.001, 100)
)

P1 = df1 %>% ggplot(aes(t, AkRk)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = cross_time1$root, colour = "blue") +
  annotate("text",
           x = cross_time1$root+20,
           y = 3,
           label = paste0("Change point = ", round(cross_time1$root,2))) +
  annotate("text",
           x = 10,
           y = 1,
           label = "y = 1",
           vjust = -0.5) +
  labs(x = NULL,title="Burn-in vs Constant")+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 14))


# wearout 
t <- seq(0.1, 86, length.out = 100)

y <- AkRk(
  t,
  betak   = 78.1,
  lambdak = 2.46e-151,
  lambda2 = 0.0250,
  pik = 0.254, pi2 = 0.505
)

df3 <- data.frame(t = t, AkRk = y)

cross_time3 <- uniroot(
  function(t) AkRk(
    t,
    betak   = 78.1,
    lambdak = 2.46e-151,
    lambda2 = 0.0250,
    pik = 0.254, pi2 = 0.505
  ) - 1,
  interval = c(1, 80)
)



P3 = df3 %>% ggplot( aes(t, AkRk)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = cross_time3$root, colour = "blue") +
  annotate("text",
           x = cross_time3$root-20,
           y = 10,
           label = paste0("Change point = ", round(cross_time3$root,2))) +
  annotate("text",
           x = 1,
           y = 1,
           label = "y = 1",
           vjust = -0.5) +
  labs(x = NULL,title="Constant vs Wear-out")+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 14))



(P1+P3)+plot_annotation(title = "Aarest : Posterior Dominance and Change Points")

############
# FRT

t <- seq(0.1,300, length.out = 100)

y <- AkRk(
  t,
  betak   = 0.740,
  lambdak = 0.0229 ,
  lambda2 = 0.00452,
  pik =0.404, pi2 = 0.351
)

df1 <- data.frame(t = t, AkRk = y)

cross_time1 <- uniroot(
  function(t) AkRk(
    t,
    betak   = 0.740,
    lambdak = 0.0229 ,
    lambda2 = 0.00452,
    pik =0.404, pi2 = 0.351
  ) - 1,
  interval = c(0.001, 300)
)

P1 = df1 %>% ggplot(aes(t, AkRk)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = cross_time1$root, colour = "blue") +
  annotate("text",
           x = cross_time1$root+80,
           y = 2,
           label = paste0("Change point = ", round(cross_time1$root,2))) +
  annotate("text",
           x = 10,
           y = 1,
           label = "y = 1",
           vjust = -0.5) +
  labs(x = NULL,title="Burn-in vs Constant")+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 14))



# wearout 
t <- seq(0.1, 300, length.out = 100)

y <- AkRk(
  t,
  betak   = 6.63,
  lambdak = 3.97e-17,
  lambda2 = 0.00452,
  pik = 0.245, pi2 = 0.351
)

df3 <- data.frame(t = t, AkRk = y)

cross_time3 <- uniroot(
  function(t) AkRk(
    t,
    betak   = 6.63,
    lambdak = 3.97e-17,
    lambda2 = 0.00452,
    pik = 0.245, pi2 = 0.351
  ) - 1,
  interval = c(1, 300)
)



P3 = df3 %>% ggplot( aes(t, AkRk)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = cross_time3$root, colour = "blue") +
  annotate("text",
           x = cross_time3$root-100,
           y =2,
           label = paste0("Change point = ", round(cross_time3$root,2))) +
  annotate("text",
           x = 10,
           y = 1,
           label = "y = 1",
           vjust = -0.5) +
  labs(x = NULL,title="Constant vs Wear-out")+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 14))



(P1+P3)+plot_annotation(title = "Device G : Posterior Dominance and Change Points")




###########
# LFP 

# FRT

t <- seq(0.1,1370, length.out = 100)

y <- AkRk(
  t,
  betak   = 0.2006064 ,
  lambdak = 3.168796e-03 ,
  lambda2 = 2.929778e-16,
  pik =0.5032814 , pi2 = 0.4967186
)

df1 <- data.frame(t = t, AkRk = y)

cross_time1 <- uniroot(
  function(t) AkRk(
    t,
    betak   = 0.2006064 ,
    lambdak = 3.168796e-03 ,
    lambda2 = 2.929778e-16,
    pik =0.5032814 , pi2 = 0.4967186
  ) - 1,
  interval = c(0.00000001, 15000000000000000)
)


AkRk(
  15000000000000000,
  betak   = 0.2006064 ,
  lambdak = 3.168796e-03 ,
  lambda2 = 2.929778e-16,
  pik =0.5032814 , pi2 = 0.4967186
)


P1 = df1 %>% ggplot(aes(t, AkRk)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = cross_time1$root, colour = "blue") +
  annotate("text",
           x = cross_time1$root+80,
           y = 2,
           label = paste0("Change point = ", round(cross_time1$root,2))) +
  annotate("text",
           x = 10,
           y = 1,
           label = "y = 1",
           vjust = -0.5) +
  labs(x = NULL,title="Burn-in vs Constant")+
  theme_minimal()+
  theme(axis.text.x = element_text(size = 14))

P1




#############


compute_empirical_hazard_surv <- function(df) {
  stopifnot(all(c("time", "event") %in% names(df)))
  library(survival)
  
  # Surv 객체 생성
  surv_obj <- Surv(time = df$time, event = df$event)
  fit <- survfit(surv_obj ~ 1)
  
  # 시간과 누적 생존율
  times <- fit$time
  surv_probs <- fit$surv
  
  # 누적 hazard (대략적 추정)
  cumhaz <- -log(surv_probs)
  
  # 시간 구간별 변화량
  delta_time <- diff(c(0, times))
  delta_hazard <- diff(c(0, cumhaz))
  hazard_rate <- delta_hazard / delta_time 
  
  return(data.frame(time=times,hazard_rate=hazard_rate))
}


haz_df <- compute_empirical_hazard_surv(df)

ggplot(haz_df, aes(x = time, y = hazard_rate)) +
  geom_point(size = 2) +
  geom_line() +
  labs(
    x = "Time",
    y = "Empirical Hazard",
    title = "Empirical Hazard Estimatel"
  ) +
  theme_minimal()+
  theme(
  axis.title.x = element_text(size = 16),
  axis.title.y = element_text(size = 16),
  axis.text.x  = element_text(size = 14),
  axis.text.y  = element_text(size = 14)
)+
  labs(x = NULL)
