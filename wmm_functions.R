library(survival)

diffB_onlyB = function(beta,event_vec,time_vec,latentZ_mat,j){
  sum(latentZ_mat[,j]*event_vec)/beta + 
    sum(latentZ_mat[,j]*event_vec*log(time_vec))-
    sum(latentZ_mat[,j]*event_vec)*sum(latentZ_mat[,j]*(time_vec^beta)*log(time_vec))/sum(latentZ_mat[,j]*(time_vec^beta))
}

barrierFunc_1 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  result =  diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=1)+(1/beta -1/(1-beta))*(bw)
  return(result)
}

barrierFunc_3 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  result =  diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=3)+bw*(1/(beta-1))
  return(result)
}

barrier_beta1 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  isna_diffbeta = function(beta) is.na(diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=1))

  if(isna_diffbeta(beta)){
    maxRange = 1e-12
  }else{
    maxRange = beta
  }
  eps = 1e-12
  if(bw==0){
    while(!isna_diffbeta(maxRange)){
      maxRange=maxRange + min(maxRange*1.01,10)
      if(diffB_onlyB(maxRange,event_vec,time_vec,latentZ_mat,j=1)*
      diffB_onlyB(1e-12,event_vec,time_vec,latentZ_mat, j=1)<0) break
    }
      result <- uniroot(function(beta) diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=1),
      interval = c(0,maxRange),tol=1e-10)
      return(result$root)
  }

  result <- uniroot(function(beta) barrierFunc_1(beta,event_vec,time_vec,latentZ_mat,bw),
  interval = c(eps,1-eps),tol=1e-10)
  return(result$root)
}

barrier_beta3 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  isna_diffbeta = function(beta) is.na(diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=3))
  
  if(isna_diffbeta(beta)){
    maxRange = 1e-12
  }else{
    maxRange = beta
  }
  eps= 1e-3
  if(bw==0){
    while(!isna_diffbeta(maxRange)&&!isna_diffbeta(eps)){
      maxRange=maxRange*1.01
      if(diffB_onlyB(maxRange,event_vec,time_vec,latentZ_mat,j=3)*
      diffB_onlyB(eps,event_vec,time_vec,latentZ_mat, j=3)<0) break
    }
    result = uniroot(function(beta) diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=3),
    interval = c(eps, maxRange),tol=1e-10)
    return(result$root)
  }

  while(!isna_diffbeta(maxRange)){
    maxRange=maxRange + min(maxRange*1.01,10)
    if(diffB_onlyB(maxRange,event_vec,time_vec,latentZ_mat,j=3)*
    diffB_onlyB(1,event_vec,time_vec,latentZ_mat, j=3)<0) break
  }
  result = uniroot(function(beta) barrierFunc_3(beta,event_vec,time_vec,latentZ_mat,bw),
  interval = c(1, maxRange),tol=1e-10)
  return(result$root)
}


barrier_beta2 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  isna_diffbeta = function(beta) is.na(diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=2))
  
  if(isna_diffbeta(beta)){
    maxRange = 1e-12
  }else{
    maxRange = beta
  }
  eps= 1e-3

  while(!isna_diffbeta(maxRange)&&!isna_diffbeta(eps)){
    maxRange=maxRange*1.01
    if(diffB_onlyB(maxRange,event_vec,time_vec,latentZ_mat,j=2)*
       diffB_onlyB(eps,event_vec,time_vec,latentZ_mat, j=2)<0) break
  }
  result = uniroot(function(beta) diffB_onlyB(beta,event_vec,time_vec,latentZ_mat, j=2),
                   interval = c(eps, maxRange),tol=1e-10)
  return(result$root)
  
}


barrier_safe_wrapper1 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  tryCatch(barrier_beta1(beta,event_vec,time_vec,latentZ_mat,bw),
   error = function(e) {beta})
}

barrier_safe_wrapper3 = function(beta,event_vec,time_vec,latentZ_mat,bw){
  tryCatch(barrier_beta3(beta,event_vec,time_vec,latentZ_mat,bw),
   error = function(e) {beta})
}

weibull_estep_annealed <- function(df, pi, lambda, beta, r = 1) {
  t     <- as.numeric(df$time)
  event <- as.numeric(df$event)
  
  n <- length(t)
  K <- length(pi)
  
  logt <- log(t)
  log_resp <- matrix(NA_real_, nrow = n, ncol = K)
  
  for (k in 1:K) {
    # Weibull component k:
    # f_k(t) = lambda_k * beta_k * t^(beta_k-1) * exp(-lambda_k * t^beta_k)
    # S_k(t) = exp(-lambda_k * t^beta_k)
    #
    # log L_ik = event_i * log f_k(t_i) + (1-event_i) * log S_k(t_i)
    #         = event_i*(log lambda_k + log beta_k + (beta_k-1)log t_i) - lambda_k t_i^beta_k
    logLik_ik <- event * (log(lambda[k]) + log(beta[k]) + (beta[k] - 1) * logt) -
      lambda[k] * (t ^ beta[k])
    
    # annealed responsibilities: gamma_ik ∝ exp( r * (log pi_k + logLik_ik) )
    log_resp[, k] <- r * (log(pi[k]) + logLik_ik)
  }
  
  # numeric stabilization and normalization across k
  rowmax <- apply(log_resp, 1, max)
  w <- exp(log_resp - rowmax)
  gamma <- w / rowSums(w)
  
  gamma
}

wmm_delta_DKL <- function(df, theta0, theta1, r) {

  # standard posterior (r=1) for theta0, theta1
  gamma0 <- weibull_estep_annealed(df, theta0$pi, theta0$lambda, theta0$beta, r = 1)
  gamma1 <- weibull_estep_annealed(df, theta1$pi, theta1$lambda, theta1$beta, r = 1)

  # annealed posterior weight (r=r) under theta0
  w <- weibull_estep_annealed(df, theta0$pi, theta0$lambda, theta0$beta, r = r)

  # per-i contribution: sum_k w_ik * (log gamma0_ik - log gamma1_ik)
  contrib <- rowSums(w * (log(gamma0) - log(gamma1)))

  mean(contrib)
}

wmm_DKL <- function(df, theta0, theta1) {
  g0 <- weibull_estep_annealed(df, theta0$pi, theta0$lambda, theta0$beta, r = 1)
  g1 <- weibull_estep_annealed(df, theta1$pi, theta1$lambda, theta1$beta, r = 1)
  mean(rowSums(g0 * (log(g0) - log(g1))))
}

wmm_bar_value <- function(beta) {
  log(beta[1]) + log(1 - beta[1]) + log(beta[3] - 1)
}

wmm_bar_diff <- function(theta_new, theta_old) {
  wmm_bar_value(theta_new$beta) - wmm_bar_value(theta_old$beta)
}

wmm_bar_value2 <- function(beta) {
  log(beta[1]) + log(1 - beta[1]) 
}

wmm_bar_diff2 <- function(theta_new, theta_old) {
  wmm_bar_value2(theta_new$beta) - wmm_bar_value2(theta_old$beta)
}

wmm_lambda_init <-function(time_vec,event_vec,beta_vec,ratio1,ratio3){
  surv_obj <- Surv(time = time_vec, event = event_vec)
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

  df_haz = data.frame(time=unique(time_vec),hazard=hazard_rate,
                      time1= beta_vec[1] * unique(time_vec)^(beta_vec[1] - 1),
                      time3= beta_vec[3] * unique(time_vec)^(beta_vec[3] - 1)
                      )
  df_haz[(nrow(df_haz)),"hazard"] = df_haz[(nrow(df_haz)-1),"hazard"]
  n = nrow(df_haz)  
  subset_df1 = df_haz %>% filter(time<max(times)*ratio1)
  subset_df3 = df_haz %>% filter(time>max(times)*ratio3)
  subset_df2 = df_haz %>% filter(time>max(times)*ratio1 & time<max(times)*ratio3 )
  fit1 <- lm(hazard ~ 0 + time1, data = subset_df1)  # '0 +'는 intercept 제외
  fit3 <- lm(hazard ~ 0 + time3, data = subset_df3)  # '0 +'는 intercept 제외
  
  # 결과 확인
  lambda_est1 <- coef(fit1)[1]
  lambda_est3 <- coef(fit3)[1]
  lambda_est2 = mean(subset_df2[,"hazard"])
  lambda_vec = c(lambda_est1,lambda_est2,lambda_est3)
  return(lambda_vec)
}

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
    #bw = bw*0.5
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
find_stationary_adap<- function(df, tol = 1e-3) {
  
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
    scale_y_continuous(labels = scales::scientific)+
    base_theme
  p3 = df %>% ggplot(aes(x=r,y=dQbeta3))+geom_point()+geom_line()+
    geom_vline(xintercept = vline_at, linetype = "dashed", color = "red") +
    labs(x = "Annealing parameter", y ="")+ggtitle( expression(Trace~of~nabla~beta[3]))+
    coord_cartesian(xlim = c(0,1))+
    scale_y_continuous(labels = scales::scientific)+
    base_theme
  p1+p3
  
}

cov_beta1 = function(df, title = NULL,
                     vline_at=NULL,
                     axis_text_y_size  = 12,
                     axis_title_y_size = 10,
                     axis_text_x_size  = 10,
                     axis_title_x_size = 10,
                     title_size        = 30
){
  if (is.null(vline_at)) {
    vline_at <- df$r[which.min(abs(df$dQbeta1))]
  }
  base_theme <- theme_bw() +
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size, angle = 0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )
  df %>% ggplot(aes(x=r, y=beta1)) + geom_point() + geom_line() +
    geom_vline(xintercept = vline_at, linetype = "dashed", color = "red") +
    labs(x = "Annealing parameter", y = "") + ggtitle(expression(Trace~of~beta[1])) +
    coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
    base_theme
}

cov_diffbeta1 = function(df, title = NULL,
                         vline_at=NULL,
                         axis_text_y_size  = 12,
                         axis_title_y_size = 10,
                         axis_text_x_size  = 10,
                         axis_title_x_size = 10,
                         title_size        = 30
){
  if (is.null(vline_at)) {
    vline_at <- df$r[which.min(abs(df$dQbeta1))]
  }
  base_theme <- theme_bw() +
    theme(
      axis.text.y  = element_text(size = axis_text_y_size),
      axis.title.y = element_text(size = axis_title_y_size, angle = 0),
      axis.text.x  = element_text(size = axis_text_x_size),
      axis.title.x = element_text(size = axis_title_x_size)
    )
  df %>% ggplot(aes(x=r, y=dQbeta1)) + geom_point() + geom_line() +
    geom_vline(xintercept = vline_at, linetype = "dashed", color = "red") +
    labs(x = "Annealing parameter", y = "") + ggtitle(expression(Trace~of~nabla~beta[1])) +
    coord_cartesian(xlim = c(0,1)) +
    scale_y_continuous(labels = scales::scientific) +
    base_theme
}

cov_pi <- function(df,
                   vline_at = NULL,
                   title = "Convergence of π"
                  ) {
  
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

wmm_bw_init_beta <- function(df, pi_init, lambda_init, beta_init, r_init,
                             tau = 0.1) {
  z_annealed <- weibull_estep_annealed(
    df, pi_init, lambda_init, beta_init, r_init
  )

  if(length(beta_init)==2){
    g1 <- diffB_onlyB(beta_init[1], df$event, df$time, z_annealed, j = 1)
    return(bw1 <- tau * abs(g1) * min(beta_init[1], 1 - beta_init[1]))
  }

  g1 <- diffB_onlyB(beta_init[1], df$event, df$time, z_annealed, j = 1)
  g3 <- diffB_onlyB(beta_init[3], df$event, df$time, z_annealed, j = 3)
  
  bw1 <- tau * abs(g1) * min(beta_init[1], 1 - beta_init[1])
  bw3 <- tau * abs(g3) * (beta_init[3] - 1)
  
  min(bw1, bw3)
}
