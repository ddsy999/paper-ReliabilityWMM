library(clue)
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(patchwork)
library(purrr)

source("wmm_functions.R")



#### ---------------------------------------------------------------------------
# Data Load
#### ---------------------------------------------------------------------------
df <- read.table("DATA\\Aarest_data.txt", header = TRUE)
#df <- read.table("DATA\\FRT_censord.txt", header = TRUE)
#df <- read.table("DATA\\LFP.txt", header = TRUE)
#df <- read.table("DATA\\SerumReversal.txt", header = TRUE)


#### ---------------------------------------------------------------------------
# Hyper-parameter (set here first)
#### ---------------------------------------------------------------------------
maxGEMiter = 1e+6
nsteps     <- 100
r_init     <- 0.1
r_end      <- 1
#bw_init    <- 1e-3
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
# Init barrier-parameter
#### ---------------------------------------------------------------------------

bw_init_wmm <- wmm_bw_init_beta(
  df = df,
  pi_init = pi_init,
  lambda_init = lambda_init,
  beta_init = beta_init,
  r_init = r_init,
  tau = 0.1
)

bw_init=bw_init_wmm

#### ---------------------------------------------------------------------------
# Train (one run each)
#### ---------------------------------------------------------------------------
verbose=TRUE
# 1) EM (standard EM: r=1, bw=0)  -- wmm_EM은 내부에서 wmm_em_at_r_bw 사용
fit_EM <- wmm_EM(df,theta=theta_init,maxGEMiter = maxGEMiter,tol=errtol,verbose=verbose)
print("end EM")
# 2) DAEM (bw=0, r schedule)
#fit_DAEM <- wmm_DAEM(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps = nsteps,r_init = r_init,r_end=r_end,tol=errtol,verbose=verbose)
#print("end DAEM")
# 3) Barrier method (r=1, bw schedule)
#fit_BM <- wmm_BM(df,theta=theta_init,maxGEMiter = maxGEMiter,nsteps=nsteps,bw_init = bw_init,bw_end = bw_end,tol = errtol,verbose=verbose)
#print("end BM")
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


df_EM <- trace_to_df_EM(fit_EM$trace)
#df_DAEM <- trace_to_df(fit_DAEM,nsteps = nsteps,r_init=r_init, K = 3)
#df_BM <- trace_to_df(fit_BM,nsteps = nsteps,r_init=r_init, K = 3)
df_DHEM <- trace_to_df(fit_DHEM,nsteps = nsteps,r_init=r_init, K = 3)
df_adapDHEM <- trace_to_df(fit_adapDHEM,nsteps = nsteps,r_init=r_init, K = 3)




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
      r = 1
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
      r = 1
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


cov_beta(df_DHEM,vline_at=find_stationary(df_DHEM)$r)/cov_diffbeta(df_DHEM,vline_at=find_stationary(df_DHEM)$r)+
  plot_annotation(
    title = expression(DHEM)
  )


cov_beta(df_adapDHEM,vline_at=find_stationary_adap(df_adapDHEM)$r)/cov_diffbeta(df_adapDHEM,vline_at=find_stationary_adap(df_adapDHEM)$r)+
  plot_annotation(
    title = expression(adapDHEM ))



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


stationary_table <- make_stationary_table(df_DHEM, df_adapDHEM)
stationary_table
tail(df_EM,1)







##########

AkRk <- function(t, betak, pik = 1, pi2 = 1, df = NULL) {
  lambdak <- sum(df$event) / sum(df$time ^ betak)
  lambda2 <- sum(df$event) / sum(df$time)

  Ak <- (pik * betak * lambdak) / (pi2 * lambda2)

  Rk <- t^(betak - 1) * exp(-lambdak * t^betak + lambda2 * t)

  AkRk <- Ak * Rk

  return(AkRk)
}


# Change point

dataName = "Adaptive DHEM"

t <- seq(0.1, max(df$time), length.out = 100)

y <- AkRk(
  t,
  betak   = stationary_table |> filter(Algorithm==dataName) |> pull(beta1),
  pik = stationary_table |> filter(Algorithm==dataName) |> pull(pi1),
  pi2 = stationary_table |> filter(Algorithm==dataName) |> pull(pi2), 
  df = df
)

df1 <- data.frame(t = t, AkRk = y)

cross_time1 <- uniroot(
  function(t) AkRk(
  t,
  betak   = stationary_table |> filter(Algorithm==dataName) |> pull(beta1),
  pik = stationary_table |> filter(Algorithm==dataName) |> pull(pi1),
  pi2 = stationary_table |> filter(Algorithm==dataName) |> pull(pi2), 
  df = df
  ) - 1,
  interval = c(0.001, max(df$time)-5)
)

P1 = df1 %>% ggplot(aes(t, AkRk)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = cross_time1$root, colour = "blue") +
  annotate("text",
           x = cross_time1$root+35,
           y = 1.5,size=5,
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
t <- seq(0.1,  max(df$time), length.out = 100)

y <- AkRk(
  t,
  betak   = stationary_table |> filter(Algorithm==dataName) |> pull(beta3),
  pik = stationary_table |> filter(Algorithm==dataName) |> pull(pi3),
  pi2 = stationary_table |> filter(Algorithm==dataName) |> pull(pi2), 
  df = df
)

df3 <- data.frame(t = t, AkRk = y)

cross_time3 <- uniroot(
  function(t) AkRk(
    t,
  betak   = stationary_table |> filter(Algorithm==dataName) |> pull(beta3),
  pik = stationary_table |> filter(Algorithm==dataName) |> pull(pi3),
  pi2 = stationary_table |> filter(Algorithm==dataName) |> pull(pi2), 
  df = df
  ) - 1,
  interval = c(1, max(df$time)-5)
)



P3 = df3 %>% ggplot( aes(t, AkRk)) +
  geom_line() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_vline(xintercept = cross_time3$root, colour = "blue") +
  annotate("text",
           x = cross_time3$root-30,
           y = 5,size=5,
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
