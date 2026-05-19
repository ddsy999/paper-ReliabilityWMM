library(ggplot2)
library(patchwork)
library(dplyr)

source("wmm_functions.R")

#### ---------------------------------------------------------------------------
# Data Load
#### ---------------------------------------------------------------------------
df <- read.table("DATA\\Aarest_data.txt", header = TRUE)
# df <- read.table("DATA\\FRT_censord.txt", header = TRUE)
# df <- read.table("DATA\\LFP.txt", header = TRUE)
# df <- read.table("DATA\\SerumReversal.txt", header = TRUE)

#### ---------------------------------------------------------------------------
# Hyper-parameter
#### ---------------------------------------------------------------------------
maxGEMiter <- 1e+3
errtol     <- 1e-9
bw_fixed   <- 1e-5  # Figure 2 오른쪽 열: 고정 barrier

#### ---------------------------------------------------------------------------
# Init-parameter
#### ---------------------------------------------------------------------------
K           <- 3
pi_init     <- rep(1 / K, K)
beta_init   <- c(0.5, 1, 2)
lambda_init <- wmm_lambda_init(df$time, df$event, beta_init,
                                ratio1 = 0.3, ratio3 = 0.7)
theta_init  <- list(beta = beta_init, pi = pi_init, lambda = lambda_init)

#### ---------------------------------------------------------------------------
# Fixed-(r, bw) EM trace
# - annealing parameter r 고정, barrier parameter bw 고정
# - pi, lambda, beta 만 매 iteration 업데이트
# - 매 iteration의 gamma, beta 기록
#### ---------------------------------------------------------------------------
run_fixed_trace <- function(df, theta, r, bw, maxiter = 1e6, tol = 1e-10) {
  pi     <- theta$pi
  lambda <- theta$lambda
  beta   <- theta$beta
  K      <- length(pi)
  N      <- nrow(df)
  t_vec  <- df$time
  e_vec  <- df$event

  gamma <- weibull_estep_annealed(df, pi, lambda, beta, r = r)

  theta_df           <- NULL
  result_latentZ_mat <- list()

  for (iter in seq_len(maxiter)) {
    new_pi    <- colSums(gamma) / N
    new_beta1 <- barrier_safe_wrapper1(beta[1], e_vec, t_vec, gamma, bw = bw)
    new_beta3 <- barrier_safe_wrapper3(beta[3], e_vec, t_vec, gamma, bw = bw)
    new_beta  <- c(new_beta1, 1, new_beta3)
    new_lam   <- sapply(seq_len(K), function(k)
      sum(gamma[, k] * e_vec) / sum(gamma[, k] * (t_vec ^ new_beta[k])))

    param_diff <- sqrt(sum((beta - new_beta)^2))
    beta <- new_beta; pi <- new_pi; lambda <- new_lam

    theta_df <- rbind(theta_df,
                      data.frame(iter  = iter,
                                 beta1 = beta[1],
                                 beta3 = beta[3]))
    result_latentZ_mat[[iter]] <- gamma

    if (param_diff < tol || iter == maxiter) break

    gamma <- weibull_estep_annealed(df, pi, lambda, beta, r = r)
  }

  list(theta_df = theta_df, result_latentZ_mat = result_latentZ_mat)
}

#### ---------------------------------------------------------------------------
# Run
#### ---------------------------------------------------------------------------
cat("r=1,   bw=bw_fixed ...\n")
fit_r1_bwf <- run_fixed_trace(df, theta_init, r = 1,   bw = bw_fixed, maxiter = maxGEMiter, tol = errtol)
cat("  iters:", length(fit_r1_bwf$result_latentZ_mat), "\n")

cat("r=0.8, bw=bw_fixed ...\n")
fit_r08_bwf <- run_fixed_trace(df, theta_init, r = 0.8, bw = bw_fixed, maxiter = maxGEMiter, tol = errtol)
cat("  iters:", length(fit_r08_bwf$result_latentZ_mat), "\n")

#### ---------------------------------------------------------------------------
# Plot builder
#### ---------------------------------------------------------------------------
make_column <- function(theta_df, result_latentZ_mat, col_title,
                        thin_to = 200) {
  n_iter <- length(result_latentZ_mat)
  n_obs  <- nrow(result_latentZ_mat[[1]])

  ratio_vec          <- sapply(result_latentZ_mat, function(z) z[1, 1] / z[1, 3])
  theta_df$log_ratio <- log10(ratio_vec)

  keep_idx <- if (n_iter > thin_to) unique(round(seq(1, n_iter, length.out = thin_to))) else seq_len(n_iter)
  sph_df <- bind_rows(lapply(keep_idx, function(i) {
    z <- result_latentZ_mat[[i]]
    data.frame(iter = i, index = seq_len(nrow(z)), z_i1 = z[, 1])
  }))

  final_val <- min(result_latentZ_mat[[n_iter]][, 1])

  btheme <- theme_minimal(base_size = 9) +
    theme(panel.grid.minor  = element_blank(),
          plot.title        = element_text(hjust = 0.5, size = 10),
          legend.key.height = unit(1.5, "cm"))

  R1 <- ggplot(theta_df, aes(iter, log_ratio, color = iter)) +
    geom_point(size = 0.4) + geom_line() +
    scale_color_gradient(low = "gray25", high = "red") +
    labs(title = expression("Trend of " * log[10](Z[11]/Z[13])), y = "", x = "") +
    btheme + theme(legend.position = "none", axis.text.y = element_text(size = 14),
                   axis.text.x = element_text(size = 14),
                   plot.title = element_text(hjust = 0.5, size = 18))

  R2 <- ggplot(theta_df, aes(iter, beta1, color = iter)) +
    geom_point(size = 0.4) + geom_line() +
    geom_hline(yintercept = tail(theta_df$beta1, 1), linetype = "dashed") +
    scale_color_gradient(low = "gray25", high = "red") +
    labs(title = expression("Convergence of "*beta[1]), y = "", x = "") +
    btheme + theme(legend.position = "none", axis.text.y = element_text(size = 16),
                   axis.text.x = element_text(size = 14),
                   plot.title = element_text(hjust = 0.5, size = 18))

  R3 <- ggplot(sph_df, aes(index, z_i1, group = iter, color = iter)) +
    geom_line(alpha = 0.3) +
    scale_color_gradient(low = "gray25", high = "red", name = "Iteration") +
    geom_hline(yintercept = final_val, linetype = "dashed") +
    annotate("text", x = n_obs * 0.85, y = final_val,
             label = formatC(final_val, format = "g", digits = 5),
             vjust = -0.5, size = 5) +
    labs(title = expression("Latent variable " * Z[i1]),
         x = paste0("Index = 1~", n_obs), y = "") +
    btheme +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5, size = 18),
          axis.text.x = element_text(size = 14),
          legend.text  = element_text(size = 10),
          legend.title = element_text(size = 12))

  (R1 + R2 + R3) +
    plot_annotation(title = col_title,
                    theme = theme(plot.title = element_text(hjust = 0.5, size = 12))) +
    plot_layout(design = "AC\nBC", widths = c(1, 1.6))
}

#### ---------------------------------------------------------------------------
# Figure 조립
#### ---------------------------------------------------------------------------
col_r1_bwf  <- make_column(fit_r1_bwf$theta_df,  fit_r1_bwf$result_latentZ_mat,  "")
col_r08_bwf <- make_column(fit_r08_bwf$theta_df, fit_r08_bwf$result_latentZ_mat, "")

col_r1_bwf
col_r08_bwf
figure2 <- col_r1_bwf | col_r08_bwf
figure2

# ggsave("figure2.pdf", figure2, width = 14, height = 5.5)
