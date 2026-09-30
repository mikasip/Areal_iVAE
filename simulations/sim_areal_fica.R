args <- commandArgs(trailingOnly = TRUE)
seed <- as.integer(args[1])
n_layers <- as.integer(args[2])
nonlin <- args[3]
nonstat_time <- as.logical(as.integer(args[4]))   # controls AR-coefficient path AND block variance in time
nonstat_space <- as.logical(as.integer(args[5]))   # controls MEM-based variance nonstationarity in space
point_type <- args[6] # "clustered" or "uniform"
spat_model <- args[7]

type <- "areal_sar_fica"
filename <- paste0("res_", type, "_", n_layers, "_", nonlin, "_",
                    nonstat_time, "_", nonstat_space, "_", point_type, "_", seed, ".RData")

library(Matrix)
library(RcppHungarian)
library(NonlinearBSS)

source("helpers/compute_moran_eigenvectors.R")
source("helpers/helpers.R")
source("helpers/areal_data_generation.R")

set.seed(seed)

# ---------------------------------------------------------------
# 1. Spatial domain: point configuration + order-1 / order-2 adjacency
# ---------------------------------------------------------------
n_s <- 200
n_t <- 300
p   <- 5   # number of latent components
s   <- 7   # observed mixed dimension

n_fine <- 3000   # denser fine layer -> smoother, more organic area boundaries

fine <- generate_fine_layer(n_fine)
seed_idx <- place_seeds_fixed(fine$pts, n_s, seed_type = point_type)
cluster_id <- assign_graph_voronoi(fine$graph, seed_idx)

# renumber clusters 1..n_s_actual (some seeds can end up unreachable/duplicated)
cluster_id <- as.integer(factor(cluster_id))
n_s_actual <- length(unique(cluster_id))

W1 <- build_areal_adjacency(fine$graph, cluster_id)

W1_list <- mat_to_nblist(W1)

deg <- rowSums(W1)
cat(sprintf("point_type=%s  mean deg=%.2f  sd deg=%.2f\n", point_type, mean(deg), sd(deg)))

W2      <- order2_neighbors(W1, cumulative = TRUE)   # order-1 U order-2, ground-truth topology
W2_list <- mat_to_nblist(W2)

# ---------------------------------------------------------------
# 2. Ground-truth MEM basis (order-2) used to drive spatial variance
#    nonstationarity, and to build the "correct" aux_data.
#    A separate order-1-only MEM basis is built for the "mismatched"
#    aux_data comparison.
# ---------------------------------------------------------------
K_mem <- 20

mem_correct_obj <- compute_MEM_vectors(W2_list, K = K_mem, positive_only = TRUE)
mem_mismatch_obj <- compute_MEM_vectors(W1_list, K = K_mem, positive_only = TRUE)

mem_correct  <- mem_correct_obj$mem    # n_s x K_mem, ground-truth (order-2) basis
mem_mismatch <- mem_mismatch_obj$mem   # n_s x K_mem, order-1-only basis (deliberately mismatched)

# ---------------------------------------------------------------
# 3. Simulate p latent components over the (n_s x n_t) areal domain.
#    Fixed spatial coupling (beta_s) across ALL components/settings,
#    so topology-induced marginal-variance heterogeneity is constant
#    across conditions and cannot confound the nonstat comparisons.
# ---------------------------------------------------------------
sigma0   <- 1

simulate_one_component <- function(mem_vals_i, beta_s, baseline_ar1) {
  # spatial variance: constant, or MEM-driven (nonstat_space)
  if (spat_model == "SAR") {
  edge_weights_W2 <- lapply(W2_list, function(nb) {
    if (length(nb) == 0) numeric(0) else rep(1 / length(nb), length(nb))
  })
  row_i <- integer(0)
  col_j <- integer(0)
  vals <- numeric(0)
  for (i in seq_len(n_s)) {
      nb <- W2_list[[i]]
  if (length(nb) == 0) next
      w_i <- edge_weights_W2[[i]]
      row_i <- c(row_i, rep(i, length(nb)))
      col_j <- c(col_j, nb)
      vals <- c(vals, -beta_s * w_i)
  }
  Bmat       <- Matrix::sparseMatrix(row_i, col_j, x = vals, dims = c(n_s, n_s))
  SAR_system <- Matrix::Diagonal(n_s) + Bmat   # time-invariant: fixed beta_s
  } else {
    R_CAR <- build_CAR_chol(as(W1, "sparseMatrix"), alpha = beta_s)  # built once, outside the per-component loop
  }


  mem_vals_cent <- mem_vals_i - min(mem_vals_i)
  mem_vals_scaled <- (mem_vals_cent) / max(mem_vals_cent) # [0, 1] range
  spat_var_min <- sqrt(runif(1, 0.1, 0.5))
  spat_var_max <- sqrt(runif(1, 3, 10))
  spat_var_range <- spat_var_max - spat_var_min
  mem_vals_scaled <- mem_vals_scaled * spat_var_range + spat_var_min
  sigma_i <- if (!nonstat_space) rep(sigma0, n_s) else sigma0 * mem_vals_scaled

  # temporal behaviour: stationary, or smooth-AR + block-variance combined (nonstat_time)
  if (!nonstat_time) {
    beta_t_t  <- rep(baseline_ar1, n_t)
    time_mult <- rep(1, n_t)
  } else {
    ar_shift <- runif(1, 0, 2 * pi)
    ar_scale <- runif(1, 1, 10)
    raw <- changing_AR1_coefs(n_t, shift = ar_shift, scale = ar_scale)
    beta_t_t  <- baseline_ar1 / 2 + baseline_ar1 / 2 * raw 
    time_mult <- block_temporal_variance(n_t, n_seg_time = 10, min_mult = 0.1, max_mult = 10)
  }

  z_mat <- matrix(0, n_s, n_t)
  for (t in seq_len(n_t)) {
    sigma_it <- sigma_i * sqrt(time_mult[t])
    eps_t <- rnorm(n_s, mean = 0, sd = sigma_it)
    if (spat_model == "SAR") {
      u_t <- as.numeric(Matrix::solve(SAR_system, eps_t))  # SAR structure applied ONLY to innovation
    } else {
      tau_it <- sigma0 * sqrt(time_mult[t])
      u_t <- sample_CAR_innovation(R_CAR, tau = tau_it, sigma_scale = sigma_i)
    }

    if (t == 1) {
      z_mat[, t] <- u_t
    } else {
      z_mat[, t] <- beta_t_t[t] * z_mat[, t - 1] + u_t     # plain AR(1) recursion, no solve() here
    }
  }
  as.numeric(z_mat)
}

n <- n_s * n_t
data <- array(NA, c(n, p))
for (i in seq_len(p)) {
  # cycle through the leading ground-truth MEM columns for spatial-variance
  # diversity across components (only used when nonstat_space = TRUE)
  mem <- mem_correct[, sample(1:K_mem, 1)]
  beta_s_i <- runif(1, 0.3, 0.9)
  ar1_coef_baseline_i <- runif(1, 0.3, 0.9)
  print(paste("Component:", i))
  print(paste("Beta:", beta_s_i))
  print(paste("AR:", ar1_coef_baseline_i))
  data[, i] <- simulate_one_component(mem, beta_s_i, ar1_coef_baseline_i)
}

data_scaled   <- data / mean(sqrt(diag(var(data))))
obs_data_full <- mix_data_over_determinated(data_scaled, s, n_layers, nonlin)
obs_data <- obs_data_full   # no missing values in this script

res_df <- data.frame(matrix(ncol = 4, nrow = 1))
ind <- 1

fica_res <- fastICA::fastICA(as.matrix(obs_data), p)
MCCfica <- absolute_mean_correlation(cor(fica_res$S, data))

res_df[ind, ] <- c("fastICA", MCCfica, n_layers, seed)

colnames(res_df) <- c("method", "MCC", "n_layers", "seed")
res_df$nonstat_time  <- nonstat_time
res_df$nonstat_space <- nonstat_space
res_df$point_type    <- point_type
res_df$nonlin        <- nonlin

res_test <- res_df

res_dir <- "results"
if (!dir.exists(res_dir)) dir.create(res_dir)
save(res_test, file = paste0(res_dir, "/", filename))
res_test