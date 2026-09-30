compute_edge_weights <- function(phi, adj_list, A, b, c) {
    n <- nrow(phi)
    weights <- vector("list", n)
    for (i in seq_len(n)) {
        neighbors <- adj_list[[i]]
        if (length(neighbors) == 0) { weights[[i]] <- numeric(0); next }
        ui <- phi[i, ]
        scores <- sapply(neighbors, function(j) {
            uj <- phi[j, ]
            # bilinear + receiver + sender
            as.numeric(t(ui) %*% A %*% uj) + sum(b * ui) + sum(c * uj)
        })
        # softmax
        scores <- scores - max(scores)  # numerical stability
        weights[[i]] <- exp(scores) / sum(exp(scores))
    }
    return(weights)
}

modulate_by_eigenvectors <- function(phi_mat, n_vec, val_min, val_max) {
    K <- ncol(phi_mat)
    phi_sub <- sample(seq_len(K), n_vec, replace = FALSE)
    phi_weights <- runif(n_vec, -1, 1)
    # Compute linear combination of selected eigenvectors
    val <- phi_mat[, phi_sub] %*% phi_weights
    # Scale to [val_min, val_max]
    val_scaled <- (val - min(val)) / (max(val) - min(val)) * (val_max - val_min) + val_min
    return(val_scaled)
}

library(deldir)
make_adjacency_delaunay <- function(locations) {
    n <- nrow(locations)
    triangulation <- deldir(locations[, 1], locations[, 2])
    edges <- triangulation$delsgs  # each row: (x1,y1,x2,y2, node1, node2)

    # Build adjacency list (1-indexed, undirected)
    adj_list <- vector("list", n)
    for (i in seq_len(n)) adj_list[[i]] <- integer(0)

    for (e in seq_len(nrow(edges))) {
        node1 <- edges[e, 5]
        node2 <- edges[e, 6]
        adj_list[[node1]] <- c(adj_list[[node1]], node2)
        adj_list[[node2]] <- c(adj_list[[node2]], node1)
    }

    # Remove duplicate neighbors (deldir can produce them)
    adj_list <- lapply(adj_list, unique)
    return(adj_list)
}

generate_points <- function(n_s, type = c("uniform", "clustered"),
                             n_clusters = 3, cluster_sd = 0.04, xy_coef) {
  type <- match.arg(type)
  if (type == "uniform") {
    coords <- matrix(runif(2 * n_s), ncol = 2)
  } else {
    centers <- matrix(runif(2 * n_clusters), ncol = 2)
    cluster_id <- sample(seq_len(n_clusters), n_s, replace = TRUE)
    coords <- centers[cluster_id, ] + matrix(c(rnorm(n_s, 0, cluster_sd), rnorm(n_s, 0, xy_coef * cluster_sd)), ncol = 2)
    coords <- pmin(pmax(coords, 0), 1)  # keep inside unit square
  }
  coords
}

build_radius_graph <- function(coords, radius) {
  D <- as.matrix(dist(coords))
  W1 <- ((D > 0) & (D < radius)) * 1

  # guard against isolated nodes: connect to nearest point if degree 0
  deg <- rowSums(W1)
  for (i in which(deg == 0)) {
    nearest <- setdiff(order(D[i, ]), i)[1]
    W1[i, nearest] <- 1
    W1[nearest, i] <- 1
  }
  W1
}

mat_to_nblist <- function(W) {
    lapply(seq_len(nrow(W)), function(i) which(W[i, ] == 1))
}

order2_neighbors <- function(W1, cumulative = TRUE) {
  W2_reach <- (W1 %*% W1) > 0
  diag(W2_reach) <- FALSE
  if (cumulative) {
    Wout <- (W1 | W2_reach) * 1
  } else {
    Wout <- (W2_reach & !as.logical(W1)) * 1
  }
  diag(Wout) <- 0
  Wout
}

get_ground_truth_mem <- function(W1, mem_index = 1, cumulative = TRUE) {
  W2 <- order2_neighbors(W1, cumulative = cumulative)
  nb2_list <- mat_to_nblist(W2)
  mem_obj <- compute_MEM_vectors(nb2_list, K = max(mem_index, 5), positive_only = TRUE)
  mem_obj$mem[, mem_index]
}

# Block-wise temporal variance multiplier: divide time into n_seg_time
# segments (nearest-center assignment), unique variance multiplier per segment
block_temporal_variance <- function(n_t, n_seg_time, min_mult = 0.2, max_mult = 5) {
  cent_times  <- sample(seq_len(n_t), n_seg_time)
  labels_time <- sapply(seq_len(n_t), function(t) which.min(abs(cent_times - t)))
  seg_mult    <- runif(n_seg_time, min_mult, max_mult)
  seg_mult[labels_time]   # length n_t, one multiplier per time point
}

changing_AR1_coefs <- function(n, shift, scale) {
  cos((2 * pi * seq_len(n) * scale) / n - shift)
}

# ---- CAR precision + Cholesky, built once per replicate ----
build_CAR_chol <- function(W_bin, alpha) {
  deg <- Matrix::rowSums(W_bin)
  Dmat <- Matrix::Diagonal(x = deg)
  Q0 <- Dmat - alpha * W_bin          # precision up to 1/tau^2 scaling
  Matrix::chol(Q0)                    # sparse upper-triangular R, t(R) %*% R = Q0
}

sample_CAR_innovation <- function(R_chol, tau, sigma_scale = NULL) {
  n <- nrow(R_chol)
  z <- rnorm(n)
  u <- tau * as.numeric(Matrix::solve(R_chol, Matrix::solve(Matrix::t(R_chol), z)))
  # Var(u) = tau^2 * Q0^{-1} = tau^2 (D - alpha*W)^{-1}, as required

  if (!is.null(sigma_scale)) u <- u * sigma_scale  # optional extra MEM-driven heterogeneity, layered on top
  u
}

simulate_simple_areal_data <- function(n_s, n_t, W1_list, mem_vals,
                                        beta_s = 0.3,
                                        stationary_sigma = TRUE, sigma0 = 1, sigma_a = 1,
                                        temporal_type = c("none", "smooth", "block"),
                                        baseline_ar1 = 0.7,
                                        ar_shift = 0, ar_scale = 3,
                                        n_seg_time = 5, block_min_mult = 0.2, block_max_mult = 5) {

  temporal_type <- match.arg(temporal_type)

  # spatially (non)stationary marginal noise sd, driven by ground-truth MEM
  sigma_i <- if (stationary_sigma) rep(sigma0, n_s) else sigma0 * exp(sigma_a * mem_vals)

  # temporal behavior: AR coefficient path (smooth) vs. constant AR + variance multiplier (block)
  if (temporal_type == "smooth") {
    raw <- changing_AR1_coefs(n_t, shift = ar_shift, scale = ar_scale)
    beta_t_t   <- 0.6 + 0.3 * raw          # rescaled into (0.3, 0.9), non-explosive
    time_mult  <- rep(1, n_t)              # no separate variance blocks in this mode
  } else if (temporal_type == "block") {
    beta_t_t   <- rep(baseline_ar1, n_t)   # AR persistence stays constant
    time_mult  <- block_temporal_variance(n_t, n_seg_time, block_min_mult, block_max_mult)
  } else {
    beta_t_t   <- rep(baseline_ar1, n_t)
    time_mult  <- rep(1, n_t)
  }

  # fixed spatial coupling (SAR), identical across all settings
  edge_weights <- lapply(W1_list, function(nb) if (length(nb) == 0) numeric(0) else rep(1 / length(nb), length(nb)))
  row_i <- integer(0); col_j <- integer(0); vals <- numeric(0)
  for (i in seq_len(n_s)) {
    nb <- W1_list[[i]]
    if (length(nb) == 0) next
    w_i <- edge_weights[[i]]
    row_i <- c(row_i, rep(i, length(nb))); col_j <- c(col_j, nb); vals <- c(vals, -beta_s * w_i)
  }
  Bmat <- Matrix::sparseMatrix(row_i, col_j, x = vals, dims = c(n_s, n_s))
  SAR_system <- Matrix::Diagonal(n_s) + Bmat  # time-invariant since beta_s is fixed

  z_mat <- matrix(0, n_s, n_t)
  for (t in seq_len(n_t)) {
    # time_mult is a variance multiplier -> take sqrt for the sd scale
    sigma_it <- sigma_i * sqrt(time_mult[t])
    eps_t <- rnorm(n_s, mean = 0, sd = sigma_it)
    rhs <- if (t == 1) eps_t else beta_t_t[t] * z_mat[, t - 1] + eps_t
    z_mat[, t] <- as.numeric(Matrix::solve(SAR_system, rhs))
  }
  as.numeric(z_mat)
}

compute_marginal_variance <- function(n_s, edge_list_spat, beta_s, sigma, edge_weights = NULL) {
  if (is.null(edge_weights)) {
    edge_weights <- lapply(edge_list_spat, function(nb) {
      if (length(nb) == 0L) numeric(0L) else rep(1 / length(nb), length(nb))
    })
  }
  row_i <- integer(0); col_j <- integer(0); vals <- numeric(0)
  for (i in seq_len(n_s)) {
    nb <- edge_list_spat[[i]]
    if (length(nb) == 0L) next
    w_i <- edge_weights[[i]]
    row_i <- c(row_i, rep(i, length(nb)))
    col_j <- c(col_j, nb)
    vals  <- c(vals, -beta_s * w_i)
  }
  B <- Matrix::sparseMatrix(row_i, col_j, x = vals, dims = c(n_s, n_s))
  IminusB <- Matrix::Diagonal(n_s) + B
  Sigma_eps <- Matrix::Diagonal(n_s, x = sigma^2)
  # Sigma_z = (I-B)^{-1} Sigma_eps (I-B)^{-T}, solved column-by-column
  IminusB_inv <- Matrix::solve(IminusB, Matrix::Diagonal(n_s))
  Sigma_z <- IminusB_inv %*% Sigma_eps %*% Matrix::t(IminusB_inv)
  diag(as.matrix(Sigma_z))
}

