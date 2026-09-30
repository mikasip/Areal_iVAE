
#' Compute Moran Eigenvector Maps from an Adjacency List
#'
#' @description
#' Constructs Moran Eigenvector Maps (MEM) from the spatial weight matrix
#' implied by \code{adj_list}. The weight matrix is taken as symmetric binary
#' (1 if neighbors, 0 otherwise) — do NOT row-standardise, as that would
#' break the symmetry required for MEM.
#'
#' @param adj_list     List of length N; each element is a 1-based integer
#'   vector of neighbor indices. Must represent an undirected graph
#'   (if i is a neighbor of j, j must be a neighbor of i).
#' @param K            Integer. Number of MEM vectors to return. Default 10.
#' @param positive_only Logical. If TRUE (default), return only eigenvectors
#'   with positive eigenvalues, i.e. those encoding positive spatial
#'   autocorrelation — the natural choice for geographic areal data where
#'   nearby areas are expected to be similar. If FALSE, return the K vectors
#'   with the largest absolute eigenvalues (includes negative autocorrelation).
#'
#' @return A list with:
#'   \item{mem}{N x K matrix of MEM vectors (columns).}
#'   \item{values}{K eigenvalues, in decreasing order. Each eigenvalue
#'     is proportional to the Moran's I of the corresponding MEM vector.}
#'   \item{W_tilde}{The doubly-centred weight matrix H*W*H (N x N dense).}
#'   \item{moran_i}{Moran's I values for each MEM vector
#'     (= eigenvalue / sum(W)).}
#'
#' @details
#' \strong{Memory note:} The centred matrix W_tilde is dense (N x N). For
#' N > ~3000 this may use substantial memory. For large N consider using
#' only a moderate K (e.g. 10-20) and setting \code{positive_only = TRUE}.
#'
#' \strong{Why not row-standardise W?} Row-standardisation (W[i,j] / d_i)
#' produces an asymmetric matrix, which breaks the eigenvector decomposition
#' that gives MEM their Moran's I interpretation. Always use the raw
#' symmetric binary (or boundary-length weighted) adjacency for MEM.
#'
#' @examples
#' \dontrun{
#' adj_list   <- make_adjacency_delaunay(locations)
#' mem_result <- compute_MEM_vectors(adj_list, K = 10)
#' aux_data   <- mem_result$mem
#' }
#' @export
compute_MEM_vectors <- function(adj_list, K = 10L, positive_only = TRUE) {
    n <- length(adj_list)
 
    # ---- Build symmetric binary adjacency matrix ----------------------------
    # Do NOT row-standardise: MEM requires symmetric W
    row_i <- integer(0); col_j <- integer(0)
    for (i in seq_len(n)) {
        for (j in adj_list[[i]]) {
            row_i <- c(row_i, i)
            col_j <- c(col_j, j)
        }
    }
    W_sparse <- Matrix::sparseMatrix(row_i, col_j, x = 1.0, dims = c(n, n))
    # Symmetrize in case adj_list has any one-directional edges
    W_sparse <- (W_sparse + Matrix::t(W_sparse)) / 2.0
    W        <- as.matrix(W_sparse)    # dense n x n (needed for centering)
 
    # ---- Doubly-centre W: W_tilde = H W H ----------------------------------
    # H = I - 11^T/n  is the centring matrix.
    # Entry (i,j) of H*W*H:
    #   = W[i,j] - rowmean[i] - rowmean[j] + grandmean
    # For symmetric W, rowmeans = colmeans, so this is symmetric.
    row_means  <- rowMeans(W)            # = colMeans(W) for symmetric W
    grand_mean <- mean(W)                # = sum(W) / n^2
    # outer(a, b, "+") gives matrix with entry (i,j) = a[i] + b[j]
    W_tilde    <- W - outer(row_means, row_means, FUN = "+") + grand_mean
 
    # ---- Eigen decomposition ------------------------------------------------
    # W_tilde is real symmetric; eigen() returns eigenvalues in DECREASING order
    eig <- eigen(W_tilde, symmetric = TRUE)
 
    # ---- Select eigenvectors ------------------------------------------------
    if (positive_only) {
        # Keep eigenvectors with positive eigenvalues:
        # positive lambda <=> positive Moran's I <=> positive autocorrelation
        pos_idx <- which(eig$values > 1e-8)
        n_pos   <- length(pos_idx)
 
        if (n_pos == 0L) stop("No MEM vectors with positive eigenvalues found.")
        if (n_pos < K) {
            warning(sprintf(
                paste0("Only %d MEM vectors with positive eigenvalues available ",
                       "(K = %d requested). Returning all %d."),
                n_pos, K, n_pos
            ))
            K <- n_pos
        }
        # Already in decreasing eigenvalue order (most positive autocorrelation first)
        sel_idx <- pos_idx[seq_len(K)]
 
    } else {
        # Return K eigenvectors with largest |eigenvalue|
        # (most extreme autocorrelation, positive or negative)
        ord     <- order(abs(eig$values), decreasing = TRUE)
        sel_idx <- ord[seq_len(min(K, n - 1L))]
    }
 
    mem_vecs <- eig$vectors[, sel_idx, drop = FALSE]  # N x K
    mem_vals <- eig$values[sel_idx]
 
    # ---- Sign convention: largest-magnitude entry positive ------------------
    for (k in seq_len(ncol(mem_vecs))) {
        if (mem_vecs[which.max(abs(mem_vecs[, k])), k] < 0)
            mem_vecs[, k] <- -mem_vecs[, k]
    }
 
    colnames(mem_vecs) <- paste0("MEM_", seq_len(ncol(mem_vecs)))
 
    # Moran's I = lambda / sum(W)  (a standard result)
    sum_W   <- sum(W)
    moran_i <- mem_vals / sum_W
 
    list(
        mem     = mem_vecs,
        values  = mem_vals,
        moran_i = moran_i,
        W_tilde = W_tilde
    )
}

