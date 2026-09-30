library(deldir)
library(sf)
library(igraph)

# ---------------------------------------------------------------
# 1. Fine-grained point layer (much denser than n_s) + Delaunay graph
# ---------------------------------------------------------------
generate_fine_layer <- function(n_fine, domain = c(0, 1, 0, 1)) {
  pts <- matrix(runif(2 * n_fine, min = domain[1], max = domain[2]), ncol = 2)
  del <- deldir::deldir(pts[, 1], pts[, 2], rw = domain)
  edges <- as.matrix(del$delsgs[, c("ind1", "ind2")])
  g <- igraph::simplify(igraph::graph_from_edgelist(edges, directed = FALSE))
  list(pts = pts, del = del, graph = g)
}

# ---------------------------------------------------------------
# 2. Seed placement controlling area-size heterogeneity.
#    Reuses your existing generate_points(): "uniform" seeds -> fairly
#    even area sizes; "clustered"/"normal" seeds -> large areas in
#    sparse regions, small areas in dense (seed-rich) regions.
# ---------------------------------------------------------------
place_seeds <- function(fine_pts, n_s, seed_type = c("uniform", "clustered")) {
  seed_type <- match.arg(seed_type)
  seed_coords <- generate_points(n_s, type = seed_type)   # your existing function
  # snap each seed to its nearest fine-layer point (so seeds are graph nodes)
  nn <- FNN::get.knnx(fine_pts, seed_coords, k = 1)$nn.index[, 1]
  unique(nn)
}

place_seeds_fixed <- function(fine_pts, n_s, seed_type = c("uniform", "clustered"), n_clusters = 10, cluster_sd = 0.04, xy_coef = 1) {
  seed_type  <- match.arg(seed_type)
  seed_coords <- generate_points(n_s, type = seed_type, n_clusters, cluster_sd, xy_coef)   # your existing function
  n_fine <- nrow(fine_pts)

  k <- min(20, n_fine)   # candidate pool size per seed; grown adaptively if needed
  nn <- FNN::get.knnx(fine_pts, seed_coords, k = k)$nn.index   # n_s x k

  assigned <- integer(n_s)
  used <- logical(n_fine)

  for (i in seq_len(n_s)) {
    cands <- nn[i, ][!used[nn[i, ]]]
    k_i <- k
    while (length(cands) == 0) {
      k_i <- k_i * 2
      if (k_i > n_fine) {
        stop("Fine layer too small relative to n_s; increase n_fine.")
      }
      nn_i  <- FNN::get.knnx(fine_pts, seed_coords[i, , drop = FALSE], k = k_i)$nn.index[1, ]
      cands <- nn_i[!used[nn_i]]
    }
    assigned[i] <- cands[1]
    used[cands[1]] <- TRUE
  }
  assigned
}

# ---------------------------------------------------------------
# 3. Graph-Voronoi assignment: each fine point -> nearest seed by
#    graph (hop) distance. This is the "principled size control" step.
# ---------------------------------------------------------------
assign_graph_voronoi <- function(g, seed_idx) {
  dist_mat <- igraph::distances(g, v = seed_idx, to = igraph::V(g), weights = NA)
  apply(dist_mat, 2, which.min)   # length = n_fine, values in 1..length(seed_idx)
}

# ---------------------------------------------------------------
# 4. Build the areal adjacency matrix COMBINATORIALLY from the fine
#    Delaunay graph (two areal units are neighbors iff some pair of
#    fine points in each cluster are Delaunay-connected). This avoids
#    the invalid-geometry issues you hit with poly2nb, since it never
#    depends on detecting shared polygon boundaries geometrically.
# ---------------------------------------------------------------
build_areal_adjacency <- function(g, cluster_id) {
  n_areas <- length(unique(cluster_id))
  edge_df <- igraph::as_data_frame(g, what = "edges")
  c1 <- cluster_id[as.integer(edge_df$from)]
  c2 <- cluster_id[as.integer(edge_df$to)]
  cross <- c1 != c2
  pairs <- unique(data.frame(a = pmin(c1[cross], c2[cross]),
                              b = pmax(c1[cross], c2[cross])))
  W <- matrix(0, n_areas, n_areas)
  W[cbind(pairs$a, pairs$b)] <- 1
  W[cbind(pairs$b, pairs$a)] <- 1
  W
}

# ---------------------------------------------------------------
# 5. (Optional, for plotting only) merge fine Voronoi cells into
#    final areal polygons via st_union, grouped by cluster.
# ---------------------------------------------------------------
build_areal_polygons <- function(fine, cluster_id) {
  tiles <- deldir::tile.list(fine$del)
  poly_list <- lapply(seq_along(tiles), function(i) {
    tl <- tiles[[i]]
    coords <- cbind(tl$x, tl$y)
    coords <- rbind(coords, coords[1, , drop = FALSE])
    sf::st_polygon(list(coords))
  })

  voronoi_sf <- sf::st_sf(
    cluster = cluster_id,
    geometry = sf::st_sfc(poly_list)
  )

  areas_list <- lapply(split(voronoi_sf, voronoi_sf$cluster), function(g) {
    sf::st_sf(
      cluster = unique(g$cluster),
      geometry = sf::st_sfc(sf::st_union(sf::st_geometry(g)))
    )
  })

  areas_sf <- do.call(rbind, areas_list)
  sf::st_make_valid(areas_sf)
}
