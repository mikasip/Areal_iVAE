
# ---------------------------------------------------------------
# Full pipeline
# ---------------------------------------------------------------
n_s    <- 50
n_fine <- round(3000/4)   # denser fine layer -> smoother, more organic area boundaries

fine <- generate_fine_layer(n_fine)
seed_idx <- place_seeds_fixed(fine$pts, n_s, seed_type = "uniform", 10, 0.04, 1)   # or "uniform"
cluster_id <- assign_graph_voronoi(fine$graph, seed_idx)

# renumber clusters 1..n_s_actual (some seeds can end up unreachable/duplicated)
cluster_id <- as.integer(factor(cluster_id))
n_s_actual <- length(unique(cluster_id))

W1 <- build_areal_adjacency(fine$graph, cluster_id)   # drop-in replacement for old W1

# check size heterogeneity achieved
area_sizes <- table(cluster_id)
cat(sprintf("n_areas=%d  size range=[%d, %d]  size CV=%.2f\n",
            n_s_actual, min(area_sizes), max(area_sizes),
            sd(area_sizes) / mean(area_sizes)))

# optional: build polygons only if you want to plot the areal units
areas_sf <- build_areal_polygons(fine, cluster_id)
plot(sf::st_geometry(areas_sf))
hist(log(area_sizes))

hist(rowSums(W1))
W2 <- order2_neighbors(W1, cumulative = TRUE)
hist(rowSums(W2))


library(sf)
library(ggplot2)
library(patchwork)

# ---- Find the area closest to the center of the domain ----
centroids <- sf::st_centroid(areas_sf)
coords_cent <- sf::st_coordinates(centroids)
dist_to_center <- sqrt((coords_cent[,1] - 0.5)^2 + (coords_cent[,2] - 0.5)^2)
central_idx <- which.min(dist_to_center)

cat("Target central area index:", central_idx, "\n")

# ---- Compute order-1 and order-2 (exclusive) neighbor sets ----
order1_nb <- which(W1[central_idx, ] == 1)

W2_exclusive <- order2_neighbors(W1, cumulative = FALSE)  # ONLY 2-hop, not 1-hop
order2_nb <- which(W2_exclusive[central_idx, ] == 1)

# sanity check: no overlap between the two sets, and neither contains the area itself
stopifnot(length(intersect(order1_nb, order2_nb)) == 0)
stopifnot(!(central_idx %in% order1_nb), !(central_idx %in% order2_nb))

# ---- Build a categorical label for each area, for each of the two plots ----
areas_sf$order1_group <- "Other"
areas_sf$order1_group[order1_nb] <- "1st-order neighbor"
areas_sf$order1_group[central_idx] <- "Target area"
areas_sf$order1_group <- factor(areas_sf$order1_group,
                                 levels = c("Target area", "1st-order neighbor", "Other"))

areas_sf$order2_group <- "Other"
areas_sf$order2_group[order2_nb] <- "2nd-order neighbor"
areas_sf$order2_group[order1_nb] <- "1st-order neighbor"
areas_sf$order2_group[central_idx] <- "Target area"
areas_sf$order2_group <- factor(areas_sf$order2_group,
                                 levels = c("Target area", "1st-order neighbor",
                                            "2nd-order neighbor", "Other"))

# ---- Consistent color scheme across both plots ----
group_colors <- c(
  "Target area"       = "red",
  "1st-order neighbor"  = "gray50",
  "2nd-order neighbor"  = "gray80",
  "Other"               = "white"
)

p1 <- ggplot(areas_sf) +
  geom_sf(aes(fill = order1_group), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = group_colors, drop = FALSE, name = NULL) +
  labs(title = "1st-order neighborhood") +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom")

p2 <- ggplot(areas_sf) +
  geom_sf(aes(fill = order2_group), color = "black", linewidth = 0.2) +
  scale_fill_manual(values = group_colors, drop = FALSE, name = NULL) +
  theme_minimal() +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.position = "bottom")

p1 + p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p2
