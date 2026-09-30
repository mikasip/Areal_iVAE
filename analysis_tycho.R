library(NonlinearBSS)
library(ggplot2)
library(xtable)
library(patchwork)
library(gridExtra)
library(fastICA)
library(ggplot2)
library(scales)
library(sf)
library(maps)
library(kernelshap)
library(dplyr)
library(patchwork)
library(tools)

source("helpers/compute_moran_eigenvectors.R")
source("helpers/helpers.R")

load("tycho_aggregated.RData")

unified_df_check <- unified_df
unified_df_check$time <- NULL
row_info <- data.frame(
  row_name = rownames(unified_df_check),
  state = sapply(strsplit(rownames(unified_df_check), "_"), `[`, 1),
  time = sapply(strsplit(rownames(unified_df_check), "_"), `[`, 2),
  stringsAsFactors = FALSE
)

n_t <- length(unique(row_info$time))
n_a <- length(unique(row_info$state))

# Get US states shapefile
us_states <- st_as_sf(map('state', plot = FALSE, fill = TRUE))
us_states$state <- toupper(sub("^([a-z]+):(.+)$", "\\2", us_states$ID))

# Check that state names in data match state names from shapefile
data_states <- unique(row_info$state)
shapefile_states <- unique(us_states$state)
unmatched_data <- setdiff(data_states, shapefile_states)
unmatched_shapefile <- setdiff(shapefile_states, data_states)

if (length(unmatched_data) > 0) {
  cat("WARNING: States in data not found in shapefile:\n")
  print(unmatched_data)
}
if (length(unmatched_shapefile) > 0) {
  cat("WARNING: States in shapefile not found in data:\n")
  print(unmatched_shapefile)
}
if (length(unmatched_data) == 0 && length(unmatched_shapefile) == 0) {
  cat("All state names match between data and shapefile.\n")
}

# Create a data frame with state-level counts at first time point
state_data <- data.frame()
for (state in unique(row_info$state)) {
  state_time_rows <- row_info$row_name[row_info$state == state & 
                                        grepl(format(first_time, "%Y-%m-%d"), row_info$row_name)]
  if (length(state_time_rows) > 0) {
    count <- unified_df_check[state_time_rows, measles_condition]
    # replace NA/zero with a small positive value so log scale can show small counts
    count_val <- ifelse(is.na(count) | count <= 0, 0.1, count)
    state_data <- rbind(state_data, 
                       data.frame(state = state, measles_count = count_val))
  }
}

# Merge with shapefile
us_states <- merge(us_states, state_data, by.x = "state", by.y = "state", all.x = TRUE)

log_data <- log(unified_df + 1)
# Proportion of missing values for each variable:
miss_data <- as.data.frame(apply(log_data, 2, function(x) sum(is.na(x)) / length(x)))
xtable(miss_data, digits = 3)

state_names <- unique(row_info$state)   # your actual 51 area names, exact spelling/case

adj_list <- list(
  AL=c("FL","GA","MS","TN"), AK=c(), AZ=c("CA","NV","NM","UT"),
  AR=c("LA","MS","MO","OK","TN","TX"), CA=c("AZ","NV","OR"),
  CO=c("KS","NE","NM","OK","UT","WY"), CT=c("MA","NY","RI"),
  DE=c("MD","NJ","PA"), DC=c("MD","VA"), FL=c("AL","GA"),
  GA=c("AL","FL","NC","SC","TN"), HI=c(),
  ID=c("MT","NV","OR","UT","WA","WY"), IL=c("IN","IA","KY","MO","WI"),
  IN=c("IL","KY","MI","OH"), IA=c("IL","MN","MO","NE","SD","WI"),
  KS=c("CO","MO","NE","OK"), KY=c("IL","IN","MO","OH","TN","VA","WV"),
  LA=c("AR","MS","TX"), ME=c("NH"), MD=c("DE","PA","VA","WV","DC"),
  MA=c("CT","NH","NY","RI","VT"), MI=c("IN","OH","WI"),
  MN=c("IA","ND","SD","WI"), MS=c("AL","AR","LA","TN"),
  MO=c("AR","IL","IA","KS","KY","NE","OK","TN"), MT=c("ID","ND","SD","WY"),
  NE=c("CO","IA","KS","MO","SD","WY"), NV=c("AZ","CA","ID","OR","UT"),
  NH=c("ME","MA","VT"), NJ=c("DE","NY","PA"), NM=c("AZ","CO","OK","TX"),
  NY=c("CT","MA","NJ","PA","VT"), NC=c("GA","SC","TN","VA"),
  ND=c("MN","MT","SD"), OH=c("IN","KY","MI","PA","WV"),
  OK=c("AR","CO","KS","MO","NM","TX"), OR=c("CA","ID","NV","WA"),
  PA=c("DE","MD","NJ","NY","OH","WV"), RI=c("CT","MA"), SC=c("GA","NC"),
  SD=c("IA","MN","MT","NE","ND","WY"),
  TN=c("AL","AR","GA","KY","MS","MO","NC","VA"), TX=c("AR","LA","NM","OK"),
  UT=c("AZ","CO","ID","NV","WY"), VT=c("MA","NH","NY"),
  VA=c("KY","MD","NC","TN","WV","DC"), WA=c("ID","OR"),
  WV=c("KY","MD","OH","PA","VA"), WI=c("IL","IA","MI","MN"),
  WY=c("CO","ID","MT","NE","SD","UT")
)

edges <- do.call(rbind, lapply(names(adj_list), function(a) {
  if (length(adj_list[[a]]) == 0) return(NULL)
  data.frame(from = a, to = adj_list[[a]])
}))

abb2name <- setNames(toupper(state.name), state.abb)
abb2name["DC"] <- "DISTRICT OF COLUMBIA"  # tweak to match your data's spelling

edges$from_name <- abb2name[edges$from]
edges$to_name   <- abb2name[edges$to]

n_s <- length(state_names)
nb_mat <- matrix(0L, n_s, n_s, dimnames = list(state_names, state_names))

for (i in seq_len(nrow(edges))) {
  a <- edges$from_name[i]; b <- edges$to_name[i]
  if (a %in% state_names && b %in% state_names) {
    nb_mat[a, b] <- 1L
    nb_mat[b, a] <- 1L
  }
}

cat("Symmetric:", isSymmetric(nb_mat), "\n")
cat("Zero diagonal:", all(diag(nb_mat) == 0), "\n")
cat("Total edges (undirected):", sum(nb_mat) / 2, "\n")
cat("Dim:", dim(nb_mat), "\n")

areas_in_log_data_order <- unique(sapply(strsplit(rownames(log_data), "_"), `[`, 1))
nb_mat_ordered <- nb_mat[areas_in_log_data_order, areas_in_log_data_order]

identical(rownames(nb_mat_ordered), areas_in_log_data_order)  # should be TRUE
isSymmetric(nb_mat_ordered)  # should be TRUE
all(diag(nb_mat_ordered) == 0)  # should be TRUE, no self-loops

# --- Convert to adjacency list (names) ---
adj_list_neighbors <- setNames(
  lapply(rownames(nb_mat_ordered), function(s) colnames(nb_mat_ordered)[nb_mat_ordered[s, ] == 1]),
  rownames(nb_mat_ordered)
)

# compute_MEM_vectors expects integer indices rather than names:
adj_list_idx <- lapply(adj_list_neighbors, function(neigh) match(neigh, rownames(nb_mat_ordered)))

mem_obj <- compute_MEM_vectors(adj_list_idx, K = 15, positive_only = TRUE)

# Plot the MEM vectors using ggplot2 and geom_sf
mem_df <- data.frame(
  state = rownames(nb_mat_ordered)
)
mem_df <- append(mem_df, as.data.frame(mem_obj$mem))

mem_plots <- lapply(seq_len(ncol(mem_obj$mem)), function(i) {
  col_name <- paste0("MEM_", i)

  us_states_i <- us_states
  us_states_i$mem_val <- mem_df[[col_name]][match(us_states_i$state, mem_df$state)]

  ggplot(us_states_i) +
    geom_sf(aes(fill = mem_val)) +
    scale_fill_viridis_c(option = "plasma") +
    labs(fill = "MEM value", title = paste("MEM Vector", i)) +
    theme_minimal()
})

# Use patchwork to combine the MEMs into 3 columns and 5 rows
combined_mem_plot <- patchwork::wrap_plots(
  mem_plots,
  ncol = 3,
  nrow = 5
)

print(combined_mem_plot)

aux_data_areal <- mem_obj$mem # n_obs x K matrix of MEM vectors (columns)
# Replicate each MEM vector for each time point, so that aux_data_areal is n_obs x K
aux_data_areal <- aux_data_areal[rep(seq_len(nrow(aux_data_areal)), times = n_t), ]

temp_spat_data <- matrix(rnorm(n_a * 2), nrow = n_a, ncol = 2)
temp_spat_data <- temp_spat_data[rep(seq_len(nrow(temp_spat_data)), times = n_t), ]

time_vec <- rep(seq_len(n_t), each = n_a)  # zero-based time points (months since first observation)

seed <- 10092026
n_epochs <- 500

# Start by selecting the number of ICs
iVAEar_fit2 <- iVAEar_radial(as.matrix(log_data), temp_spat_data, 
                            time_vec, 2, n_s = n_a, aux_data = aux_data_areal, 
                            spatial_basis = c(0), ar_order = 2, epochs = n_epochs, seed = seed,
                            batch_size = 64, error_dist = "gaussian", get_elbo = TRUE)
iVAEar_fit3 <- iVAEar_radial(as.matrix(log_data), temp_spat_data, 
                            time_vec, 3, n_s = n_a, aux_data = aux_data_areal, 
                            spatial_basis = c(0), ar_order = 2, epochs = n_epochs, seed = seed,
                            batch_size = 64, error_dist = "gaussian", get_elbo = TRUE)
iVAEar_fit4 <- iVAEar_radial(as.matrix(log_data), temp_spat_data, 
                            time_vec, 4, n_s = n_a, aux_data = aux_data_areal, 
                            spatial_basis = c(0), ar_order = 2, epochs = n_epochs, seed = seed,
                            batch_size = 64, error_dist = "gaussian", get_elbo = TRUE)
iVAEar_fit5 <- iVAEar_radial(as.matrix(log_data), temp_spat_data, 
                            time_vec, 5, n_s = n_a, aux_data = aux_data_areal, 
                            spatial_basis = c(0), ar_order = 2, epochs = n_epochs, seed = seed,
                            batch_size = 64, error_dist = "gaussian", get_elbo = TRUE)
iVAEar_fit6 <- iVAEar_radial(as.matrix(log_data), temp_spat_data, 
                            time_vec, 6, n_s = n_a, aux_data = aux_data_areal,
                            spatial_basis = c(0), ar_order = 2, epochs = n_epochs, seed = seed,
                            batch_size = 64, error_dist = "gaussian", get_elbo = TRUE)

elbos <- data.frame(
  n_ic = 2:6,
  elbo = c(as.numeric(iVAEar_fit2$elbo), as.numeric(iVAEar_fit3$elbo), as.numeric(iVAEar_fit4$elbo), as.numeric(iVAEar_fit5$elbo), as.numeric(iVAEar_fit6$elbo))
)

ggplot(elbos, aes(x = n_ic, y = elbo)) +
  geom_line() +
  geom_point() +
  labs(x = "Latent Dimension", y = "ELBO") +
  theme_minimal()

head(iVAEar_fit$elbo)

pAIC <- function(elbo, p) {
  -2 * elbo + 2 * p
}
which.min(apply(elbos, 1, function(row) pAIC(row["elbo"], row["n_ic"])))

# Final model:
p <- 6
cor(iVAEar_fit6$IC)

# Interpretations of the ICs using kernelshap
library(kernelshap)

X <- as.data.frame(iVAEar_fit6$IC)
n <- nrow(X)
set.seed(10092026)
bg_X <- X[1:200, ]
explainer2 <- kernelshap(iVAEar_fit6, X, bg_X = bg_X, pred_fun = function(object, X) {
  pred <- predict(object, newdata = as.matrix(X), IC_to_data = TRUE)
  return(pred)
})

explainer2$baseline
shap_vars2 <- data.frame(matrix(NA, ncol = p, nrow = ncol(log_data)))
rownames(shap_vars2) <- colnames(log_data)
colnames(shap_vars2) <- sapply(1:p, function(i) paste0("IC", i))
i <- 1
for (l in explainer2$S) {
  shap_vars2[i, ] <- apply(l, 2, function(x) mean(abs(x)))
  i <- i + 1
}
shap_vars_scaled2 <- sweep(shap_vars2, 1, rowSums(shap_vars2), "/")
shap_vars_scaled2
colMeans(shap_vars_scaled2)
col_ord <- order(-colMeans(shap_vars_scaled2))
shap_vars_scaled2 <- shap_vars_scaled2[, col_ord]

# Convert to xtable
df_fmt <- as.data.frame(
  lapply(shap_vars_scaled2, function(col) { 
    if (is.numeric(col)) {
      sapply(col, function(x) {
        if (x > 0.2) {
          sprintf("\\textbf{%.3f}", x)
        } else {
          sprintf("%.3f", x)
        }
      })
    } else {
      col
    }
  })
)

rownames(df_fmt) <- rownames(shap_vars_scaled2)
tab <- xtable(df_fmt, digits = 3)
print(tab, sanitize.text.function = identity, include.rownames = TRUE)
xtable(data.frame(t(colMeans(shap_vars_scaled2))), digits = 3)


# --- 1. Attach IC values to area/time metadata ---
ic_matrix <- iVAEar_fit6$IC
rownames(ic_matrix) <- rownames(log_data)

row_meta <- data.frame(
  row_name = rownames(log_data),
  area = sapply(strsplit(rownames(log_data), "_"), `[`, 1),
  time = as.Date(sapply(strsplit(rownames(log_data), "_"), function(x) x[length(x)])),
  stringsAsFactors = FALSE
)

ic_df <- cbind(row_meta, as.data.frame(ic_matrix))
colnames(ic_df)[-(1:3)] <- paste0("IC", seq_len(ncol(ic_matrix)))

head(ic_df)

us_states_valid <- st_make_valid(us_states)
us_states_valid <- us_states_valid[!(us_states_valid$state %in% c("ALASKA","HAWAII")), ]

make_ic_plots <- function(ic_idx, title_idx, area_to_plot, ic_df, us_states_valid, add_title = FALSE, prefix = "IC") {
  ic_col <- paste0(prefix, ic_idx)

  df <- ic_df
  df$month <- as.integer(format(df$time, "%m"))
  df$season <- dplyr::case_when(
    df$month %in% c(6, 7, 8)  ~ "Summer",
    df$month %in% c(12, 1, 2) ~ "Winter",
    TRUE                      ~ NA_character_
  )

  seasonal_avg <- df %>%
    filter(!is.na(season)) %>%
    group_by(area, season) %>%
    summarise(ic_mean = mean(.data[[ic_col]], na.rm = TRUE), .groups = "drop")

  shared_limits <- range(seasonal_avg$ic_mean, na.rm = TRUE)

  # the highlighted area's outline geometry, reused for both season maps
  highlight_geom <- us_states_valid %>% filter(state == area_to_plot)

  make_season_map <- function(season_name, add_title = FALSE) {
    df_season <- seasonal_avg %>% filter(season == season_name)
    us_states_season <- us_states_valid %>%
      left_join(df_season, by = c("state" = "area"))

    gplot <- ggplot(us_states_season) +
      geom_sf(aes(fill = ic_mean), color = "white", linewidth = 0.1) +
      geom_sf(data = highlight_geom, fill = NA, color = "red", linewidth = 1) +
      scale_fill_viridis_c(na.value = "lightgray", limits = shared_limits, name = paste0(prefix, title_idx)) +
      theme_minimal() +
      theme(axis.text = element_blank(), axis.ticks = element_blank(),
            plot.title = element_text(size = 10))
    if (add_title) gplot <- gplot + labs(title = season_name)
    return(gplot)
  }

  summer_plot <- make_season_map("Summer", add_title = add_title)
  winter_plot <- make_season_map("Winter", add_title = add_title)

  ic_area <- df[df$area == area_to_plot, c("time", ic_col)]
  temporal_plot <- ggplot(ic_area, aes(x = time, y = .data[[ic_col]])) +
    geom_line() +
    geom_point(size = 0.6) +
    ylab(paste0(prefix, title_idx))
    theme_minimal() +
    theme(plot.title = element_text(size = 10))
  if (add_title) temporal_plot <- temporal_plot + labs(title = toTitleCase(tolower(area_to_plot)))

  list(summer = summer_plot, winter = winter_plot, temporal = temporal_plot)
}

area_to_plot <- "MISSOURI"   # fixed across ICs for easy visual comparison; swap in your chosen area

all_plots <- list()

for (i in 1:6) {
  ic_idx <- col_ord[i]
  p <- make_ic_plots(ic_idx, i, area_to_plot, ic_df, us_states_valid, (i == 1))
  all_plots <- c(all_plots, list(p$summer, p$winter, p$temporal))
}

full_grid <- wrap_plots(all_plots, ncol = 3, nrow = 6, widths = c(0.35, 0.35, 0.3))
full_grid


# --------------------- fastICA --------------------------#
p <- 6
log_data_complete <- na.omit(log_data)
nrow(log_data)
nrow(log_data_complete)
res_fica <- fastICA(as.matrix(log_data_complete), p)

head(res_fica$S %*% res_fica$A)
loading_df <- as.matrix(t(res_fica$A))
colnames(loading_df) <- sapply(1:p, function(i) paste0("IC", i))
rownames(loading_df) <- colnames(log_data)

loading_df_norm <- sweep(loading_df, 1, apply(loading_df, 1, function(x) sum(abs(x))), "/")

tab <- xtable(loading_df_norm, digits = 3)
print(tab, sanitize.text.function = identity, include.rownames = TRUE)

area_to_plot <- "MISSOURI"   # fixed across ICs for easy visual comparison; swap in your chosen area
all_plots <- list()

row_meta_complete <- row_meta[complete.cases(log_data), ]
ic_df_fica <- cbind(row_meta_complete, as.data.frame(res_fica$S))
colnames(ic_df_fica)[-(1:3)] <- paste0("IC", seq_len(ncol(log_data)))

for (i in 1:6) {
  ic_idx <- i
  p <- make_ic_plots(ic_idx, i, area_to_plot, ic_df_fica, us_states_valid, (i == 1))
  all_plots <- c(all_plots, list(p$summer, p$winter, p$temporal))
}

full_grid <- wrap_plots(all_plots, ncol = 3, nrow = 6, widths = c(0.35, 0.35, 0.3))
full_grid


# --------------------- PCA --------------------------#
log_data_complete <- na.omit(log_data)
nrow(log_data)
nrow(log_data_complete)
res_pca <- prcomp(as.matrix(log_data_complete))

loading_df <- as.matrix(res_pca$rotation)

tab <- xtable(loading_df, digits = 3)
print(tab, sanitize.text.function = identity, include.rownames = TRUE)

area_to_plot <- "MISSOURI"   # fixed across ICs for easy visual comparison; swap in your chosen area
all_plots <- list()

row_meta_complete <- row_meta[complete.cases(log_data), ]
ic_df_pca <- cbind(row_meta_complete, as.data.frame(-res_pca$x))
colnames(ic_df_pca)[-(1:3)] <- paste0("PC", seq_len(ncol(log_data)))

for (i in 1:6) {
  ic_idx <- i
  p <- make_ic_plots(ic_idx, i, area_to_plot, ic_df_pca, us_states_valid, (i == 1), prefix = "PC")
  all_plots <- c(all_plots, list(p$summer, p$winter, p$temporal))
}

full_grid <- wrap_plots(all_plots, ncol = 3, nrow = 6, widths = c(0.35, 0.35, 0.3))
full_grid
