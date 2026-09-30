library(ggplot2)
library(patchwork)
library(httpgd)
library(xtable)
hgd()
# Simulation results

sim_plot <- function(data_names, xlab, ylab, title, y_cut = 0.6, nonstat_time = FALSE, nonstat_space = FALSE, data_type = "uniform") {
  res_all <- data.frame(matrix(NA, ncol = 9, nrow = 0))
  names(res_all) <- c("model", "MCC", "n_layers", "seed")
  for (j in 1:3) {
    for (i in seq(1, 200)) {
      for (data_name in data_names) {
            filename <- paste0("/Users/mikasipila/Documents/areal_data/simulations/results/", 
                data_name, "_", j, "_elu_", nonstat_time, "_", nonstat_space, "_", data_type, "_", i, ".RData")
            tryCatch(
                {
                    load(filename)
                    res_all <- rbind(res_all, res_test)
                },
                error = function(e) { },
                warning = function(e) { }
            )
      }
    }
  }
  print(dim(res_all))
  names(res_all)
  res_all <- na.omit(res_all)

  colnames(res_all) <- c("model", "MCC", "n_layers", "seed")
  print(head(res_all))

  summary_df <- do.call(rbind, lapply(split(res_all, list(res_all$model, res_all$n_layers), drop = TRUE), function(df) {
    x <- as.numeric(df$MCC)
    data.frame(
      model = unique(as.character(df$model)),
      n_layers = unique(as.numeric(df$n_layers)),
      mean_MCC = mean(x, na.rm = TRUE),
      sd_MCC = if (length(x) > 1) sd(x, na.rm = TRUE) else 0,
      stringsAsFactors = FALSE
    )
  }))
  summary_df <- summary_df[order(summary_df$model, summary_df$n_layers), ]
  summary_df$mean_sd <- sprintf("%.3f (%.3f)", summary_df$mean_MCC, summary_df$sd_MCC)

  tab_wide <- reshape(
    summary_df[, c("model", "n_layers", "mean_sd")],
    idvar = "model",
    timevar = "n_layers",
    direction = "wide"
  )
  names(tab_wide) <- gsub("^mean_sd\\.", "n_layers_", names(tab_wide))
  tab_wide <- tab_wide[order(tab_wide$model), ]
  rownames(tab_wide) <- tab_wide$model
  tab_wide$model <- NULL

  tab_xt <- xtable(
    tab_wide,
    caption = if (title == "") "Mean MCC across replications (SD in brackets)" else paste0(title, ": mean MCC across replications (SD in brackets)"),
    label = paste0("tab:", gsub("[^A-Za-z0-9_]", "_", data_name)),
    digits = 3
  )

  p <- ggplot(res_all, aes(x = model, y = as.numeric(MCC), fill = as.factor(n_layers), group = interaction(as.factor(n_layers), as.factor(model)))) +
    geom_boxplot() +
    ylim(c(0, 1)) +
    coord_cartesian(ylim = c(y_cut, 1)) +
    labs(x = xlab, y = ylab, fill = "Number of mixing layers", title = title) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      strip.text.y = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 12)
    ) +
    theme(legend.position = "bottom")
  return(list(plot = p, table = tab_xt))
}

#p1 <- sim_plot("res_ar1_blocks", "", "", "Setting 1")
p1_0 <- sim_plot(c("res_areal_sar", "res_areal_sar_fica"), "", "MCC", "Setting A1", nonstat_time = FALSE, nonstat_space = FALSE)
p2_0 <- sim_plot(c("res_areal_sar", "res_areal_sar_fica"), "", "MCC", "Setting A2", nonstat_time = FALSE, nonstat_space = TRUE)
p3_0 <- sim_plot(c("res_areal_sar", "res_areal_sar_fica"), "", "MCC", "Setting A3", nonstat_time = TRUE, nonstat_space = TRUE)
p1_1 <- sim_plot(c("res_areal_sar", "res_areal_sar_fica"), "", "MCC", "Setting B1", nonstat_time = FALSE, nonstat_space = FALSE, data_type = "clustered")
p2_1 <- sim_plot(c("res_areal_sar", "res_areal_sar_fica"), "", "MCC", "Setting B2", nonstat_time = FALSE, nonstat_space = TRUE, data_type = "clustered")
p3_1 <- sim_plot(c("res_areal_sar", "res_areal_sar_fica"), "", "MCC", "Setting B3", nonstat_time = TRUE, nonstat_space = TRUE, data_type = "clustered")


# Combine plots into a 2x3 layout with a shared legend at the bottom using patchwork
## Combine plots into a 2x3 layout with a shared legend at the bottom using patchwork
# use plot_annotation to set a shared theme to avoid '&' method issues
p_combined <- (p1_0$plot + p2_0$plot + p3_0$plot + p1_1$plot + p2_1$plot + p3_1$plot) +
  plot_layout(ncol = 3, nrow = 2, guides = "collect") & theme(legend.position = "bottom")
p_combined

# Print tables if needed
p1_0$table
p1_1$table

p2_0$table
p2_1$table

p3_0$table
p3_1$table

