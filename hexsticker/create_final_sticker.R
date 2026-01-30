# Create final CU Boulder hex sticker with black background
library(ggplot2)
library(hexSticker)

source("saddle_function.R")

# CU Boulder brand colors
cu_gold <- "#CFB87C"
cu_gold_accessible <- "#8D7334"
cu_black <- "#000000"

# Generate 2 optimization paths: left and right of saddle
cat("Generating 2 optimization paths from left and right of saddle...\n")

starting_points <- list(
  c(-0.05, 0.1),   # Left of saddle
  c(0.05, 0.1)     # Right of saddle
)

# Collect all paths
all_paths <- list()

devtools::load_all()

for (i in seq_along(starting_points)) {
  x0 <- starting_points[[i]]
  cat("Path", i, "from (", x0[1], ",", x0[2], "): ")

  result <- arcopt(
    x0 = x0,
    fn = saddle_fn,
    gr = saddle_gr,
    hess = saddle_hess,
    control = list(trace = 3)
  )

  # Combine initial point with trajectory
  path_matrix <- result$trace$x
  path <- rbind(x0, path_matrix)

  all_paths[[i]] <- data.frame(
    x = path[, 1],
    y = path[, 2],
    path_id = i
  )

  cat(nrow(path), "points, converged to (",
      round(result$par[1], 2), ",", round(result$par[2], 2), ")\n")
}

# Combine all paths into one dataframe
paths_df <- do.call(rbind, all_paths)

cat("\nTotal paths:", length(all_paths), "\n")
cat("Total points:", nrow(paths_df), "\n")

# Create grid for contour plot
grid <- expand.grid(
  x = seq(-1.8, 1.8, length.out = 300),
  y = seq(-1.2, 1.2, length.out = 300)
)
grid$z <- saddle_fn_vectorized(grid$x, grid$y)

# Saddle point and minimum values
min_val <- min(grid$z)  # -0.25 at minima
saddle_val <- 0          # at saddle

cat("\nMinimum function value:", min_val, "\n")
cat("Saddle point function value:", saddle_val, "\n")

# Show 3 contours: saddle point + 2 lowest around minima
contour_levels <- c(-0.24, -0.22, 0)
cat("Showing contours at levels:", contour_levels, "\n")

# Greyscale colors: white for smallest (minima), grey for largest (saddle)
contour_colors <- c("white", "grey70", "grey40")  # white for -0.24, lighter for -0.22, darker for 0

# Function to create sticker variant
create_final_sticker <- function(path_color = cu_gold, hex_fill = "black", filename_suffix = "") {

  # Create contour plot with greyscale on black background
  p <- ggplot(grid, aes(x = x, y = y, z = z)) +
    # Add each contour level with its specific color
    geom_contour(breaks = -0.24, color = "white", linewidth = 0.4) +
    geom_contour(breaks = -0.22, color = "grey70", linewidth = 0.4) +
    geom_contour(breaks = 0, color = "grey40", linewidth = 0.4) +
    # Add paths (no arrows, no points)
    geom_path(
      data = paths_df, aes(x = x, y = y, z = NULL, group = path_id),
      color = path_color, linewidth = 0.8, lineend = "round"
    ) +
    # Expand to fill the plot area completely
    coord_cartesian(xlim = c(-1.8, 1.8), ylim = c(-1.2, 1.2), expand = FALSE) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      plot.margin = margin(0, 0, 0, 0)
    )

  filename <- paste0("arcopt_final", filename_suffix, ".png")

  sticker(
    p,
    package = "arcopt",
    p_size = 20,
    p_color = path_color,  # CU gold text
    p_y = 1.4,
    s_x = 1,
    s_y = 0.85,
    s_width = 1.8,
    s_height = 1.4,
    h_fill = hex_fill,
    h_color = path_color,  # CU gold hex border
    h_size = 1.5,
    filename = filename,
    dpi = 300
  )

  cat("Created:", filename, "\n")
}

cat("\n=== Creating final CU Boulder hex stickers ===\n\n")

# Black background variants
cat("Black background variants:\n")
cat("1. CU Gold paths and text (black background)\n")
create_final_sticker(
  path_color = cu_gold,
  hex_fill = "black",
  filename_suffix = "_gold_black"
)

cat("2. Accessible CU Gold paths and text (black background)\n")
create_final_sticker(
  path_color = cu_gold_accessible,
  hex_fill = "black",
  filename_suffix = "_darkgold_black"
)

# Transparent background variants
cat("\nTransparent background variants:\n")
cat("3. CU Gold paths and text (transparent background)\n")
create_final_sticker(
  path_color = cu_gold,
  hex_fill = "transparent",
  filename_suffix = "_gold_transparent"
)

cat("4. Accessible CU Gold paths and text (transparent background)\n")
create_final_sticker(
  path_color = cu_gold_accessible,
  hex_fill = "transparent",
  filename_suffix = "_darkgold_transparent"
)

cat("\n=== Summary ===\n")
cat("Created 4 final CU Boulder branded hex sticker variants.\n")
cat("2 with black background, 2 with transparent background.\n")
cat("Greyscale contours (white=minima, grey=saddle).\n")
cat("CU Gold thin paths (no arrows, no points).\n")
cat("CU Gold text and hex border.\n")
cat("Clean, professional design.\n")
