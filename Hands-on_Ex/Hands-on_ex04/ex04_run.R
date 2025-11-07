## Hands-on Exercise 04 runner (Chapter 8: Spatial Weights & Applications)
## One-call script: source("Hands-on_Ex/Hands-on_ex04/ex04_run.R")

message("[ex04] Start...")

## 0) Packages ---------------------------------------------------------------
pkg_needed <- c("sf", "spdep", "tmap", "tidyverse")
installed <- rownames(installed.packages())
for (p in pkg_needed) {
  if (!(p %in% installed)) {
    message(sprintf("[ex04] Installing missing package: %s", p))
    install.packages(p, dependencies = TRUE)
  }
}
invisible(lapply(pkg_needed, library, character.only = TRUE))

tmap_mode("plot")

## 1) Paths ------------------------------------------------------------------
base_dir <- file.path("Hands-on_Ex", "Hands-on_ex04")
geo_path <- file.path(base_dir, "data", "geospatial", "Hunan.shp")
csv_path <- file.path(base_dir, "data", "aspatial", "Hunan_2012.csv")
fig_dir <- file.path(base_dir, "figures")
out_dir <- file.path(base_dir, "outputs")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

stopifnot(file.exists(geo_path))
stopifnot(file.exists(csv_path))

## 2) Data ingest & join -----------------------------------------------------
message("[ex04] Reading geospatial: ", geo_path)
hunan <- sf::st_read(geo_path, quiet = TRUE)
message(sprintf("[ex04] hunan: n=%s, cols=%s, crs=%s", nrow(hunan), ncol(hunan), st_crs(hunan)$input))

message("[ex04] Reading CSV: ", csv_path)
tbl <- readr::read_csv(csv_path, show_col_types = FALSE)
message(sprintf("[ex04] csv: n=%s, cols=%s", nrow(tbl), ncol(tbl)))

if (!"County" %in% names(hunan)) {
  stop("[ex04] 'County' field not found in shapefile attributes. Please ensure the shapefile has a 'County' column for joining.")
}
if (!all(c("County", "GDPPC") %in% names(tbl))) {
  stop("[ex04] CSV must contain 'County' and 'GDPPC' fields.")
}

hunan <- dplyr::left_join(hunan, tbl, by = "County")
if (nrow(hunan) != 88) warning("[ex04] Joined features not 88. Actual: ", nrow(hunan))

## 3) Quick GDPPC choropleth -------------------------------------------------
message("[ex04] Export GDPPC choropleth...")
tm1 <- tm_shape(hunan) +
  tm_polygons("GDPPC", palette = "YlOrRd", style = "quantile", title = "GDPPC") +
  tm_layout(frame = FALSE, legend.outside = TRUE)

tmap::tmap_save(tm1, filename = file.path(fig_dir, "ex04_gdppc.png"), width = 2000, height = 1400, units = "px")

## 4) Neighbors (Queen/Rook) -------------------------------------------------
message("[ex04] Computing Queen/Rook neighbors...")
nb_q <- spdep::poly2nb(hunan, queen = TRUE)
nb_r <- spdep::poly2nb(hunan, queen = FALSE)

card_q <- spdep::card(nb_q)
card_r <- spdep::card(nb_r)

summarize_nb <- function(nb, label) {
  cards <- spdep::card(nb)
  data.frame(
    type = label,
    n_obs = length(nb),
    total_links = sum(cards),
    mean_neighbors = mean(cards),
    min_neighbors = min(cards),
    max_neighbors = max(cards),
    islands = sum(cards == 0),
    stringsAsFactors = FALSE
  )
}

## 5) Coordinates & lines for plotting --------------------------------------
coords <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(hunan)))

to_lines_sf <- function(nb, coords) {
  ## Prefer spdep::nb2lines(as_sf=TRUE); fallback manual if needed
  ln <- try(spdep::nb2lines(nb, coords, as_sf = TRUE), silent = TRUE)
  if (inherits(ln, "try-error")) {
    segs <- lapply(seq_along(nb), function(i) {
      if (length(nb[[i]]) == 0) return(NULL)
      do.call(sf::st_sfc, lapply(nb[[i]], function(j) {
        sf::st_linestring(rbind(coords[i, ], coords[j, ]))
      }))
    })
    segs <- do.call(c, segs)
    sf::st_as_sf(data.frame(id = seq_along(segs)), geometry = segs)
  } else ln
}

lines_q <- to_lines_sf(nb_q, coords)
lines_r <- to_lines_sf(nb_r, coords)
pts <- sf::st_as_sf(as.data.frame(coords), coords = c("X","Y"), crs = 4326)

## 6) Plot adjacency networks ------------------------------------------------
message("[ex04] Export Queen/Rook network plots...")
tm_q <- tm_shape(hunan) + tm_polygons(col = "white", border = "grey80") +
  tm_shape(lines_q) + tm_lines(col = "steelblue", lwd = 1) +
  tm_shape(pts) + tm_dots(size = 0.03, col = "darkred") +
  tm_layout(title = "Queen contiguity", frame = FALSE)

tm_r <- tm_shape(hunan) + tm_polygons(col = "white", border = "grey80") +
  tm_shape(lines_r) + tm_lines(col = "seagreen", lwd = 1) +
  tm_shape(pts) + tm_dots(size = 0.03, col = "darkred") +
  tm_layout(title = "Rook contiguity", frame = FALSE)

tmap::tmap_save(tm_q, filename = file.path(fig_dir, "ex04_queen.png"), width = 2000, height = 1400, units = "px")
tmap::tmap_save(tm_r, filename = file.path(fig_dir, "ex04_rook.png"), width = 2000, height = 1400, units = "px")

## 7) Distance: nearest neighbor, fixed 62 km, K=6 ---------------------------
message("[ex04] Computing distance-based neighbors...")
kn1 <- spdep::knearneigh(coords, k = 1, longlat = TRUE)
nb_kn1 <- spdep::knn2nb(kn1)
nn1_dists <- spdep::nbdists(nb_kn1, coords, longlat = TRUE)
max_nn_km <- max(unlist(nn1_dists))
writeLines(sprintf("max_nearest_neighbor_km,%.6f", max_nn_km), con = file.path(out_dir, "ex04_nn_maxdist.txt"))

nb_d62 <- spdep::dnearneigh(coords, d1 = 0, d2 = 62, longlat = TRUE)
nb_kn6 <- spdep::knn2nb(spdep::knearneigh(coords, k = 6, longlat = TRUE))

lines_d62 <- to_lines_sf(nb_d62, coords)
lines_kn6 <- to_lines_sf(nb_kn6, coords)

tm_d62 <- tm_shape(hunan) + tm_polygons(col = "white", border = "grey80") +
  tm_shape(lines_d62) + tm_lines(col = "orange", lwd = 1) +
  tm_shape(pts) + tm_dots(size = 0.03, col = "black") +
  tm_layout(title = "Fixed distance (<= 62 km)", frame = FALSE)

tm_kn6 <- tm_shape(hunan) + tm_polygons(col = "white", border = "grey80") +
  tm_shape(lines_kn6) + tm_lines(col = "purple", lwd = 1) +
  tm_shape(pts) + tm_dots(size = 0.03, col = "black") +
  tm_layout(title = "KNN (k = 6)", frame = FALSE)

tmap::tmap_save(tm_d62, filename = file.path(fig_dir, "ex04_dist62.png"), width = 2000, height = 1400, units = "px")
tmap::tmap_save(tm_kn6, filename = file.path(fig_dir, "ex04_knn6.png"), width = 2000, height = 1400, units = "px")

## 8) Summaries --------------------------------------------------------------
sum_q   <- summarize_nb(nb_q,   "queen")
sum_r   <- summarize_nb(nb_r,   "rook")
sum_d62 <- summarize_nb(nb_d62, "dist62km")
sum_kn6 <- summarize_nb(nb_kn6, "knn6")
nb_summ <- dplyr::bind_rows(sum_q, sum_r, sum_d62, sum_kn6)
readr::write_csv(nb_summ, file.path(out_dir, "ex04_neighbors_summary.csv"))

## 9) Row-standardized lag, binary-sum lag ----------------------------------
message("[ex04] Computing lags (mean/sum) and window stats...")
rswm_q <- spdep::nb2listw(nb_q, style = "W", zero.policy = TRUE)
lag_mean_q <- spdep::lag.listw(rswm_q, hunan$GDPPC, zero.policy = TRUE)

bw_q <- spdep::nb2listw(nb_q, style = "B", zero.policy = TRUE)
lag_sum_q <- spdep::lag.listw(bw_q, hunan$GDPPC, zero.policy = TRUE)

nb_q_self <- spdep::include.self(nb_q)
rswm_q_self <- spdep::nb2listw(nb_q_self, style = "W", zero.policy = TRUE)
win_avg <- spdep::lag.listw(rswm_q_self, hunan$GDPPC, zero.policy = TRUE)

bw_q_self <- spdep::nb2listw(nb_q_self, style = "B", zero.policy = TRUE)
win_sum <- spdep::lag.listw(bw_q_self, hunan$GDPPC, zero.policy = TRUE)

hunan$lag_mean_q <- lag_mean_q
hunan$lag_sum_q  <- lag_sum_q
hunan$win_avg    <- win_avg
hunan$win_sum    <- win_sum

## 10) Lag comparison plots --------------------------------------------------
tm_lag1 <- tm_shape(hunan) + tm_polygons("GDPPC", palette = "Blues", style = "quantile", title = "GDPPC") + tm_layout(frame = FALSE)
tm_lag2 <- tm_shape(hunan) + tm_polygons("lag_mean_q", palette = "Greens", style = "quantile", title = "Lag mean (W)") + tm_layout(frame = FALSE)
lag_grid <- tmap::tmap_arrange(tm_lag1, tm_lag2, ncol = 2)
tmap::tmap_save(lag_grid, filename = file.path(fig_dir, "ex04_lag_compare.png"), width = 2200, height = 1200, units = "px")

tm_win1 <- tm_shape(hunan) + tm_polygons("win_avg", palette = "Oranges", style = "quantile", title = "Window avg (incl. self)") + tm_layout(frame = FALSE)
tm_win2 <- tm_shape(hunan) + tm_polygons("win_sum", palette = "Reds", style = "quantile", title = "Window sum (incl. self)") + tm_layout(frame = FALSE)
win_grid <- tmap::tmap_arrange(tm_win1, tm_win2, ncol = 2)
tmap::tmap_save(win_grid, filename = file.path(fig_dir, "ex04_window_compare.png"), width = 2200, height = 1200, units = "px")

## 11) Case: Anxiang ---------------------------------------------------------
idx_anx <- which(hunan$County == "Anxiang")
if (length(idx_anx) == 1) {
  nb_ids <- nb_q[[idx_anx]]
  nb_names <- hunan$County[nb_ids]
  nb_gdppc <- hunan$GDPPC[nb_ids]
  anx_lag_mean <- lag_mean_q[idx_anx]
  case_df <- tibble::tibble(
    county = "Anxiang",
    neighbor_index = nb_ids,
    neighbor_name = nb_names,
    neighbor_GDPPC = nb_gdppc,
    lag_mean_q = rep(anx_lag_mean, length(nb_ids))
  )
  readr::write_csv(case_df, file.path(out_dir, "ex04_case_anxiang.csv"))
} else {
  warning("[ex04] 'Anxiang' not found uniquely in County. Skipping case export.")
}

message("[ex04] Done. Figures in ", fig_dir, "; tables in ", out_dir)

