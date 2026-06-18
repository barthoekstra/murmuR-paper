# Animated map of simulated autumn nocturnal bird migration over the IBM data
# domain (window set by win_start/win_end; default 15 Sep - 15 Nov 2018).
#
# Features
#   * Realistic Earth basemap. With a MapTiler key (env MAPTILER_API_KEY) the
#     day side uses "Satellite Plain" (satellite-v4) and the night side
#     "Satellite Plain Dark" (satellite-v4-dark), blended across the solar
#     terminator. Without a key it falls back to Esri World Imagery + a smooth
#     navy night-gradient overlay.
#   * Realistic moving solar terminator (suncalc); smooth day -> dusk -> night.
#   * ERA5 850 hPa wind field as vectors (the IBM's flight level), coloured by
#     speed; linearly interpolated between hours.
#   * Simulated migrants linearly interpolated to FRAME_STEP_MIN resolution so
#     motion is smooth at sub-hourly cadence.
#   * Birds bound for the Netherlands light up (highlight colour) the moment
#     their track crosses NL; their source/take-off regions are emphasised
#     (SOURCE_MODE: "glow" soft origin halos, or "trails" longer dark trails).
#
# Data (read-only): tracks processed/tracks/seasons/2018_autumn.RDS (full
# domain); winds ERA5 pressure-level 850 hPa from the original ibm-ml tree.
# Bird states: -1 inactive (never drawn); 0 grounded; 1 flying.
#
# MODE "design" -> short review sequence (incl. sub-hour frames) + preview GIF;
# MODE "full"   -> every FRAME_STEP_MIN frame of the window -> MP4.

suppressPackageStartupMessages({
  library(data.table); library(sf); library(terra); library(maptiles)
  library(ggplot2); library(rnaturalearth); library(ncdf4); library(suncalc)
  library(ragg); library(gifski); library(patchwork); library(systemfonts)
})

# Register Open Sans (downloaded into figures/fonts) for ragg text rendering.
register_font("Open Sans",
  plain = "figures/fonts/OpenSans-Regular.ttf",
  bold  = "figures/fonts/OpenSans-Bold.ttf")
FONT <- "Open Sans"

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
MODE          <- Sys.getenv("ANIM_MODE", "design")   # design | clip | full

data_path     <- "/data/birdcloudstorage-tvm/ibm-ml-refactor/data/"
wind_dir      <- "/data/birdcloudstorage-tvm/ibm-ml/data/ecmwf-era5-pressurelevels/2018"
wind_level    <- 850L                       # hPa (== murmuR::ibm_pressurelevel)
tracks_file   <- file.path(data_path, "processed/tracks/seasons/2018_autumn.RDS")
cache_file    <- "figures/.cache_movers_sep15_nov15.RDS"
nl_cache      <- "figures/.cache_nl_polygon.RDS"
esri_cache    <- "figures/.cache_basemap_esri.RDS"
mt_cache      <- "figures/.cache_basemap_maptiler.RDS"

# Balanced wider 16:9 view (lon full; lat cropped to ~10 deg around NL).
view          <- c(xmin = -10, ymin = 48, xmax = 20, ymax = 58)
radar_lon     <- 5.1381; radar_lat <- 51.8369

win_start     <- as.POSIXct("2018-09-15 00:00", tz = "UTC")
win_end       <- as.POSIXct("2018-11-15 23:00", tz = "UTC")
FRAME_STEP_MIN<- 10L

wind_stride   <- 5L
arrow_scale   <- 0.045                       # deg per (m/s) (850 hPa is faster)

trail_K       <- 3L                          # hours of fading trail (short, quick fade)
trail_K_nl    <- 6L                          # longer trail for NL birds in "trails" mode
trail_gamma   <- 1.9                          # >1 = fade out quicker
SUBSAMPLE_FRAC<- 0.5
subsample_seed<- 42L
SOURCE_MODE   <- "off"                        # "off" | "glow" | "trails"

night_darken  <- 0.60                         # multiply the dark basemap (smaller = darker night)
term_step     <- 0.5
maptiler_key  <- Sys.getenv("MAPTILER_API_KEY")
sat_zoom      <- 6L

phen_cache    <- "figures/.cache_phenology_2018autumn.RDS"

# Smooth preview clip (MODE = "clip"): a full dusk -> night -> dawn peak night.
clip_from     <- "2018-10-31 22:00"
clip_to       <- "2018-11-01 02:00"

# NL polygon crop (European NL only)
nl_bbox       <- st_bbox(c(xmin = 3, ymin = 50, xmax = 8, ymax = 54), crs = 4326)

design_hours <- as.POSIXct(c(
  "2018-10-18 12:00", "2018-10-18 17:00",
  "2018-10-18 19:00", "2018-10-18 19:20", "2018-10-18 19:40",  # sub-hour interp
  "2018-10-18 22:00", "2018-10-19 06:00"
), tz = "UTC")

# ---------------------------------------------------------------------------
# Palette
# ---------------------------------------------------------------------------
pal <- list(
  border = "#cdd6e0",
  arrow_low = "#88b0d6", arrow_high = "#eaf6ff", arrow_alpha = 0.5,
  trail = "#ff7a2f", fly_core = "#ffe07a", fly_halo = "#ff9d33",
  hl_trail = "#36e0ff", hl_core = "#ffffff", hl_halo = "#36e0ff",
  ground = "#9fb0c8", radar = "#ff4d4d",
  src_glow = "#ffd27f",
  text = "white", text_dim = "grey75"
)
# Smooth night gradient (Esri fallback): solar-altitude -> rgba
term_alt_stops <- c( 6,    1,    -3,    -8,    -14,   -40)
term_col_stops <- c("#0a142800","#ff9a4d22","#b5642e4d","#243a6699","#0d1a3ad9","#070e22f2")

# ---------------------------------------------------------------------------
# tidx <-> datetime
# ---------------------------------------------------------------------------
dts  <- seq.POSIXt(as.POSIXct("2018-08-01 00:00", tz = "UTC"),
                   as.POSIXct("2018-12-01 00:00", tz = "UTC") - 3600, by = "1 hour")
tmap <- data.table(tidx = seq_along(dts), datetime = dts); setkey(tmap, datetime)
tidx_of <- function(dt) tmap[.(dt), tidx]

# ---------------------------------------------------------------------------
# Movers subset + subsample
# ---------------------------------------------------------------------------
if (!file.exists(cache_file)) {
  message("Building movers cache ...")
  ta <- tidx_of(win_start); tb <- tidx_of(win_end)
  s  <- as.data.table(readRDS(tracks_file))
  movers <- s[tidx >= ta & tidx <= tb & state == 1L, unique(bird)]
  sub <- s[bird %in% movers & tidx >= (ta - 9L) & tidx <= tb & state != -1L,
           .(tidx, bird, x, y, state)]
  sub <- merge(sub, tmap, by = "tidx"); setkey(sub, tidx, bird)
  saveRDS(list(sub = sub, movers = movers), cache_file); rm(s); gc()
}
C   <- readRDS(cache_file); sub <- C$sub
set.seed(subsample_seed)
keep_birds <- sample(C$movers, ceiling(length(C$movers) * SUBSAMPLE_FRAC))
sub <- sub[bird %in% keep_birds]; setkey(sub, tidx, bird)

# Per-bird sub-hour phase offset in [0,1): desynchronises the hourly quantisation
# so departures (and the aloft count) spread smoothly across each hour instead
# of all snapping to :00. Each bird's whole timeline is shifted by its phase.
set.seed(subsample_seed + 1L)
bird_phase <- data.table(bird = keep_birds, phi = runif(length(keep_birds)))
setkey(bird_phase, bird)

# ---------------------------------------------------------------------------
# NL crossing + source origins (precompute once)
# ---------------------------------------------------------------------------
nl_poly <- suppressWarnings(st_union(st_crop(st_sf(geometry = readRDS(nl_cache)), nl_bbox)))
pts_all <- unique(sub[, .(bird, tidx, x, y)])
inside  <- st_intersects(st_as_sf(pts_all, coords = c("x", "y"), crs = 4326),
                         nl_poly, sparse = FALSE)[, 1]
pts_all[, in_nl := inside]
cross_info <- pts_all[in_nl == TRUE, .(crossed_tidx = min(tidx)), by = bird]
setkey(cross_info, bird)
nl_crossers <- cross_info$bird
origins <- sub[bird %in% nl_crossers][order(bird, tidx),
              .(ox = x[1], oy = y[1], born = tidx[1]), by = bird]
cat(sprintf("Drawn movers: %d | NL-crossers: %d\n",
            length(keep_birds), length(nl_crossers)))

# ---------------------------------------------------------------------------
# Basemap
# ---------------------------------------------------------------------------
view_sfc <- st_as_sfc(st_bbox(c(xmin = view[["xmin"]], ymin = view[["ymin"]],
                                xmax = view[["xmax"]], ymax = view[["ymax"]]), crs = 4326))
fetch_tiles_array <- function(provider_obj) {
  r <- get_tiles(view_sfc, provider = provider_obj, zoom = sat_zoom,
                 crop = TRUE, cachedir = "/tmp/tilecache")
  r <- terra::project(r, "EPSG:4326")
  list(arr = terra::as.array(r) / 255,
       ext = as.vector(terra::ext(r)))   # xmin,xmax,ymin,ymax
}
USE_MAPTILER <- nzchar(maptiler_key)
if (USE_MAPTILER) {
  if (!file.exists(mt_cache)) {
    message("Fetching MapTiler satellite (day + dark) ...")
    day_p  <- create_provider("mt_day",
      sprintf("https://api.maptiler.com/maps/satellite-v4/256/{z}/{x}/{y}.jpg?key=%s", maptiler_key),
      citation = "(c) MapTiler (c) OpenStreetMap contributors")
    dark_p <- create_provider("mt_dark",
      sprintf("https://api.maptiler.com/maps/satellite-v4-dark/256/{z}/{x}/{y}.jpg?key=%s", maptiler_key),
      citation = "(c) MapTiler (c) OpenStreetMap contributors")
    day  <- fetch_tiles_array(day_p)
    dark <- fetch_tiles_array(dark_p)
    saveRDS(list(day = day, dark = dark), mt_cache)
  }
  BM <- readRDS(mt_cache)
  base_ext <- BM$day$ext
} else {
  if (!file.exists(esri_cache)) {
    message("Fetching Esri World Imagery basemap ...")
    e <- fetch_tiles_array("Esri.WorldImagery")
    saveRDS(e, esri_cache)
  }
  BM <- readRDS(esri_cache)
  base_ext <- BM$ext
}
# Pixel lon/lat for terminator blending (image rows N->S, cols W->E)
img_dim <- if (USE_MAPTILER) dim(BM$day$arr) else dim(BM$arr)
px_lon  <- seq(base_ext[1], base_ext[2], length.out = img_dim[2])
px_lat  <- seq(base_ext[4], base_ext[3], length.out = img_dim[1])   # top=north

borders <- ne_download(scale = 10, category = "cultural",
            type = "admin_0_boundary_lines_land", returnclass = "sf")
borders <- suppressWarnings(st_crop(st_make_valid(borders),
            st_bbox(c(xmin = view[["xmin"]], ymin = view[["ymin"]],
                      xmax = view[["xmax"]], ymax = view[["ymax"]]), crs = 4326)))
radar_pt <- st_sf(geometry = st_sfc(st_point(c(radar_lon, radar_lat)), crs = 4326))

# ---------------------------------------------------------------------------
# Solar terminator (coarse grid for Esri gradient; pixel weight for blend)
# ---------------------------------------------------------------------------
term_grid <- as.data.table(expand.grid(
  lon = seq(view[["xmin"]], view[["xmax"]], by = term_step),
  lat = seq(view[["ymin"]], view[["ymax"]], by = term_step)))
terminator_grid <- function(dt) {
  p <- getSunlightPosition(data = data.frame(date = dt, lat = term_grid$lat, lon = term_grid$lon))
  d <- copy(term_grid); d[, alt := p$altitude * 180 / pi]; d
}
# Per-pixel daylight weight (1 day -> 0 night) via a coarse solar-altitude
# grid bilinearly upsampled to the basemap pixel grid.
wt_lon <- seq(min(px_lon), max(px_lon), length.out = 90)   # increasing
wt_lat <- seq(min(px_lat), max(px_lat), length.out = 70)   # increasing
wt_g   <- expand.grid(lon = wt_lon, lat = wt_lat)          # lon fastest
cx0 <- findInterval(px_lon, wt_lon, all.inside = TRUE)
cy0 <- findInterval(px_lat, wt_lat, all.inside = TRUE)     # px_lat is N->S
fcx <- (px_lon - wt_lon[cx0]) / (wt_lon[cx0 + 1] - wt_lon[cx0])
fcy <- (px_lat - wt_lat[cy0]) / (wt_lat[cy0 + 1] - wt_lat[cy0])
n_pxlon <- length(px_lon); n_pxlat <- length(px_lat)
fcy_row <- matrix(1 - fcy, nrow = n_pxlon, ncol = n_pxlat, byrow = TRUE)
fcy_rowc<- matrix(fcy,     nrow = n_pxlon, ncol = n_pxlat, byrow = TRUE)
daylight_weight <- function(dt) {
  p <- getSunlightPosition(data = data.frame(date = dt, lat = wt_g$lat, lon = wt_g$lon))
  W <- matrix(pmin(1, pmax(0, (p$altitude * 180 / pi + 6) / 12)),
              nrow = length(wt_lon), ncol = length(wt_lat))   # W[lon, lat]
  # separable bilinear: interpolate along lon (rows), then along lat (cols)
  Wc <- W[cx0, , drop = FALSE] * (1 - fcx) + W[cx0 + 1, , drop = FALSE] * fcx  # [px_lon, wt_lat]
  out_t <- Wc[, cy0, drop = FALSE] * fcy_row + Wc[, cy0 + 1, drop = FALSE] * fcy_rowc
  t(out_t)                                              # [px_lat (N->S), px_lon]
}
blend_basemap <- function(dt) {
  wa <- array(daylight_weight(dt), dim = c(img_dim[1], img_dim[2], 3))
  grDevices::as.raster(BM$day$arr * wa + (BM$dark$arr * night_darken) * (1 - wa))
}

# ---------------------------------------------------------------------------
# Winds: 850 hPa, decimated, memoised per hour, linearly interpolated.
# The window can span several monthly ERA5 files, so open every month it
# touches and build one datetime -> (file, time index) map.
# ---------------------------------------------------------------------------
win_months <- sort(unique(format(seq(win_start, win_end, by = "day"), "%Y%m")))
wind_files <- file.path(wind_dir, sprintf("%s.nc", win_months))
ncs <- lapply(wind_files, nc_open)
nc1 <- ncs[[1]]
nc_lon  <- nc1$dim[[grep("^lon", names(nc1$dim), ignore.case = TRUE)[1]]]$vals
nc_lat  <- nc1$dim[[grep("^lat", names(nc1$dim), ignore.case = TRUE)[1]]]$vals
nc_lev  <- nc1$dim[[grep("level", names(nc1$dim), ignore.case = TRUE)[1]]]$vals
lev_idx <- which(nc_lev == wind_level)
time_map <- rbindlist(lapply(seq_along(ncs), function(fi) {
  nc <- ncs[[fi]]
  tn <- grep("time", names(nc$dim), ignore.case = TRUE, value = TRUE)[1]
  torig <- as.POSIXct(sub(".*since *", "", nc$dim[[tn]]$units), tz = "UTC")
  data.table(datetime = torig + nc$dim[[tn]]$vals * 3600, fi = fi,
             ti = seq_along(nc$dim[[tn]]$vals))
}))
setkey(time_map, datetime)
lon_keep <- seq(1, length(nc_lon), by = wind_stride)
lat_keep <- seq(1, length(nc_lat), by = wind_stride)
wgrid    <- CJ(j = lat_keep, i = lon_keep); wgrid[, `:=`(lon = nc_lon[i], lat = nc_lat[j])]
.wcache  <- new.env()
wind_grid <- function(tidx) {
  key <- as.character(tidx)
  if (!is.null(.wcache[[key]])) return(.wcache[[key]])
  m <- time_map[.(dts[tidx])]
  if (is.na(m$fi)) m <- time_map[which.min(abs(as.numeric(datetime) - as.numeric(dts[tidx])))]
  u <- ncvar_get(ncs[[m$fi]], "u", start = c(1, 1, lev_idx, m$ti), count = c(-1, -1, 1, 1))
  v <- ncvar_get(ncs[[m$fi]], "v", start = c(1, 1, lev_idx, m$ti), count = c(-1, -1, 1, 1))
  d <- copy(wgrid); d[, u := u[cbind(i, j)]][, v := v[cbind(i, j)]]
  .wcache[[key]] <- d; d
}
wind_interp <- function(dt) {
  t0 <- tidx_of(as.POSIXct(trunc(dt, "hours"), tz = "UTC"))
  f  <- as.numeric(difftime(dt, dts[t0], units = "hours"))
  g0 <- wind_grid(t0); g1 <- if (f > 0) wind_grid(t0 + 1L) else g0
  d <- copy(g0); d[, u := g0$u + f * (g1$u - g0$u)][, v := g0$v + f * (g1$v - g0$v)]
  d[, ws := sqrt(u^2 + v^2)]
  d[, lon_end := lon + (u * arrow_scale) / cos(lat * pi / 180)][, lat_end := lat + v * arrow_scale]
  d[]
}

# ---------------------------------------------------------------------------
# Bird state at an arbitrary (interpolated) time
# ---------------------------------------------------------------------------
frame_birds <- function(dt) {
  # continuous tidx of this (sub-hourly) frame; tidx 1 == dts[1]
  t_h <- 1 + as.numeric(difftime(dt, dts[1], units = "hours"))
  # per-bird shifted time -> per-bird floor hour t0 and fraction f
  B <- copy(bird_phase)
  B[, tau := t_h - phi]
  B[, t0 := as.integer(floor(tau))]
  B[, f := tau - t0]

  # current-hour state/position (only birds active at their own t0 are present)
  qcur <- B[, .(bird, tidx = t0, f, bt0 = t0)]
  a <- sub[qcur, on = .(bird, tidx), nomatch = 0L,
           .(bird, x, y, state, f = i.f, t0 = i.bt0)]
  nxt <- sub[a[, .(bird, tidx = t0 + 1L)], on = .(bird, tidx), nomatch = NA,
             .(bird, xb = x, yb = y)]
  a <- merge(a, nxt, by = "bird", all.x = TRUE)
  a[!is.na(xb), `:=`(x = x + f * (xb - x), y = y + f * (yb - y))]
  fly <- a[state == 1L]

  # Smooth trails: every vertex is the bird's position interpolated at the SAME
  # fraction f, so the whole trail slides forward continuously (no hourly
  # stutter) and stays continuous across hour boundaries.
  K  <- trail_K
  Fd <- fly[, .(bird, t0, f)]
  vts <- vector("list", K + 1L)
  for (j in 0:K) {
    xl <- sub[Fd[, .(bird, tidx = t0 - j, f)], on = .(bird, tidx), nomatch = 0L,
              .(bird, xl = x, yl = y, f = i.f)]
    xh <- sub[Fd[, .(bird, tidx = t0 - j + 1L)], on = .(bird, tidx), nomatch = 0L,
              .(bird, xh = x, yh = y)]
    m <- merge(xl, xh, by = "bird")
    m[, `:=`(vx = xl + f * (xh - xl), vy = yl + f * (yh - yl), j = j)]
    vts[[j + 1L]] <- m[, .(bird, j, vx, vy)]
  }
  V <- rbindlist(vts); setorder(V, bird, j)
  V[, `:=`(vx1 = data.table::shift(vx, -1L), vy1 = data.table::shift(vy, -1L),
           j1 = data.table::shift(j, -1L)), by = bird]
  seg <- V[!is.na(j1) & j1 == j + 1L]
  seg[, age := j][, w := pmax(0, (1 - age / (K + 1))) ^ trail_gamma]

  list(fly = fly, seg = seg, t0 = as.integer(floor(t_h)),
       aloft = round(nrow(fly) / SUBSAMPLE_FRAC))   # scale subsample -> est. true aloft
}

# Source glow: origins of NL-birds that have departed by t0.
source_layer <- function(t0) origins[born <= t0]

# ---------------------------------------------------------------------------
# Inset: seasonal phenology of birds aloft with a moving current-time bar
# ---------------------------------------------------------------------------
PHEN <- as.data.table(readRDS(phen_cache))   # datetime, n_aloft over full season
# light smoothing so the curve reads as phenology, not hourly spikes
PHEN[, n_smooth := frollmean(n_aloft, 5, align = "center", fill = NA)]
PHEN[is.na(n_smooth), n_smooth := n_aloft]
phen_ymax <- max(PHEN$n_smooth)

inset_plot <- function(dt) {
  ggplot(PHEN, aes(datetime, n_smooth)) +
    annotate("rect", xmin = win_start, xmax = win_end, ymin = 0, ymax = phen_ymax,
      fill = "#ffffff", alpha = 0.07) +
    geom_area(fill = "#ffb24d", alpha = 0.35) +
    geom_line(colour = "#ffd27f", linewidth = 0.35) +
    geom_vline(xintercept = as.numeric(as.POSIXct(dt, tz = "UTC")),
      colour = "#ffffff", linewidth = 0.6) +
    scale_x_datetime(date_labels = "%b", date_breaks = "1 month", expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(title = "Birds aloft (autumn season)", x = NULL, y = NULL) +
    theme_minimal(base_size = 9, base_family = FONT) +
    theme(
      plot.background  = element_rect(fill = scales::alpha("black", 0.42), colour = NA),
      panel.background = element_rect(fill = NA, colour = NA),
      panel.grid       = element_blank(),
      plot.title       = element_text(colour = "white", size = 8, hjust = 0),
      axis.text.x      = element_text(colour = "grey80", size = 7),
      axis.text.y      = element_blank(),
      plot.margin      = margin(4, 6, 2, 6))
}

# ---------------------------------------------------------------------------
# Frame plot
# ---------------------------------------------------------------------------
make_frame <- function(dt) {
  w  <- wind_interp(dt)
  fb <- frame_birds(dt)
  src<- source_layer(fb$t0)
  date_str <- paste0(as.integer(format(dt, "%d")), " ",
                     format(dt, "%B %Y  %H:%M"), " UTC")   # e.g. 18 October 2018  22:00 UTC

  p <- ggplot()
  if (USE_MAPTILER) {
    p <- p + annotation_raster(blend_basemap(dt),
      xmin = base_ext[1], xmax = base_ext[2], ymin = base_ext[3], ymax = base_ext[4],
      interpolate = TRUE)
  } else {
    p <- p + annotation_raster(grDevices::as.raster(BM$arr),
      xmin = base_ext[1], xmax = base_ext[2], ymin = base_ext[3], ymax = base_ext[4],
      interpolate = TRUE)
  }
  # source-region glow (origin/take-off areas of NL-bound birds; builds over time)
  if (SOURCE_MODE == "glow" && nrow(src) > 0)
    p <- p +
      geom_point(data = src, aes(ox, oy), colour = pal$src_glow, size = 5.5,
                 alpha = 0.07, shape = 16) +
      geom_point(data = src, aes(ox, oy), colour = pal$src_glow, size = 2.6,
                 alpha = 0.13, shape = 16) +
      geom_point(data = src, aes(ox, oy), colour = pal$src_glow, size = 1.0,
                 alpha = 0.30, shape = 16)
  # (grounded / non-migrating birds are intentionally NOT drawn)
  # Esri night gradient (skipped for MapTiler blend)
  if (!USE_MAPTILER) {
    tg <- terminator_grid(dt)
    p <- p + geom_raster(data = tg, aes(lon, lat, fill = alt), interpolate = TRUE) +
      scale_fill_gradientn(colours = term_col_stops, values = scales::rescale(term_alt_stops),
        limits = range(term_alt_stops), oob = scales::squish, guide = "none")
  }
  p <- p +
    geom_sf(data = borders, colour = pal$border, linewidth = 0.15, alpha = 0.45) +
    geom_segment(data = w, aes(lon, lat, xend = lon_end, yend = lat_end, colour = ws),
      arrow = arrow(length = unit(0.028, "inches"), type = "closed"),
      linewidth = 0.25, alpha = pal$arrow_alpha) +
    scale_colour_gradient(low = pal$arrow_low, high = pal$arrow_high,
      limits = ws_rng, oob = scales::squish, guide = "none") +
    # flight trails (short, quick fade)
    geom_segment(data = fb$seg, aes(vx, vy, xend = vx1, yend = vy1, alpha = w),
      colour = pal$trail, linewidth = 0.24, lineend = "round") +
    scale_alpha_identity() +
    # flyers: halo + bright core (single warm colour)
    geom_point(data = fb$fly, aes(x, y), colour = pal$fly_halo,
      size = 1.5, alpha = 0.16, shape = 16) +
    geom_point(data = fb$fly, aes(x, y), colour = pal$fly_core,
      size = 0.5, alpha = 0.95, shape = 16) +
    geom_sf(data = radar_pt, shape = 21, size = 1.9, stroke = 0.6,
      colour = pal$radar, fill = NA) +
    coord_sf(xlim = view[c("xmin", "xmax")], ylim = view[c("ymin", "ymax")],
      expand = FALSE, crs = 4326) +
    # --- top-left info panel ---
    annotate("rect", xmin = view[["xmin"]], xmax = view[["xmin"]] + 10.6,
      ymin = view[["ymax"]] - 3.35, ymax = view[["ymax"]],
      fill = "black", alpha = 0.40) +
    annotate("text", x = view[["xmin"]] + 0.55, y = view[["ymax"]] - 0.55,
      label = date_str, hjust = 0, vjust = 1,
      colour = pal$text, family = FONT, fontface = "bold", size = 5.9) +
    annotate("text", x = view[["xmin"]] + 0.55, y = view[["ymax"]] - 1.5,
      label = "Simulated nocturnal bird migration · 850 hPa wind",
      hjust = 0, vjust = 1, colour = pal$text_dim, family = FONT, size = 3.3) +
    annotate("text", x = view[["xmin"]] + 0.55, y = view[["ymax"]] - 2.28,
      label = format(fb$aloft, big.mark = ","), hjust = 0, vjust = 1,
      colour = pal$fly_core, family = FONT, fontface = "bold", size = 8.4) +
    annotate("text", x = view[["xmin"]] + 0.62, y = view[["ymax"]] - 3.02,
      label = "migrants aloft", hjust = 0, vjust = 1,
      colour = pal$text_dim, family = FONT, size = 3.5) +
    labs(x = NULL, y = NULL) + theme_void() +
    theme(plot.background = element_rect(fill = "black", colour = NA),
          plot.margin = margin(1, 1, 1, 1))

  # compose with the seasonal phenology inset (bottom-right)
  p + inset_element(inset_plot(dt), left = 0.635, bottom = 0.04,
                    right = 0.99, top = 0.30, align_to = "panel")
}

# ---------------------------------------------------------------------------
# Wind colour range across the window (consistent colour)
# ---------------------------------------------------------------------------
rng_hours <- seq.POSIXt(win_start, win_end, by = "3 hours")
ws_rng <- range(vapply(rng_hours, function(h) range(wind_grid(tidx_of(h))[, sqrt(u^2 + v^2)],
                                                    na.rm = TRUE), numeric(2)))

# ---------------------------------------------------------------------------
# Render
# ---------------------------------------------------------------------------
midlat   <- mean(view[c("ymin", "ymax")])
geo_aspc <- (view[["xmax"]] - view[["xmin"]]) * cos(midlat * pi / 180) /
            (view[["ymax"]] - view[["ymin"]])
DPI <- 150
W_PX <- 1920L
H_PX <- 2L * round(W_PX / geo_aspc / 2)        # force EVEN height (libx264/yuv420p)
W_IN <- W_PX / DPI; H_IN <- H_PX / DPI
FPS  <- 30L                                  # playback frame rate
ffmpeg_bin <- path.expand("~/bin/ffmpeg")
out_mp4    <- "figures/migration_sep15_nov15_2018.mp4"
cat(sprintf("Basemap: %s | view aspect %.3f -> %.0fx%.0f px | wind range [%.1f,%.1f]\n",
            ifelse(USE_MAPTILER, "MapTiler blend", "Esri+gradient"),
            geo_aspc, W_IN * DPI, H_IN * DPI, ws_rng[1], ws_rng[2]))

render_design <- function() {
  dir.create("figures/anim_frames", showWarnings = FALSE)
  files <- character(length(design_hours))
  for (k in seq_along(design_hours)) {
    f <- sprintf("figures/anim_frames/design_%02d.png", k)
    agg_png(f, width = W_IN, height = H_IN, units = "in", res = DPI)
    print(make_frame(design_hours[k])); dev.off()
    files[k] <- f; cat("wrote", f, "\n")
  }
  gifski(files, "figures/preview_daynight.gif",
         width = round(W_IN * DPI), height = round(H_IN * DPI), delay = 0.7, progress = FALSE)
  cat("wrote figures/preview_daynight.gif\n")
}

render_full <- function() {
  all_t <- seq.POSIXt(win_start, win_end, by = sprintf("%d min", FRAME_STEP_MIN))
  # Pre-read every hourly wind grid into .wcache BEFORE forking so the parallel
  # workers never touch the (non-fork-safe) netCDF handle.
  for (tt in tidx_of(win_start):tidx_of(win_end)) wind_grid(tt)
  fdir <- "figures/anim_frames/full"; unlink(fdir, recursive = TRUE); dir.create(fdir, recursive = TRUE)
  ncores <- max(1L, min(12L, parallel::detectCores() - 2L))
  cat(sprintf("Full render: %d frames on %d cores\n", length(all_t), ncores))
  t0 <- Sys.time()
  res <- parallel::mclapply(seq_along(all_t), function(k) {
    f <- sprintf("%s/frame_%05d.png", fdir, k)
    agg_png(f, width = W_IN, height = H_IN, units = "in", res = DPI)
    print(make_frame(all_t[k])); dev.off()
    file.exists(f)
  }, mc.cores = ncores, mc.preschedule = TRUE)
  ok <- sum(vapply(res, isTRUE, logical(1)))
  cat(sprintf("rendered %d/%d frames in %.0fs (%.2fs/frame)\n", ok, length(all_t),
              as.numeric(difftime(Sys.time(), t0, units = "secs")),
              as.numeric(difftime(Sys.time(), t0, units = "secs")) / length(all_t)))
  if (ok < length(all_t)) stop("Some frames failed to render; not encoding.")
  status <- system2(ffmpeg_bin, c("-y", "-framerate", FPS, "-i",
    file.path(fdir, "frame_%05d.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p",
    "-movflags", "+faststart", out_mp4))
  cat(sprintf("ffmpeg exit %d -> %s\n", status, out_mp4))
}

render_clip <- function(from, to, out = "figures/preview_clip.mp4", fps = 12L) {
  ts <- seq.POSIXt(as.POSIXct(from, tz = "UTC"), as.POSIXct(to, tz = "UTC"),
                   by = sprintf("%d min", FRAME_STEP_MIN))
  cdir <- "figures/anim_frames/clip"; unlink(cdir, recursive = TRUE); dir.create(cdir, recursive = TRUE)
  for (k in seq_along(ts)) {
    f <- sprintf("%s/c_%04d.png", cdir, k)
    agg_png(f, width = W_IN, height = H_IN, units = "in", res = DPI)
    print(make_frame(ts[k])); dev.off()
  }
  system2(ffmpeg_bin, c("-y", "-framerate", fps, "-i", file.path(cdir, "c_%04d.png"),
    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-movflags", "+faststart", out))
  cat(sprintf("wrote %s (%d frames)\n", out, length(ts)))
}

if (MODE == "design") render_design()
if (MODE == "clip")   render_clip(clip_from, clip_to)
if (MODE == "full")   render_full()
invisible(lapply(ncs, nc_close))
