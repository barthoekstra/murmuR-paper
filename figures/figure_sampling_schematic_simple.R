# Simplified three-panel "dynamic sampling" schematic figure.
#
# Left panel  -- North-to-south active night: simulated migrants in flight;
#                en-route sampling reads a NARROW, source-directed footprint
#                of ERA5 cells along inbound flight paths.
# Middle panel -- Suppressed night: no flying simulants, so the predictor
#                falls back to grounded migrants across the catchment, giving
#                a BROAD footprint.
# Right panel -- ENE-to-WSW active night, close in time to the first two;
#                sampling reads a NARROW, source-directed footprint of ERA5
#                cells along the inbound flight paths converging on the radar.
#
# Dynamically sampled ERA5 cells share one pale-yellow highlight, while the
# always-sampled radar core remains blue. This removes the temporal colour
# gradient while retaining the distinction between the two sampling areas.
#
# The script only READS murmuR-paper input data; nothing is modified.
# Set SAMPLING_SCHEMATIC_CHRONOLOGICAL=true to select and chronologically
# arrange a consecutive three-night no-migration / north-origin / east-origin
# window. The default retains the original selection and panel order.

suppressPackageStartupMessages({
  library(murmuR)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(sf)
  library(stars)
  library(ggplot2)
  library(patchwork)
  library(rnaturalearth)
})

options(murmuR.data_path = Sys.getenv(
  "SAMPLING_SCHEMATIC_DATA_PATH",
  unset = "/data/birdcloudstorage-tvm/ibm-ml-refactor/data/"
))

# ---------------------------------------------------------------------------
# Parameters (the brief's recommended values; printed below)
# ---------------------------------------------------------------------------
radar         <- "nlhrw"          # Herwijnen, inland
season        <- "autumn"
year          <- as.integer(Sys.getenv(
  "SAMPLING_SCHEMATIC_YEAR", unset = "2018"
))                                  # 2018 is retained as the default
                                  # in the local data tree and matches the
                                  # existing VPTS example figure narrative.
field_label   <- "850 hPa wind speed (m/s)"
flight_plevel <- ibm_pressurelevel       # 850 hPa -- the IBM's assumed flight level
chronological <- identical(
  tolower(Sys.getenv("SAMPLING_SCHEMATIC_CHRONOLOGICAL", unset = "false")),
  "true"
)

# Nocturnal hours we treat as "the night" (solar elevation < -6 deg).
night_solar_thresh <- -0.105      # radians (nautical twilight)

# Geometry constants (PULLED from murmuR::constants -- never hardcoded).
radar_lon   <- radars_named[[radar]][["lon"]]
radar_lat   <- radars_named[[radar]][["lat"]]
inner_dist  <- radar_dist     # 100 km radar buffer: always-include catchment
                              # core AND the polygon `annotate_radar_vpts()`
                              # uses to sample ERA5 vertical profiles every
                              # hour, regardless of bird state.
wedge_dist  <- if (season == "autumn") 400000 else 400000 # Herwijnen: 400 km
wedge_bear  <- if (season == "autumn") 40 else 220
wedge_width <- 65

# Discover which ERA5 pressure-level months we have on disk for this year.
# The night picks below are restricted to these months so the shared field
# underlay (850 hPa wind) is always available for both panels.
era5_dir         <- murmuR_path(sprintf("ecmwf-era5-pressurelevels/%d", year))
era5_files_all   <- list.files(era5_dir, pattern = "^\\d{6}\\.nc$",
                               full.names = TRUE)
era5_months      <- as.integer(substr(basename(era5_files_all), 5, 6))
season_avail_mons <- sort(intersect(era5_months, months_seasons[[season]]))
if (length(season_avail_mons) == 0L) {
  stop("No ERA5 pressure-level files for ", season, " ", year,
       " in ", era5_dir)
}

cat(sprintf(
  "Parameters: radar=%s | season=%s | year=%d | field=%s\n",
  radar, season, year, field_label
))

# ---------------------------------------------------------------------------
# Load tracks and identify nights
# ---------------------------------------------------------------------------
birds  <- simulated_birds(radar, year, season)
phenol <- load_phenology(radar, year, season)

# Build a tidx <-> datetime mapping for the season.
startmonth <- min(months_seasons[[season]])
endmonth   <- max(months_seasons[[season]])
dts_all <- seq.POSIXt(
  from = as.POSIXct(sprintf("%d-%02d-01 00:00:00", year, startmonth), tz = "UTC"),
  to   = as.POSIXct(sprintf("%d-%02d-01 00:00:00", year, endmonth + 1), tz = "UTC") -
           lubridate::hours(1),
  by   = "1 hour"
)
tidx_map <- data.table(tidx = seq_along(dts_all), datetime = dts_all)
setkey(tidx_map, tidx)
birds <- merge(birds, tidx_map, by = "tidx")

# Tag each tidx with a night_num using solar elevation.
phenol_n <- copy(phenol)
phenol_n[, is_night := solar_elev < night_solar_thresh]
setorder(phenol_n, tidx)
phenol_n[, transition := is_night != shift(is_night, fill = is_night[1])]
phenol_n[, grp := cumsum(transition)]
phenol_n[is_night == TRUE, night_num := .GRP, by = grp]
birds <- merge(birds, phenol_n[, .(tidx, night_num)], by = "tidx", all.x = TRUE)

# Per-night counts -- enable the north / suppressed / east-ENE choices.
night_summary <- birds[!is.na(night_num), .(
  date       = as.Date(min(datetime)),
  n_flying   = sum(state == 1L & radar == TRUE & migration == TRUE),
  n_grounded = sum(state == 0L),
  median_flight_direction = {
    dirs <- direction[state == 1L & radar == TRUE & migration == TRUE]
    if (length(dirs)) median(dirs, na.rm = TRUE) else NA_real_
  },
  prop_north_origin = {
    dirs <- direction[state == 1L & radar == TRUE & migration == TRUE]
    dirs <- dirs[is.finite(dirs)]
    if (length(dirs)) mean(dirs >= 135 & dirs < 225) else NA_real_
  },
  prop_east_origin = {
    dirs <- direction[state == 1L & radar == TRUE & migration == TRUE]
    dirs <- dirs[is.finite(dirs)]
    if (length(dirs)) mean(dirs >= 225 & dirs <= 315) else NA_real_
  }
), by = night_num][order(night_num)]

# Restrict choices to nights that fall in months with available ERA5 data
# (so the shared field underlay works for both panels).
night_summary[, mon := lubridate::month(date)]
avail_nights <- night_summary[mon %in% season_avail_mons]

if (chronological) {
  # Classify substantial flights by where they arrive from. Heading sectors
  # are deliberately broad; a night must have a strict majority in one sector.
  # Among all qualifying three-night windows, use the one maximizing the
  # weaker of the two migration nights.
  peak_flying <- max(avail_nights$n_flying)
  avail_nights[, night_class := fifelse(
    n_flying == 0L & n_grounded > 0L, "none",
    fifelse(
      n_flying >= 0.25 * peak_flying & prop_north_origin > 0.5, "north",
      fifelse(
        n_flying >= 0.25 * peak_flying & prop_east_origin > 0.5,
        "east", "other"
      )
    )
  )]

  requested_dates <- c(
    none = Sys.getenv("SAMPLING_SCHEMATIC_NONE_DATE"),
    north = Sys.getenv("SAMPLING_SCHEMATIC_NORTH_DATE"),
    east = Sys.getenv("SAMPLING_SCHEMATIC_EAST_DATE")
  )
  if (all(nzchar(requested_dates))) {
    requested_dates <- as.Date(requested_dates)
    selected <- avail_nights[match(requested_dates, date)]
    if (nrow(selected) != 3L || anyNA(selected$night_num)) {
      stop("One or more requested chronological dates are unavailable")
    }
    expected_classes <- c("none", "north", "east")
    if (!identical(selected$night_class, expected_classes)) {
      stop(
        "Requested dates classify as ",
        paste(selected$night_class, collapse = ", "),
        "; expected none, north, east"
      )
    }
    suppressed_night <- selected[1, night_num]
    active_night <- selected[2, night_num]
    directional_night <- selected[3, night_num]
  } else {
    chronological_windows <- rbindlist(lapply(
      seq_len(max(nrow(avail_nights) - 2L, 0L)),
      function(i) {
        z <- avail_nights[i:(i + 2L)]
        if (!all(diff(as.integer(z$date)) == 1L) ||
            !setequal(z$night_class, c("none", "north", "east"))) {
          return(NULL)
        }
        data.table(
          start_date = z$date[1],
          score = min(z[night_class %in% c("north", "east"), n_flying]),
          north_night = z[night_class == "north", night_num],
          suppressed_night = z[night_class == "none", night_num],
          directional_night = z[night_class == "east", night_num]
        )
      }
    ), fill = TRUE)
    if (nrow(chronological_windows) == 0L) {
      stop("No qualifying consecutive three-night window found")
    }
    setorder(chronological_windows, -score, start_date)
    active_night <- chronological_windows[1, north_night]
    suppressed_night <- chronological_windows[1, suppressed_night]
    directional_night <- chronological_windows[1, directional_night]
  }
  active_mon <- night_summary[night_num == active_night, mon]
} else {
  # Active night: max number of flying-over-radar bird-hours.
  active_night <- avail_nights[which.max(n_flying), night_num]
  active_mon   <- night_summary[night_num == active_night, mon]

  # Suppressed night: nearest night to active_night with n_flying == 0 and
  # n_grounded > 0, IN THE SAME ERA5 MONTH for a fair shared field.
  suppressed_candidates <- avail_nights[
    n_flying == 0L & n_grounded > 0L & mon == active_mon
  ]
  suppressed_night <- suppressed_candidates[
    which.min(abs(night_num - active_night)), night_num
  ]

  # Directional active night: a substantial movement heading between southwest
  # and west (225--270 degrees), i.e. arriving from the east/northeast.
  directional_candidates <- avail_nights[
    n_flying >= 0.25 * max(avail_nights$n_flying) &
      mon == active_mon &
      median_flight_direction >= 225 & median_flight_direction <= 270
  ]
  if (nrow(directional_candidates) == 0L) {
    stop("No substantial east/ENE-origin migration night found in the same ERA5 month")
  }
  directional_candidates[, night_gap := abs(night_num - active_night)]
  setorder(directional_candidates, night_gap, -n_flying)
  directional_night <- directional_candidates[1, night_num]
}

era5_file <- file.path(
  era5_dir, sprintf("%d%02d.nc", year, active_mon)
)

cat(sprintf(
  paste0("Chosen nights -- north-origin: night_num=%d ",
         "(%s, n_flying=%d, heading=%.1f deg) | suppressed: night_num=%d ",
         "(%s, n_flying=%d, n_grounded=%d) | east/ENE-origin: night_num=%d ",
         "(%s, n_flying=%d, heading=%.1f deg)\n"),
  active_night,
  format(night_summary[night_num == active_night, date]),
  night_summary[night_num == active_night, n_flying],
  night_summary[night_num == active_night, median_flight_direction],
  suppressed_night,
  format(night_summary[night_num == suppressed_night, date]),
  night_summary[night_num == suppressed_night, n_flying],
  night_summary[night_num == suppressed_night, n_grounded],
  directional_night,
  format(night_summary[night_num == directional_night, date]),
  night_summary[night_num == directional_night, n_flying],
  night_summary[night_num == directional_night, median_flight_direction]
))

active_tidxs     <- birds[night_num == active_night,     unique(tidx)]
suppressed_tidxs <- birds[night_num == suppressed_night, unique(tidx)]
directional_tidxs_all <- birds[night_num == directional_night, unique(tidx)]
active_times     <- tidx_map[tidx %in% active_tidxs,     datetime]
suppressed_times <- tidx_map[tidx %in% suppressed_tidxs, datetime]
directional_times_all <- tidx_map[
  tidx %in% directional_tidxs_all, datetime
]

# Autumn nights lengthen across September. Use the same UTC hours as the
# original active example in the later directional panel, so differences in
# tracks and footprints cannot be attributed to a longer plotted window.
reference_clock_hours <- format(active_times, "%H:%M")
directional_times <- directional_times_all[
  format(directional_times_all, "%H:%M") %in% reference_clock_hours
]
directional_tidxs <- tidx_map[datetime %in% directional_times, tidx]

cat(sprintf(
  "Plot windows -- north-origin: %s | suppressed: %s | directional: %s\n",
  paste(format(range(active_times), "%H:%M"), collapse = "--"),
  paste(format(range(suppressed_times), "%H:%M"), collapse = "--"),
  paste(format(range(directional_times), "%H:%M"), collapse = "--")
))

# ---------------------------------------------------------------------------
# Geometry: radar marker, 50 km always-include circle, catchment wedge
# ---------------------------------------------------------------------------
radar_pt <- st_sfc(st_point(c(radar_lon, radar_lat)), crs = 4326)

inner_buf <- radar_pt |>
  st_transform(3857) |>
  st_buffer(inner_dist) |>
  st_transform(4326)

wedge <- buffer_wedge(
  st_as_sf(data.frame(geometry = radar_pt)),
  radius       = wedge_dist,
  degree       = wedge_bear,
  degree_width = wedge_width
)

# Panel extent: catchment wedge bbox + small margin. Slightly larger south
# margin so the bottom-left counts have breathing room from the 50 km radar
# circle.
panel_bbox <- st_bbox(wedge)
margin_deg <- 0.5
south_extra <- 1.0
panel_bbox <- c(
  xmin = unname(panel_bbox["xmin"]) - margin_deg,
  ymin = unname(panel_bbox["ymin"]) - margin_deg - south_extra,
  xmax = unname(panel_bbox["xmax"]) + margin_deg,
  ymax = unname(panel_bbox["ymax"]) + margin_deg
)

# ---------------------------------------------------------------------------
# Load ERA5 field (10 m wind speed) and crop to panel extent
# ---------------------------------------------------------------------------
read_era5_wind <- function(file, times, bbox, plevel) {
  uv <- read_ncdf(
    file,
    var = c("u", "v"),
    proxy = FALSE
  )
  # Standardise dimension names: time and level.
  d <- names(st_dimensions(uv))
  if ("valid_time"     %in% d) uv <- st_set_dimensions(uv, "valid_time",
                                                       names = "time")
  if ("pressure_level" %in% d) uv <- st_set_dimensions(uv, "pressure_level",
                                                       names = "level")
  st_crs(uv) <- 4326
  crop_bbox <- st_bbox(
    setNames(c(bbox[["xmin"]], bbox[["ymin"]],
               bbox[["xmax"]], bbox[["ymax"]]),
             c("xmin", "ymin", "xmax", "ymax")),
    crs = 4326
  )
  uv <- st_crop(uv, crop_bbox)
  lv <- st_get_dimension_values(uv, "level")
  lvi <- which(lv == plevel)
  if (length(lvi) == 0L) stop("Pressure level ", plevel,
                              " not in ERA5 file ", file)
  uv <- uv[, , , lvi, ]                                  # subset on level
  tv <- st_get_dimension_values(uv, "time")
  keep <- which(tv %in% times)
  if (length(keep) == 0L) {
    stop("None of the requested datetimes are in ERA5 file ", file)
  }
  uv <- uv[, , , , keep]                                 # subset on time
  # Mean over the night (drops time dimension); return night-mean u and v.
  u_bar <- st_apply(uv["u"], c("longitude", "latitude"),
                    function(x) mean(x, na.rm = TRUE))
  v_bar <- st_apply(uv["v"], c("longitude", "latitude"),
                    function(x) mean(x, na.rm = TRUE))
  list(u = u_bar, v = v_bar)
}

wind_active     <- read_era5_wind(era5_file, active_times,     panel_bbox,
                                  flight_plevel)
wind_suppressed <- read_era5_wind(era5_file, suppressed_times, panel_bbox,
                                  flight_plevel)
wind_directional <- read_era5_wind(era5_file, directional_times, panel_bbox,
                                   flight_plevel)

# Build a long-format data frame of (lon, lat, u, v, ws) per cell, with the
# vector decimated to every Nth cell so the arrow field stays legible.
arrow_stride <- 3L  # plot every 3rd cell in each axis (~0.75 deg spacing)

vector_df <- function(uv, label, stride) {
  du <- as.data.frame(uv$u, xy = TRUE)
  dv <- as.data.frame(uv$v, xy = TRUE)
  names(du)[1:2] <- c("lon", "lat")
  names(dv)[1:2] <- c("lon", "lat")
  d <- merge(du, dv, by = c("lon", "lat"))
  names(d)[3:4] <- c("u", "v")
  if (inherits(d$u, "units")) d$u <- units::drop_units(d$u)
  if (inherits(d$v, "units")) d$v <- units::drop_units(d$v)
  d$ws <- sqrt(d$u ^ 2 + d$v ^ 2)
  lons <- sort(unique(d$lon))
  lats <- sort(unique(d$lat))
  keep_lon <- lons[seq(1, length(lons), by = stride)]
  keep_lat <- lats[seq(1, length(lats), by = stride)]
  d <- d[d$lon %in% keep_lon & d$lat %in% keep_lat, ]
  d$panel <- label
  d
}

vd_active     <- vector_df(wind_active,     "active",     arrow_stride)
vd_suppressed <- vector_df(wind_suppressed, "suppressed", arrow_stride)
vd_directional <- vector_df(wind_directional, "directional", arrow_stride)

field_limits <- range(c(vd_active$ws, vd_suppressed$ws, vd_directional$ws),
                      na.rm = TRUE, finite = TRUE)
# Three round breaks inside the observed range; expand limits to wrap them
# so the colour bar starts and ends on the labelled values.
field_breaks <- pretty(field_limits, n = 3)
field_breaks <- field_breaks[field_breaks >= field_limits[1] - 0.5 &
                             field_breaks <= field_limits[2] + 0.5]
if (length(field_breaks) > 3L) {
  field_breaks <- field_breaks[round(seq(1, length(field_breaks),
                                         length.out = 3))]
}
field_limits <- range(c(field_limits, field_breaks))
cat(sprintf("Shared %d hPa wind-speed range: [%.2f, %.2f] m/s | breaks: %s\n",
            flight_plevel, field_limits[1], field_limits[2],
            paste(field_breaks, collapse = ", ")))

# Convert each (u, v) vector to a geographic arrow whose visual length scales
# with wind speed but is bounded to roughly one ERA5 cell at the max so the
# field stays uncluttered.
arrow_scale <- 0.05  # degrees per (m/s); 15 m/s -> 0.75 deg

vd_active$lon_end     <- vd_active$lon     + vd_active$u     * arrow_scale
vd_active$lat_end     <- vd_active$lat     + vd_active$v     * arrow_scale
vd_suppressed$lon_end <- vd_suppressed$lon + vd_suppressed$u * arrow_scale
vd_suppressed$lat_end <- vd_suppressed$lat + vd_suppressed$v * arrow_scale
vd_directional$lon_end <- vd_directional$lon +
  vd_directional$u * arrow_scale
vd_directional$lat_end <- vd_directional$lat +
  vd_directional$v * arrow_scale

# ---------------------------------------------------------------------------
# Sampled footprint (the differentiator) -- snap bird positions to ERA5 cells
# ---------------------------------------------------------------------------
lon_centers <- sort(st_get_dimension_values(wind_active$u, "longitude"))
lat_centers <- sort(st_get_dimension_values(wind_active$u, "latitude"))
cell_dlon   <- diff(lon_centers)[1]
cell_dlat   <- diff(lat_centers)[1]

snap_to_cells <- function(positions) {
  if (nrow(positions) == 0L) {
    return(data.table(lon_c = numeric(0), lat_c = numeric(0)))
  }
  ix <- findInterval(positions$x, lon_centers - cell_dlon / 2,
                     all.inside = FALSE)
  iy <- findInterval(positions$y, lat_centers - cell_dlat / 2,
                     all.inside = FALSE)
  ix[ix < 1 | ix > length(lon_centers)] <- NA_integer_
  iy[iy < 1 | iy > length(lat_centers)] <- NA_integer_
  cells <- data.table(lon_c = lon_centers[ix], lat_c = lat_centers[iy])
  unique(cells[!is.na(lon_c) & !is.na(lat_c)])
}

# Active footprints: en-route cells along the inbound flights of birds that
# pass over the radar during each active night. Mirrors
# calculate_flight_conditions(): all (tidx_i, bird) rows of contributing
# birds for tidx_i <= max(night_tidxs), excluding each bird's first row
# (the stopover/takeoff position).
flight_sample_for_night <- function(night, night_tidxs) {
  ids <- birds[
    night_num == night & tidx %in% night_tidxs & state == 1L &
      radar == TRUE & migration == TRUE,
    unique(bird)
  ]
  positions <- birds[
    bird %in% ids & tidx <= max(night_tidxs)
  ][order(bird, tidx)][, if (.N > 1L) .SD[2:.N], by = bird]
  # Restrict to en-route (flight) positions, not stopover hops along the path.
  positions <- positions[state == 1L]
  list(ids = ids, positions = positions)
}

active_sample <- flight_sample_for_night(active_night, active_tidxs)
directional_sample <- flight_sample_for_night(
  directional_night, directional_tidxs
)
flying_ids <- active_sample$ids
directional_flying_ids <- directional_sample$ids
flight_positions <- active_sample$positions
directional_flight_positions <- directional_sample$positions

active_cells <- snap_to_cells(flight_positions)
active_cells$panel <- "active"
directional_cells <- snap_to_cells(directional_flight_positions)
directional_cells$panel <- "directional"

# Suppressed footprint: grounded non-migrating birds across the catchment
# during the suppressed night (the "fallback" sample).
grounded_positions <- birds[
  night_num == suppressed_night & state == 0L & migration == FALSE
]
suppressed_cells <- snap_to_cells(grounded_positions)
suppressed_cells$panel <- "suppressed"

# Footprint cells -> sf polygons (one 0.25 deg cell per row).
cells_to_polygons <- function(cells) {
  if (nrow(cells) == 0L) return(st_sf(geometry = st_sfc(crs = 4326)))
  pl <- lapply(seq_len(nrow(cells)), function(i) {
    cx <- cells$lon_c[i]; cy <- cells$lat_c[i]
    st_polygon(list(rbind(
      c(cx - cell_dlon / 2, cy - cell_dlat / 2),
      c(cx + cell_dlon / 2, cy - cell_dlat / 2),
      c(cx + cell_dlon / 2, cy + cell_dlat / 2),
      c(cx - cell_dlon / 2, cy + cell_dlat / 2),
      c(cx - cell_dlon / 2, cy - cell_dlat / 2)
    )))
  })
  st_sf(geometry = st_sfc(pl, crs = 4326))
}
active_cells_sf     <- cells_to_polygons(active_cells)
suppressed_cells_sf <- cells_to_polygons(suppressed_cells)
directional_cells_sf <- cells_to_polygons(directional_cells)

# Always-sampled ERA5 cells: every cell that intersects the 100 km radar
# buffer. These are read every hour by annotate_radar_vpts() and contribute
# the rdr_* radar-profile predictors regardless of bird state. We subtract
# them from each regime's dynamic footprint so each yellow cell is a cell
# the dynamic sampler reaches *in addition to* the always-sampled core.
buf_bbox <- st_bbox(inner_buf)
candidate_cells <- as.data.table(expand.grid(
  lon_c = lon_centers, lat_c = lat_centers
))
candidate_cells <- candidate_cells[
  lon_c >= buf_bbox[["xmin"]] - cell_dlon &
  lon_c <= buf_bbox[["xmax"]] + cell_dlon &
  lat_c >= buf_bbox[["ymin"]] - cell_dlat &
  lat_c <= buf_bbox[["ymax"]] + cell_dlat
]
candidate_polys <- cells_to_polygons(candidate_cells)
hits <- lengths(st_intersects(candidate_polys, inner_buf)) > 0L
always_cells    <- candidate_cells[hits]
always_cells_sf <- candidate_polys[hits, ]

always_keys <- paste(always_cells$lon_c, always_cells$lat_c)
strip_always <- function(cells) {
  cells <- copy(cells)
  cells[, key := paste(lon_c, lat_c)]
  cells[!key %in% always_keys, .(lon_c, lat_c)]
}
active_dyn_cells     <- strip_always(active_cells)
suppressed_dyn_cells <- strip_always(suppressed_cells)
directional_dyn_cells <- strip_always(directional_cells)
active_dyn_sf        <- cells_to_polygons(active_dyn_cells)
suppressed_dyn_sf    <- cells_to_polygons(suppressed_dyn_cells)
directional_dyn_sf   <- cells_to_polygons(directional_dyn_cells)

cat(sprintf("Footprint cells -- always (radar core): %d\n", nrow(always_cells)))
cat(sprintf("  active   total %d = %d always + %d dynamic\n",
            nrow(active_cells), nrow(always_cells),
            nrow(active_dyn_cells)))
cat(sprintf("  suppressed total %d = %d always + %d dynamic\n",
            nrow(suppressed_cells), nrow(always_cells),
            nrow(suppressed_dyn_cells)))
cat(sprintf("  east/ENE total %d = %d always + %d dynamic\n",
            nrow(directional_cells), nrow(always_cells),
            nrow(directional_dyn_cells)))

# ---------------------------------------------------------------------------
# Bird positions / tracks for visual reference
# ---------------------------------------------------------------------------
# Active panel: inbound flight track segments, one polyline per flying bird,
# trimmed to the active night's tidxs.
active_track_lines <- flight_positions[tidx %in% active_tidxs][order(bird, tidx)]
active_track_sf <- active_track_lines |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  group_by(bird) |>
  summarise(do_union = FALSE) |>
  st_cast("LINESTRING")

directional_track_lines <- directional_flight_positions[
  tidx %in% directional_tidxs
][order(bird, tidx)]
directional_track_sf <- directional_track_lines |>
  st_as_sf(coords = c("x", "y"), crs = 4326) |>
  group_by(bird) |>
  summarise(do_union = FALSE) |>
  st_cast("LINESTRING")

# Suppressed panel: grounded bird positions (drop duplicates per bird).
suppressed_pts_sf <- unique(grounded_positions[, .(bird, x, y)]) |>
  st_as_sf(coords = c("x", "y"), crs = 4326)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")

panel_xlim <- c(panel_bbox["xmin"], panel_bbox["xmax"])
panel_ylim <- c(panel_bbox["ymin"], panel_bbox["ymax"])

# Simple geographic scale bar drawn directly with annotate(): horizontal
# line of `length_km`, end ticks, "<n> km" label. Positioned in the
# bottom-right of the panel via fractional coordinates.
add_scalebar <- function(p, bbox, length_km = 100,
                         x_right_frac = 0.97, y_frac = 0.04) {
  lat_mid   <- mean(c(bbox[["ymin"]], bbox[["ymax"]]))
  m_per_deg <- 111320 * cos(lat_mid * pi / 180)
  bar_dlon  <- length_km * 1000 / m_per_deg
  x_right   <- bbox[["xmin"]] + x_right_frac *
                 (bbox[["xmax"]] - bbox[["xmin"]])
  x_left    <- x_right - bar_dlon
  y_bar     <- bbox[["ymin"]] + y_frac *
                 (bbox[["ymax"]] - bbox[["ymin"]])
  tick_h    <- 0.025 * (bbox[["ymax"]] - bbox[["ymin"]])

  p +
    annotate("segment", x = x_left, xend = x_right,
             y = y_bar, yend = y_bar,
             colour = "black", linewidth = 0.6) +
    annotate("segment", x = x_left, xend = x_left,
             y = y_bar - tick_h, yend = y_bar + tick_h,
             colour = "black", linewidth = 0.6) +
    annotate("segment", x = x_right, xend = x_right,
             y = y_bar - tick_h, yend = y_bar + tick_h,
             colour = "black", linewidth = 0.6) +
    annotate("text", x = (x_left + x_right) / 2,
             y = y_bar + tick_h * 1.6,
             label = paste0(length_km, " km"),
             size = 2.8, vjust = 0, hjust = 0.5,
             colour = "grey15", fontface = "bold")
}

make_panel <- function(vec_data, dyn_cells_sf, always_cells_sf,
                       overlay_sf, overlay_type,
                       title, subtitle, counts) {
  p <- ggplot() +
    geom_sf(data = world, fill = "grey97", colour = "grey55",
            linewidth = 0.25) +
    geom_sf(data = wedge,     fill = NA, colour = "black",
            linewidth = 0.6, linetype = "dashed") +
    geom_sf(data = always_cells_sf,
            aes(fill = "Always sampled (radar core)"),
            colour = "#3a7a8c", alpha = 0.55, linewidth = 0.15) +
    geom_sf(data = dyn_cells_sf,
            aes(fill = "ERA5 cells used this night"),
            colour = "#b58900", alpha = 0.55, linewidth = 0.15) +
    geom_sf(data = inner_buf, fill = NA, colour = "black", linewidth = 0.4)

  if (overlay_type == "lines" && nrow(overlay_sf) > 0L) {
    p <- p + geom_sf(data = overlay_sf,
                     colour = "#111111", linewidth = 0.15, alpha = 0.6)
  } else if (overlay_type == "points" && nrow(overlay_sf) > 0L) {
    p <- p + geom_sf(data = overlay_sf,
                     colour = "#111111", size = 0.45, alpha = 0.7,
                     shape = 16)
  }

  p <- p +
    geom_segment(
      data = vec_data,
      aes(x = lon, y = lat, xend = lon_end, yend = lat_end, colour = ws),
      arrow = arrow(length = unit(0.05, "inches"), type = "closed"),
      linewidth = 0.45
    ) +
    geom_sf(data = radar_pt, colour = "black", fill = "white",
            shape = 21, size = 2.4, stroke = 0.7)

  # Catchy counts clustered in the bottom-left corner with a "Sampled:"
  # header above. Numbers right-aligned at a shared column; descriptions
  # left-aligned at a shared column to the right, on the same baseline.
  # Header sits above the cluster, aligned with the left edge of the
  # widest number. Added LAST so all of this sits on top of the wind
  # arrows and tracks.
  n_c         <- nrow(counts)
  vstep       <- 0.48
  y_bottom    <- panel_bbox[["ymin"]] + 0.35
  row_y       <- y_bottom + (n_c - seq_len(n_c)) * vstep
  x_num_right <- panel_bbox[["xmin"]] + 1.30
  x_desc_left <- x_num_right + 0.10
  y_header    <- max(row_y) + 0.50
  x_header    <- panel_bbox[["xmin"]] + 0.22
  x_box_left  <- panel_bbox[["xmin"]] + 0.12
  box_width   <- if (any(grepl("^grounded", counts$small))) 5.65 else 4.95
  x_box_right <- panel_bbox[["xmin"]] + box_width
  y_box_bottom <- min(row_y) - 0.30
  y_box_top    <- y_header + 0.30

  p <- p +
    annotate(
      "rect", xmin = x_box_left, xmax = x_box_right,
      ymin = y_box_bottom, ymax = y_box_top,
      fill = "white", alpha = 0.90,
      colour = "grey55", linewidth = 0.25
    ) +
    annotate("text", x = x_header, y = y_header, label = "Sampled:",
             hjust = 0, vjust = 0.5, fontface = "bold",
             size = 3.4, colour = "grey15")
  for (i in seq_len(n_c)) {
    p <- p +
      annotate("text", x = x_num_right, y = row_y[i], label = counts$big[i],
               hjust = 1, vjust = 0.5, fontface = "bold",
               size = 6.0, colour = counts$big_col[i]) +
      annotate("text", x = x_desc_left, y = row_y[i], label = counts$small[i],
               hjust = 0, vjust = 0.5, fontface = "bold",
               size = 3.0, colour = counts$big_col[i])
  }

  p +
    scale_colour_viridis_c(
      option = "mako", direction = -1,
      name = field_label, limits = field_limits,
      breaks = field_breaks,
      labels = field_breaks,
      oob = scales::squish,
      guide = guide_colorbar(
        order = 1,
        direction = "horizontal",
        barwidth = unit(7.8, "lines"), barheight = unit(0.7, "lines"),
        title.position = "top", title.hjust = 0
      )
    ) +
    scale_fill_manual(
      name   = NULL,
      breaks = c("Always sampled (radar core)",
                 "ERA5 cells used this night"),
      values = c("Always sampled (radar core)" = "#b8d8e3",
                 "ERA5 cells used this night" = "#fff7c2"),
      guide = guide_legend(
        order = 2,
        override.aes = list(
          colour    = c("#3a7a8c", "#b58900"),
          linewidth = 0.3,
          alpha     = 0.55
        ),
        keywidth = unit(0.8, "lines"), keyheight = unit(0.8, "lines")
      )
    ) +
    coord_sf(xlim = panel_xlim, ylim = panel_ylim, expand = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid       = element_line(colour = "grey85", linewidth = 0.2),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border     = element_rect(fill = NA, colour = "grey40",
                                      linewidth = 0.4),
      plot.title       = element_text(face = "bold", size = 11.5),
      plot.subtitle    = element_text(colour = "grey30", size = 8.5,
                                      lineheight = 1.1),
      legend.title     = element_text(size = 8.5),
      legend.text      = element_text(size = 8),
      axis.text        = element_text(size = 7, colour = "grey40")
    )
}

active_dates <- range(active_times)
supp_dates   <- range(suppressed_times)
directional_dates <- range(directional_times)
fmt_night <- function(r) {
  same_month <- format(r[1], "%Y-%m") == format(r[2], "%Y-%m")
  if (same_month) {
    sprintf("%s-%s %s, %s-%s UTC",
            format(r[1], "%-d"), format(r[2], "%-d"),
            format(r[1], "%b %Y"),
            format(r[1], "%H:%M"), format(r[2], "%H:%M"))
  } else {
    sprintf("%s, %s - %s, %s UTC",
            format(r[1], "%-d %b %Y"), format(r[1], "%H:%M"),
            format(r[2], "%-d %b %Y"), format(r[2], "%H:%M"))
  }
}

# Cohorts are identified by their preferred migration direction (`pref_dir`).
# Show the number and share of sampled birds from each cohort rather than just
# the number of represented cohorts. This makes a token second cohort visually
# distinct from a genuinely mixed nightly sample.
cohort_map <- rbindlist(lapply(names(bird_routes[[season]]), function(x) {
  data.table(
    pref_dir = bird_routes[[season]][[x]]$endogenous_heading,
    cohort = sub("^cohort", "cohort ", x)
  )
}))

cohort_annotation_rows <- function(sampled_birds) {
  z <- unique(sampled_birds[!is.na(pref_dir), .(bird, pref_dir)])
  z <- z[, .(n_birds = uniqueN(bird)), by = pref_dir]
  z <- merge(z, cohort_map, by = "pref_dir", all.x = TRUE, sort = FALSE)
  z[, cohort_order := match(pref_dir, cohort_map$pref_dir)]
  setorder(z, cohort_order)
  z[, share := n_birds / sum(n_birds)]
  data.frame(
    big = format(z$n_birds, big.mark = ","),
    small = sprintf(
      "%s (%s°, %.0f%%)", z$cohort, format(z$pref_dir), 100 * z$share
    ),
    big_col = rep("grey10", nrow(z))
  )
}

active_cohort_rows <- cohort_annotation_rows(
  birds[bird %in% flying_ids, .(bird, pref_dir)]
)
suppressed_cohort_rows <- cohort_annotation_rows(
  grounded_positions[, .(bird, pref_dir)]
)
directional_cohort_rows <- cohort_annotation_rows(
  birds[bird %in% directional_flying_ids, .(bird, pref_dir)]
)

counts_active <- rbind(data.frame(
  big     = c(format(length(flying_ids), big.mark = ","),
              format(nrow(active_dyn_cells), big.mark = ",")),
  small   = c("simulated migrants",
              "ERA5 cells"),
  big_col = rep("grey10", 2)
), active_cohort_rows)[c(1, seq_len(nrow(active_cohort_rows)) + 2, 2), ]
counts_suppressed <- rbind(data.frame(
  big     = c("0",
              format(nrow(suppressed_pts_sf), big.mark = ","),
              format(nrow(suppressed_dyn_cells), big.mark = ",")),
  small   = c("simulated migrants",
              "grounded simulated migrants",
              "ERA5 cells"),
  big_col = rep("grey10", 3)
), suppressed_cohort_rows)[
  c(1:2, seq_len(nrow(suppressed_cohort_rows)) + 3, 3),
]
counts_directional <- rbind(data.frame(
  big     = c(format(length(directional_flying_ids), big.mark = ","),
              format(nrow(directional_dyn_cells), big.mark = ",")),
  small   = c("simulated migrants",
              "ERA5 cells"),
  big_col = rep("grey10", 2)
), directional_cohort_rows)[
  c(1, seq_len(nrow(directional_cohort_rows)) + 2, 2),
]

active_title <- if (chronological) {
  "Migration predominantly from the north"
} else {
  "Southwest-oriented cohort dominates the flight"
}
directional_title <- if (chronological) {
  "Migration predominantly from the east"
} else {
  "West-oriented cohort dominates a mixed flight"
}

p_active <- make_panel(
  vec_data        = vd_active,
  dyn_cells_sf    = active_dyn_sf,
  always_cells_sf = always_cells_sf,
  overlay_sf      = active_track_sf,
  overlay_type    = "lines",
  title           = active_title,
  subtitle        = fmt_night(active_dates),
  counts          = counts_active
) +
  theme(legend.position = "none")

p_suppressed <- make_panel(
  vec_data        = vd_suppressed,
  dyn_cells_sf    = suppressed_dyn_sf,
  always_cells_sf = always_cells_sf,
  overlay_sf      = suppressed_pts_sf,
  overlay_type    = "points",
  title           = "No migration: simulated birds grounded",
  subtitle        = fmt_night(supp_dates),
  counts          = counts_suppressed
) +
  theme(
    legend.position         = "inside",
    legend.position.inside  = c(0.98, 0.02),
    legend.justification    = c("right", "bottom"),
    legend.box              = "vertical",
    legend.spacing.y        = unit(0.1, "lines"),
    legend.background       = element_blank(),
    legend.box.background   = element_rect(fill = "white", colour = "grey70",
                                           linewidth = 0.2),
    legend.box.margin       = margin(3, 6, 3, 6),
    legend.box.just         = "left",
    legend.margin           = margin(0, 0, 0, 0)
  )

p_directional <- make_panel(
  vec_data        = vd_directional,
  dyn_cells_sf    = directional_dyn_sf,
  always_cells_sf = always_cells_sf,
  overlay_sf      = directional_track_sf,
  overlay_type    = "lines",
  title           = directional_title,
  subtitle        = fmt_night(directional_dates),
  counts          = counts_directional
) +
  theme(legend.position = "none")

if (chronological) {
  panels <- list(
    active = p_active,
    suppressed = p_suppressed,
    directional = p_directional
  )
  panel_nights <- c(
    active = active_night,
    suppressed = suppressed_night,
    directional = directional_night
  )
  fig <- wrap_plots(panels[order(panel_nights)], nrow = 1)
} else {
  fig <- p_active | p_suppressed | p_directional
}

# ---------------------------------------------------------------------------
# Acceptance checks (warn loudly if any fail)
# ---------------------------------------------------------------------------
check <- function(cond, msg) {
  if (!isTRUE(cond)) warning("ACCEPTANCE CHECK FAILED: ", msg, call. = FALSE)
}
check(nrow(active_cells) < nrow(suppressed_cells),
      "n_cells(active) should be << n_cells(suppressed)")
check(night_summary[night_num == active_night, n_flying] > 0,
      "t_active must have n_flying > 0")
check(night_summary[night_num == suppressed_night, n_flying] == 0 &&
      night_summary[night_num == suppressed_night, n_grounded] > 0,
      "t_suppressed must have n_flying == 0 and n_grounded > 0")
check(night_summary[night_num == directional_night, n_flying] > 0,
      "t_directional must have n_flying > 0")
if (!chronological) {
  check(night_summary[night_num == directional_night,
                      median_flight_direction] >= 225 &&
        night_summary[night_num == directional_night,
                      median_flight_direction] <= 270,
        "t_directional must have a southwest-to-west median heading")
}
check(abs(directional_night - active_night) <= 7,
      "directional example should be within one week of the active example")
if (chronological) {
  selected_dates <- sort(night_summary[
    night_num %in% c(active_night, suppressed_night, directional_night), date
  ])
  check(length(selected_dates) == 3L &&
        diff(range(as.integer(selected_dates))) <= 6L,
        "chronological panels must come from a window of at most seven nights")
  check(uniqueN(night_summary[
    night_num %in% c(active_night, suppressed_night, directional_night), mon
  ]) == 1L, "selected nights must use the same ERA5 month")
  check(night_summary[night_num == active_night, prop_north_origin] > 0.5,
        "north-origin night must have a majority of north-origin flight-hours")
  check(night_summary[night_num == directional_night, prop_east_origin] > 0.5,
        "east-origin night must have a majority of east-origin flight-hours")
}
check(length(directional_times) == length(active_times) &&
      identical(format(directional_times, "%H:%M"),
                format(active_times, "%H:%M")),
      "directional and original active panels must use identical UTC hours")
catchment_union <- st_union(wedge, inner_buf)
n_act_inside <- sum(lengths(st_intersects(active_dyn_sf, catchment_union)) > 0L)
n_sup_inside <- sum(lengths(st_intersects(suppressed_dyn_sf, catchment_union)) > 0L)
n_dir_inside <- sum(lengths(st_intersects(directional_dyn_sf,
                                          catchment_union)) > 0L)
cat(sprintf(paste0("Dynamic footprint cells inside catchment -- active: ",
                   "%d/%d (%.0f%%) | suppressed: %d/%d (%.0f%%) | ",
                   "east/ENE: %d/%d (%.0f%%)\n"),
            n_act_inside, nrow(active_dyn_sf),
            100 * n_act_inside / max(nrow(active_dyn_sf), 1),
            n_sup_inside, nrow(suppressed_dyn_sf),
            100 * n_sup_inside / max(nrow(suppressed_dyn_sf), 1),
            n_dir_inside, nrow(directional_dyn_sf),
            100 * n_dir_inside / max(nrow(directional_dyn_sf), 1)))

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------
out_dir <- Sys.getenv("SAMPLING_SCHEMATIC_OUT_DIR", unset = "figures")
out_suffix <- if (chronological) "_chronological" else ""
out_tag <- Sys.getenv("SAMPLING_SCHEMATIC_TAG")
if (nzchar(out_tag)) out_suffix <- paste0(out_suffix, "_", out_tag)
out_pdf <- file.path(
  out_dir, paste0("sampling_schematic_simple", out_suffix, ".pdf")
)
out_png <- file.path(
  out_dir, paste0("sampling_schematic_simple", out_suffix, ".png")
)

ggsave(out_pdf, fig, width = 18.25, height = 6.5, units = "in")
ggsave(out_png, fig, width = 18.25, height = 6.5, units = "in", dpi = 300)

cat(sprintf("Wrote %s and %s\n", out_pdf, out_png))
