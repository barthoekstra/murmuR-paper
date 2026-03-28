library(murmuR)
options(murmuR.data_path = "/data/birdcloudstorage-tvm/ibm-ml-refactor/data/")

library(dplyr)
library(purrr)
library(sf)
library(geosphere)        # for buffer_wedge_metric() internals

# --- params ---
wedge_center <- 40
wedge_width  <- 65

# helper: wrap a single POINT sfg into a 1-row sf for buffer_wedge_metric()
.point_sf <- function(pt, crs) {
  sf::st_as_sf(tibble::tibble(id = 1), geometry = sf::st_sfc(pt, crs = crs))
}

# --- your input ---
radars <- tibble::tibble(
  name   = c("Herwijnen", "Herwijnen", "Den Helder", "Den Helder"),
  lat    = c(51.8371, 51.8371, 52.9528, 52.9528),
  lon    = c(5.138,   5.138,   4.7906,  4.7906),
  season = c("autumn", "spring", "autumn", "spring"),
  dist   = c(400000,  400000,   650000,  400000)  # season-specific radius (m)
)

radars_sf <- radars %>%
  st_as_sf(coords = c("lon","lat"), crs = 4326, remove = FALSE)

radars_buf <- radars_sf %>%
  st_transform(3857) %>%
  mutate(geometry = st_buffer(geometry, 100000)) %>%  # per-row radius
  st_transform(4326)

radars_wedge <- radars_sf %>%
  mutate(
    bearing = if_else(season == "autumn", wedge_center, (wedge_center + 180) %% 360),
    wedge_geom = st_sfc(
      pmap(list(geometry, dist, bearing), ~{
        w <- buffer_wedge_metric(
          point        = .point_sf(..1, st_crs(radars_sf)),
          radius       = ..2,         # meters (season-specific)
          degree       = ..3,         # per-row bearing
          degree_width = wedge_width
        )
        st_geometry(w)[[1]]           # sfg POLYGON
      }),
      crs = 4326
    )
  ) %>%
  st_set_geometry("wedge_geom") %>%
  select(name, lon, lat, dist, season)

world <- rnaturalearth::ne_countries(scale = "large", returnclass = "sf")

ggplot() +
  geom_sf(data = world, fill = "white") +
  geom_sf(data = radars_wedge, aes(color = name, linetype = season), fill = NA, linewidth = 0.6) +
  geom_sf(data = radars_buf,   aes(color = name, linetype = season), fill = NA, linewidth = 0.5) +
  geom_sf(data = radars_sf,    aes(color = name),                     size = 1) +
  scale_color_discrete(name = "Radar") +
  scale_linetype_discrete(name = "Season", labels = c(autumn = "Autumn", spring = "Spring")) +
  xlim(-5, 20) + ylim(48, 60) +
  coord_sf() +
  theme_minimal(base_size = 14)

ggsave(murmuR_path("models/plots/catchment_areas.png"), width = 10, height = 7)