# Map for Collection Sites of Chapter 3 Data

# load libraries
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

#---------
# Making three different maps, one for each lake and corresponding rivers, to
# plot sea lamprey collection locations.
#---------

# setting up map
world <- rnaturalearth::ne_countries(scale = "large",
                                     returnclass = "sf")
n_america <- world %>% 
  dplyr::filter(adm0_a3 %in% c("CAN", "USA"))

lakes <- rnaturalearth::ne_download(scale = 10, 
                                    type = 'lakes', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269)


# LAKE HURON
lake_huron <- data.frame(
  river = c("Carp River", "Cheboygan River"),
  latitude = c(46.0129,45.6561),
  longitude = c(-84.4133,-84.4653),
  phase = c("adult", "adult")
)
# plot
ggplot() +
  geom_sf(data = n_america,
          mapping = aes(geometry = geometry),
          color = "black",
          fill = "lightgray")  +
  geom_sf(data = lakes,
          mapping = aes(geometry = geometry),
          color = "black",
          fill = "lightblue") +
  geom_point(data = lake_huron, aes(x = longitude, y = latitude, color = phase), size = 3) +
  coord_sf(ylim = c(43, 46.5),
           xlim = c(-85, -81),
           expand = TRUE) +
  labs(x = "Longitude",
       y = "Latitude") 



# LAKE SUPERIOR
lake_superior <- data.frame(
  river = c("Misery", "Falls", "Firesteel" ),
  latitude = c(46.99722,46.75688,46.93466),
  longitude = c(-88.98111,-88.45708,-89.19653),
  phase = "adult"
)
# plot
ggplot() +
  geom_sf(data = n_america,
          mapping = aes(geometry = geometry),
          color = "black",
          fill = "lightgray")  +
  geom_sf(data = lakes,
          mapping = aes(geometry = geometry),
          color = "black",
          fill = "lightblue") +
  geom_point(data = lake_superior, aes(x = longitude, y = latitude, color = phase), size = 3) +
  coord_sf(ylim = c(46, 49),
           xlim = c(-92, -84),
           expand = TRUE) +
  labs(x = "Longitude",
       y = "Latitude") 


# LAKE CHAMPLAIN
lake_champlain <- data.frame(
  river = c("Mallets", "Sunderland", "Pond", "Mullen",
            "Great Chazy", "Winooski", "Morpion"),
  latitude = c(44.58050,44.50727,44.557001,44.1274,44.93237,44.5304,45.072),
  longitude = c(-73.1839,-73.2132,-73.19100,-73.4630,-73.38699,-73.2745,-73.0971),
  phase = "adult"
)
# plot
ggplot() +
  geom_sf(data = n_america,
          mapping = aes(geometry = geometry),
          color = "black",
          fill = "lightgray")  +
  geom_sf(data = lakes,
          mapping = aes(geometry = geometry),
          color = "black",
          fill = "lightblue") +
  geom_point(data = lake_champlain, aes(x = longitude, y = latitude, color = phase), size = 3) +
  coord_sf(ylim = c(43.7, 45.1),
           xlim = c(-73.8, -72.9),
           expand = TRUE) +
  labs(x = "Longitude",
       y = "Latitude") +
  scale_x_continuous(breaks = c(-73.8, -73.3, -72.9))



