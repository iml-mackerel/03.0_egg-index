root<- ""
library(leaflet)
library(sf)
library(ggplot2)
library(webshot2)
lookup <- read.delim(paste0(root,"data/lookup_station_egg.txt")) # position and depth of stations.

# ---- Prepare station and polygon data ----
# 1. Convert your station to 'sf' 
lookup_sf <- st_as_sf(lookup, coords = c("longitude", "latitude"), crs = 4326, remove=F)

m<- leaflet(options = leafletOptions(
  dragging = FALSE,
  zoomControl = FALSE,
  scrollWheelZoom = FALSE,
  doubleClickZoom = FALSE,
  boxZoom = FALSE,
  keyboard = FALSE,
  minZoom = 7.5,  # Set the same value for min and max zoom!
  maxZoom = 7.5
))  %>%
  addProviderTiles(providers$Esri.OceanBasemap) %>%   
  
  addCircleMarkers(
    data = lookup_sf,
    lng = ~longitude,
    lat = ~latitude,
    radius = 10,
    fillColor =    "white",  
    fillOpacity = 1,
    stroke = TRUE,
    color =   "black",           # White border for added contrast
    weight = 5                   # Thicker border
  ) %>%
  
  # Add scale bar
  addScaleBar(
    position = "bottomleft",
    options = scaleBarOptions(metric = TRUE, imperial = TRUE)
  ) 

# 1. Save the map as HTML
htmlwidgets::saveWidget(m, file = paste0(root,"img/2024/area_of_interest.html"), selfcontained = TRUE)

# 2. Screenshot the HTML as PNG
webshot(
  paste0(root,"img/2024/area_of_interest.html"), 
  file = paste0(root,"img/2024/area_of_interest.png"), 
  vwidth = 1000/1.2,  
  vheight = 800/1.2,
  zoom = 4      
)

