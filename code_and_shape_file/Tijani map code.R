## For Tijani to plot map and seedz wards
## code I used to plot study area
library(rgdal)
library(ggplot2)
library(ggsn) # cant remember why I needed this package
library(ggmap)
library(sf)

## seedz.wards data file attached

# seedz.wards <- readOGR("/Users/gemmachaters/Dropbox/PHD/GIS data/KiliArushaManyaraGIS",
#                                       "seedz.wards")
seedz.wards <- readOGR("/Users/s2120283/OneDrive - University of Edinburgh/PhD projects/movementSimulation/code_and_shape_file", "seedz.wards")
## this is the code I used to get the google maps tile. I'm not sure if the key is unique, you might need
## to get another key, I followed instructions online, I had to change my settings withing my google account to allow
## the api to be accessed or something but it was fairly straightforward. think i used a youtube video or something!
ward.poly <- sf::st_read("2012 Wards Shapefiles/regions/Regions.shp") %>% 
  filter(Region_Nam == "Manyara"|
         Region_Nam == "Arusha"|
         Region_Nam == "Kilimanjaro")

region.cord <-
  data.frame(
    Region = c("Arusha", "Kilimanjaro", "Manyara"),
    long = c(36, 37.63, 36.83),
    lat = c(-3, -3.75, -4.5))


ggplot(data = region.cord) +
    geom_point(aes(x = long, y = lat, color = Region),
            show.legend = FALSE)
  
  
register_google(key = "AIzaSyBzvwlcWRB2SjtqFGZuOzcKSGAUzxbVwZo", write = TRUE)
TZmap <- ggmap::get_googlemap(c(35.5, -5.3), zoom = 6)
ggmap(TZmap, extent = "device")

# study.regions <- 
#   seedz.wards[seedz.wards$Region_Nam == "Manyara" | 
#                 seedz.wards$Region_Nam == "Arusha" | 
#                 seedz.wards$Region_Nam == "Kilimanjaro",]

study.regions <-
  subset(seedz.wards, Region_Nam == "Manyara" | Region_Nam == "Arusha" | Region_Nam == "Kilimanjaro")

stdy.reg <- st_read("2012 Wards Shapefiles/regions/Regions.shp")

stdy.reg <- 
  stdy.reg[stdy.reg$Region_Nam == "Manyara" | 
             stdy.reg$Region_Nam == "Arusha" | 
             stdy.reg$Region_Nam == "Kilimanjaro", ]


#pdf("sim.output/studyArea.pdf")
ggmap(TZmap, extent = "device") +
  geom_polygon(data=fortify(study.regions,),
               aes(x=long, y=lat, group=group), fill = "#471164FF") +
  geom_point(data = region.cord, aes(x = long, y = lat), 
             size = 1.5, shape = 21, fill = "lightgray", color = "black") +
  # geom_point(data = ward.poly, aes(x = long, y = lat, color = Region_Nam),
  #            show.legend = FALSE, size = 0.05) +
  annotate("text", x = 36, y = -3.1, label = "Arusha", size = 2, color = "white")+
  annotate("text", x = 37.53, y = -3.75, label = "Kilimanjaro", size = 2, color = "white", angle = -55)+
  annotate("text", x = 36.83, y = -4.6, label = "Manyara", size = 2, color = "white")+
  theme_void()
#dev.off()  

#+ ## unhash if you want a title
  #  ggtitle("Arusha and Manyara study regions in northern Tanzania\nhighlighted in purple") ## edit title if using one
  
  show_col(viridis_pal()(25)) ## the colour I chose to highlight the study regions was dark purple. feel free to change the fill = "#471164FF"
