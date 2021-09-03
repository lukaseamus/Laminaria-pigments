#################################################################
##### Project: Pigment composition of Laminaria detritus    #####
##### Script purpose: Plotting a Notheast Atlantic base map #####
##### Author: Luka Seamus Wright                            #####
#################################################################

#### Load data ####
require(rworldmap)
Europe <- getMap(resolution = "high")

#### Set theme ####
require(ggplot2)
mytheme <- theme(panel.background = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "sans"))
#### Plot ####
ggplot() + 
  coord_map(xlim = c(-30, 30), ylim = c(27, 69)) +
  geom_polygon(data = Europe, mapping = aes(x = long, y = lat, group = group),
               fill = "#5B5B5E", colour = "#ffffff", size = 0.1) +
  geom_segment(aes(x = -28.8, xend = -17, y = 30, yend = 30)) +
  geom_text(aes(x = -22.9, y = 29, label = "1000 km"), size = 5) +
  mytheme

#### Clean up ####
detach(package:rworldmap)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")
