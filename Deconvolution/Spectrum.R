##################################################################
##### Project: Pigment composition of Laminaria detritus     #####
##### Script purpose: Annotated spectrum to show methods     #####
##### Author: Luka Seamus Wright                             #####
##################################################################

#### 1.    Data preparation ####
#### 1.1   Load functions and data ####
source("~/Desktop/Plymouth University/Dissertation/Pigments/Deconvolution/pigment.function.R")
s <- read.csv("~/Desktop/Plymouth University/Dissertation/Pigments/Deconvolution/Dry.csv", header = T)
s <- subset(s, wavelength >= 400 & wavelength <= 700)[, -1]
w <- 400:700

#### 1.2   Calculate unweighted pigment spectra ####
uwps <- pigment.basis(w)
matplot(w, uwps, type = "l")

#### 1.3   Calculate pigment weights ####
weights <- pigment.fit(w, s)
pw <- weights$m[[1]]$x[8:15] # extract pigment weights for sample d1A

#### 1.4   Calculate weighted pigment spectra ####
wps <- uwps %*% diag(pw)
colnames(wps) <- colnames(uwps)
wps <- data.frame(wps, wavelength = w)

#### 1.5   Calculate summed fitted spectra ####
sps <- pigment.spectrum(weights)
sbs <- background.spectrum(weights)
ss <- fitted.spectrum(weights)

#### 1.6   Add raw data and fitted spectra to dataframe ####
wps$raw <- s$d1A # raw absorbance data for sample d1A
wps$pigments <- sps[,1] # summed pigment spectrum for sample d1A
wps$background <- sbs[,1] # summed background spectrum for sample d1A
wps$total <- ss[,1] # total summed spectrum (pigment + background) for sample d1A

#### 1.7   Remove unrepresented pigments ####
wps$bb.Car <- NULL
wps$Phe.a <- NULL
wps$Viola <- NULL

#### 1.8   Create dataframe for plotting ####
# summed background and pigment spectra are of little interest
# because they are only components of the total summed spectrum
# and therefore not included in this dataframe
df <- data.frame(wavelength = wps$wavelength,
                 absorbance = with(wps, c(Chl.a, Chl.c1, Chl.c2, Fuco, Zea)),
                 fitted.abs = wps$total,
                 raw.abs = wps$raw,
                 group = factor(rep(colnames(wps)[1:5], each = 301)))

# data are now ready to be plotted with ggplot2!

#### 2.    Data visualisation ####
require(ggplot2)
require(ggspectra)

#### 2.1   Customise theme ####
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .4, .2, .2),"cm"),
                 axis.line = element_line(),
                 axis.title = element_text(size = 15),
                 axis.text = element_text(size = 12, colour = "black"),
                 axis.ticks.length = unit(.25, "cm"),
                 axis.ticks = element_line(colour = "black"),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title = element_blank(),
                 text = element_text(family = "Helvetica Neue"))

#### 2.2   Plot ####
spec <- ggplot(data = df, aes(x = wavelength)) + 
  geom_point(aes(y = raw.abs), shape = 16, size = 2, alpha = 0.04) +
  geom_line(aes(y = fitted.abs)) +
  geom_line(aes(y = absorbance, colour = group, lty = group)) +
  scale_colour_manual(values = c("#0d98ba", "#50c878", "#50c878", "#ffa500", "#fedf00"),
                      labels = c(expression("Chlorophyll"*italic(" a")), expression("Chlorophyll"*italic(" c")[1]),
                                 expression("Chlorophyll"*italic(" c")[2]), "Fucoxanthin", "Zeaxanthin"),
                      guide = guide_legend()) +
  scale_linetype_manual(values = c(1, 1, 5, 1, 1),
                      labels = c(expression("Chlorophyll"*italic(" a")), expression("Chlorophyll"*italic(" c")[1]),
                                 expression("Chlorophyll"*italic(" c")[2]), "Fucoxanthin", "Zeaxanthin"),
                      guide = guide_legend()) +
  ylab("Absorbance (arbitrary units)") +
  xlab("Wavelength (nm)") +
  theme(legend.position = c(.885, .8)) +
  scale_x_continuous(expand = c(0,0), breaks = seq(400, 700, by = 50)) +
  scale_y_continuous(expand = c(0,0)) +
  wl_guide(ymax = 2, ymin = 1.99) +
  coord_cartesian(xlim = c(400, 700), ylim = c(0, 2)) +
  mytheme

spec # dimensions: 4 x 7 in

#### 3.   Clean up ####
detach(package:ggspectra)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")
