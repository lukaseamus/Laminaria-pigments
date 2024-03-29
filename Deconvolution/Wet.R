##################################################################
##### Project: Pigment composition of Laminaria detritus     #####
##### Script purpose: Spectral deconvolution of wet extracts #####
##### Author: Luka Seamus Wright                             #####
##################################################################

# Modified from Thrane, J.E., Kyle, M., Striebel, M., Haande, S., 
# Grung, M., Rohrlack, T., Andersen, T., 2015. Spectrophotometric 
# analysis of pigments: a critical assessment of a high-throughput 
# method for analysis of algal pigment mixtures by spectral deconvolution. 
# PLoS One 10, e0137645. 10.1371/journal.pone.0137645. 

source("~/PATH/pigment.function.R")
pigment.objects <- ls()
w <- 400:700

# Wet spectra of kelps collected at Mount Batten 
# Wavelengths in one nm increments
# Columns are replicates of different species and individuals
s <- read.csv("~/PATH/Wet.csv", header = T)
# Use wavelengths from 400-700 nm
s <- subset(s, wavelength >= 400 & wavelength <= 700)[, -1]
head(s)

# Fit crude spectrum as weighted sum of Gaussian peak spectra
fit <- pigment.fit(w, s, m.bias = 0, s.bias = 1)

### Root mean squared deviances 
summary(sqrt(sapply(fit$m, deviance) / 301))

# Extract fitted spectra
fitted <- fitted.spectrum(fit)
background <- background.spectrum(fit)
pigment <- pigment.spectrum(fit)

# Pathlength in crystal cuvette (z[cm] = V[cm^3]/A[cm^2]) 
z = 1

# Calculate pigment concentrations in extract
mg.L <- pigment.concentration(fit, pathl = z)

# Read sample weights
mass <- read.csv("~/PATH/wet.mass.csv", header = T)

# Samples weighed ~1 g (see wet.mass.csv) and extraction volume at 25 ml; dilution is 1.
ug.g  <- sweep(mg.L, 2, 25 / mass$g, "*") # --> ug pigment / g dry tissue

# Write concentrations to csv file
write.csv(ug.g, "~/PATH/wet.pigments.csv")
