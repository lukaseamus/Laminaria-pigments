##################################################################
##### Project: Pigment composition of Laminaria detritus     #####
##### Script purpose: Analysis of environmental data         #####
##### Author: Luka Seamus Wright                             #####
##################################################################

#### 1.    Data preparation ####
#### 1.1   Load data ####
E <- read.csv("~/PATH/Environmental.csv")

#### 1.2   Rename variables ####
temp <- E$temp
lux <- E$lux
day <- E$day # unique "day" measure for each sample (i.e. day 0, 0.1, 0.2 etc.)
d <- E$d # grouped measurements (i.e. day 0, 1, 2 etc.)

#### 2.    Data exploration ####
#### 2.1   Temperature ####
## Temperature
par(mfrow = c(2,1), mar = c(2,2,2,1))
hist(temp) # slight right skew
boxplot(temp, horizontal = T) # but quite normal
mean(temp) # mean = 13.49867°C
sd(temp)/length(temp) # se = 0.0002040612°C
aggregate(temp ~ daytime, mean, data = E)

#### 2.2   Lux ####
hist(lux) # extreme right skew
boxplot(lux, horizontal = T)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
mean(lux*0.019) # mean = 0.3633402 μmol photons m-2 s-1
sd(lux*0.019)/length(lux) # se = 0.0004602058 μmol photons m-2 s-1
range(lux*0.019) # range = 0-29.45 μmol photons m-2 s-1
aggregate(lux*0.019 ~ daytime, mean, data = E) # day mean = 0.5439438 μmol photons m-2 s-1
aggregate(lux*0.019 ~ daytime, function(x){sd(x)/length(x)}, data = E) # day se = 0.0008289701 μmol photons m-2 s-1
aggregate(lux*0.019 ~ daytime, range, data = E) # day range = 0-29.45 μmol photons m-2 s-1

#### 3.    Data analysis ####
#### 3.1   Linear model ####
#### 3.1.1 Test model fit ####
m1 <- lm(temp ~ day)

par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m1) # homogeneity ok for a simple model,
# but clearly temperature is oscillating
par(mfrow = c(1,2), mar = c(2,2,2,1))

hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1)) # quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.1.2 Interpret model ####
require(car)
Anova(m1) # Type II sums of squares test
# Response: temp
#           Sum Sq   Df F value    Pr(>F)    
# day       424.22    1    1688 < 2.2e-16 ***
# Residuals 768.52 3058 

summary(m1) # y = 0.04x + 12.83

#### 3.2   Additive model ####
#### 3.2.1 Test model fit ####
require(mgcv)
m2 <- gam(temp ~ s(day, k = -1, fx = F, bs = "cr"))

plot(resid(m2) ~ fitted(m2)) # homogeneity is quite good
abline(0,0)

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m2))
qqnorm(resid(m2))
qqline(resid(m2)) # right-skewed
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.2.2 Goodness-of-fit of gamma ####
require(fitdistrplus)
gamma <- fitdist(temp, "gamma") 
norm <- fitdist(temp, "norm") 

par(mfrow = c(1,2), mar = c(2,2,2,1))
denscomp(list(gamma, norm), 
         legendtext = c("Gamma", "Normal"), 
         fitlty = 1)
cdfcomp(list(gamma, norm), 
        legendtext = c("Gamma", "Normal"), 
        fitlty = 1)
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

gofstat(list(gamma, norm), fitnames = c("Gamma", "Normal"))
# gamma distribution fits somewhat better than the normal distribution

#### 3.2.3 Test model fit ####
m3 <- gam(temp ~ s(day, k = -1, fx = F, bs = "cr"),
          family = Gamma(link = "log"))

plot(resid(m3) ~ fitted(m3)) # homogeneity is quite good
abline(0,0)

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3)) # still right-skewed
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# proceed with m2

#### 3.2.4 Interpret model ####
summary(m2)
# Approximate significance of smooth terms:
# edf Ref.df     F p-value    
# s(day) 8.887  8.996 592.1  <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.635   Deviance explained = 63.6%
# GCV = 0.14272  Scale est. = 0.14226   n = 3060

#### 4.    Data visualisation ####
#### 4.1   Calculate model predictions ####
fit <- data.frame(predict(m1, interval = "confidence"))
E$l.fit <- fit$fit
E$l.hi <- fit$upr
E$l.lo <- fit$lwr

fit <- data.frame(predict(m2, se = T, type = "response"))
E$a.fit <- fit$fit
E$a.hi <- E$a.fit + 1.96 * fit$se.fit
E$a.lo <- E$a.fit - 1.96 * fit$se.fit

#### 4.2  Define theme ####
require(ggplot2)
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.2, .3, .2, .2),"cm"),
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

#### 4.3  Plot ####
p <- ggplot(data = E, aes(x = day)) +
        geom_point(aes(y = temp), shape = 16, size = 2, alpha = 0.04) +
        geom_line(aes(y = l.fit), size = 0.5) +
        geom_ribbon(aes(ymin = l.lo, ymax = l.hi), alpha = 0.5) +
        geom_line(aes(y = a.fit), size = 0.5, colour = "#0d98ba") +
        geom_ribbon(aes(ymin = a.lo, ymax = a.hi), alpha = 0.5, fill = "#0d98ba") +
        ylab("Temperature (°C)") +
        xlab("Time (d)") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0), breaks = seq(0, 35, by = 5)) +
        coord_cartesian(ylim = c(12, 16), xlim = c(0, 35)) +
        mytheme

p # dimsenisons: 4 x 7 in

#### 5.    Clean up ####
detach(package:ggplot2)
detach(package:fitdistrplus)
detach(package:mgcv)
detach(package:car)
rm(list = ls())
graphics.off()
cat("\014")
