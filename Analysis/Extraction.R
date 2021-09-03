##################################################################
##### Project: Pigment composition of Laminaria detritus     #####
##### Script purpose: Comparison of extraction techniques    #####
##### Author: Luka Seamus Wright                             #####
##################################################################

#### 1.    Data preparation ####
#### 1.1   Load data ####
E <- read.csv("~/Desktop/Plymouth University/Dissertation/Pigments/Data/Extraction.csv")

#### 1.2   Rename variables ####
wet <- E$wet # wet extraction concentration (μg gWW-1)
dry <- E$dry # dry extraction concentration (μg gDW-1)
pig <- E$pigment

#### 2.    Data analysis ####
#### 2.1   Determine fixed components ####
m1 <- lm(log10(wet) ~ log10(dry) * pig) 
# log10 transformation of both axes is necessary because concentrations are so 
# different between different pigments 
drop1(m1, test = "F") # retain interaction

#### 2.2   Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m1, pch = pig)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m1) ~ log10(dry)) # variance is smaller for low concentrations
boxplot(resid(m1) ~ pig) # variance differs between pigments
# homogeneity could be improved

hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1))
# quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 2.3   Determine variance structure ####
require(nlme)
m2 <- gls(log10(wet) ~ log10(dry) * pig,
          weights = varIdent(form = ~1|pig))
m3 <- gls(log10(wet) ~ log10(dry) * pig,
          weights = varExp(form = ~log10(dry)|pig))
anova(m3, m2) # continue with m3
m3 <- gls(log10(wet) ~ log10(dry) * pig,
          weights = varExp(form = ~log10(dry)|pig),
          method = "ML")

#### 2.4   Determine fixed components ####
drop1(m3, test = "Chisq") # retain interaction
m3 <- gls(log10(wet) ~ log10(dry) * pig,
          weights = varExp(form = ~log10(dry)|pig),
          method = "REML")

#### 2.5  Test model fit ####
plot(m3, pch = pig)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m3, type = "normalized") ~ log10(dry))
boxplot(resid(m3, type = "normalized") ~ pig)
# variance no longer differs between concentrations or pigments
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3))
# normality is similar
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m3 is chosen as the optimal model

#### 2.6  Interpret model ####
require(car)
Anova(m3, type = 3) # Type III sums of squares test
# Response: log10(wet)
#                Df   Chisq Pr(>Chisq)    
# (Intercept)     1  2.4371   0.118499    
# log10(dry)      1  7.7127   0.005483 ** 
# pig             2 19.4865  5.869e-05 ***
# log10(dry):pig  2 22.8033  1.118e-05 ***

summary(m3)
# log10(wet) = 1.367063*log10(dry) - 2.295305

pig <- factor(pig, levels = c("Chl.a", "Chl.c", "Caro"))
m3 <- gls(log10(wet) ~ log10(dry) * pig,
          weights = varExp(form = ~log10(dry)|pig),
          method = "REML")
Anova(m3, type = 3)
# Response: log10(wet)
#                Df   Chisq Pr(>Chisq)    
# (Intercept)     1  8.3531    0.00385 ** 
# log10(dry)      1 23.6667  1.145e-06 ***
# pig             2 19.4865  5.869e-05 ***
# log10(dry):pig  2 22.8033  1.118e-05 ***
        
summary(m3)
# log10(wet) = 1.511816*log10(dry) - 2.733094

pig <- factor(pig, levels = c("Chl.c", "Chl.a", "Caro"))
m3 <- gls(log10(wet) ~ log10(dry) * pig,
          weights = varExp(form = ~log10(dry)|pig),
          method = "REML")
Anova(m3, type = 3)
# Response: log10(wet)
#                Df   Chisq Pr(>Chisq)    
# (Intercept)     1 15.0861  0.0001027 ***
# log10(dry)      1  0.0487  0.8253418    
# pig             2 19.4865  5.869e-05 ***
# log10(dry):pig  2 22.8033  1.118e-05 ***

summary(m3)
# log10(wet) = 0.032633*log10(dry) + 1.125841

#### 3.    Data visualisation ####
#### 3.1   Calculate model predictions ####
pig <- factor(pig, levels = c("Chl.a", "Chl.c", "Caro"))
m3 <- gls(log10(wet) ~ log10(dry) * pig,
          weights = varExp(form = ~log10(dry)|pig),
          method = "REML")

E$fit <- predict(m3) # predicted values from m3
modmat <-  model.matrix(formula(m3)[-2]) # [-2] removes response variable from the formula
int <- diag(modmat %*% vcov(m3) %*% t(modmat))
E$lower <- with(E, fit - 1.96*sqrt(int))
E$upper <- with(E, fit + 1.96*sqrt(int))


#### 3.2   Customise theme ####
require(ggplot2)
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 plot.margin = unit(c(.3, .3, .2, .2),"cm"),
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

#### 3.3   Plot ####
wdp <- ggplot(data = E, aes(log10(dry), log10(wet))) +
       geom_point(aes(colour = pig), size = 1) +
       geom_line(aes(log10(dry), fit, colour = pig, lty = pig), size = 0.5) +
       geom_ribbon(aes(ymin = lower, ymax = upper, fill = pig),
                   alpha = .5) +
       scale_colour_manual(values = c("#0d98ba", "#50c878", "#ffae42"),
                      labels = c(expression("Chlorophyll"*italic(" a")*"       log"[10]*"(y) = 1.51log"[10]*"(x) — 2.73"),
                                 expression("Chlorophyll"*italic(" c")*"       log"[10]*"(y) = 0.03log"[10]*"(x) + 1.13"),
                                 expression("Carotenoids"*"         log"[10]*"(y) = 1.37log"[10]*"(x) — 2.3")),
                      guide = guide_legend()) +
       scale_fill_manual(values = c("#0d98ba", "#50c878", "#ffae42"),
                      labels = c(expression("Chlorophyll"*italic(" a")*"       log"[10]*"(y) = 1.51log"[10]*"(x) — 2.73"),
                                 expression("Chlorophyll"*italic(" c")*"       log"[10]*"(y) = 0.03log"[10]*"(x) + 1.13"),
                                 expression("Carotenoids"*"         log"[10]*"(y) = 1.37log"[10]*"(x) — 2.3")),
                      guide = guide_legend()) +
       scale_linetype_manual(values = c(1, 5, 1),
                             guide = F) +
       ylab(expression("Wet extract concentration ("*mu*"g g"^-1*")")) +
       xlab(expression("Dry extract concentration ("*mu*"g g"^-1*")")) +
       scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(expression("10"^1), expression("10"^2), 
                                                             expression("10"^3), expression("10"^4)), 
                          expand = c(0,0)) +
       scale_y_continuous(breaks = c(0, 1, 2, 3), labels = c(expression("10"^0), expression("10"^1), 
                                                             expression("10"^2), expression("10"^3)), 
                          expand = c(0,0)) +
       coord_cartesian(ylim = c(0, 3), xlim = c(1, 4)) +
       theme(legend.position = c(.32, .9)) +
       mytheme

wdp # print plot (dimensions: 4 x 7 in)

#### 4.   Clean up ####
detach(package:car)
detach(package:nlme)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")
