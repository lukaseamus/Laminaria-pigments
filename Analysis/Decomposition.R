##############################################################
##### Project: Pigment composition of Laminaria detritus #####
##### Script purpose: Analysis of decomposition          #####
##### Author: Luka Seamus Wright                         #####
##############################################################

#### 1.   Data preparation ####
#### 1.1  Load data ####
D <- read.csv("~/Desktop/Plymouth University/Dissertation/Pigments/Data/Decomposition.csv")

#### 1.2  Reorder levels of species factor ####
D <- within(D,{
  species <- factor(species, levels = c("o","h","d"))
})

#### 1.3  Rename variables ####
loss <- D$loss # biomass loss (g d-1)
sp <- D$species
age <- D$age
bag <- D$bag

#### 2.   Data analysis ####
#### 2.1  Determine random components ####
require(lme4)
m1 <- lm(loss ~ sp * age) # fixed effects model
m2 <- lmer(loss ~ sp * age + (age|bag), REML = F) # mixed effects model
anova(m2, m1) # continue with m1

#### 2.2  Determine fixed components ####
m3 <- update(m1, .~. - sp : age) # remove interaction
anova(m3, m1) # m3 is the better model

#### 2.3  Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m3, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m3) ~ age) # residual variance slightly increases with age
plot(resid(m3) ~ sp) # and barely with species
# overall quite homogenous

hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3))
# overall perfectly normal

par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)
# m3 is chosen as the optimal model based on its good homogeneity
# and perfect normality

#### 2.4  Interpret model ####
require(car)
Anova(m3, type = 2) # Type II Sums of Squares test
# Response: loss
#           Sum Sq Df F value    Pr(>F)    
# sp        15.778  2   6.497   0.00247 ** 
# age       60.436  1  49.771 6.543e-10 ***
# Residuals 93.499 77 

summary(m3)
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  4.38378    0.42138  10.403 2.47e-16 ***
# sph         -1.02929    0.29991  -3.432 0.000967 ***
# spd         -0.80099    0.29991  -2.671 0.009229 ** 
# age         -0.11010    0.01561  -7.055 6.54e-10 ***

sp <- factor(sp, levels = c("h", "d", "o"))
m3 <- lm(loss ~ sp + age)
summary(m3)
#               Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  3.35450    0.42138   7.961 1.21e-11 ***
# spd          0.22830    0.29991   0.761 0.448852    
# spo          1.02929    0.29991   3.432 0.000967 ***
# age         -0.11010    0.01561  -7.055 6.54e-10 ***

sp <- factor(sp, levels = c("o", "h", "d"))


#### 3.   Descriptive statistics ####
require(psych)
loss.stat <- describeBy(loss, sp, mat = T, digits = 10)

#### 4.   Clean up ####
detach(package:car)
detach(package:lme4)
detach(package:psych)
rm(list = ls())
graphics.off()
cat("\014")
