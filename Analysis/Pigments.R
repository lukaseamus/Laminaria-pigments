##################################################################
##### Project: Pigment composition of Laminaria detritus     #####
##### Script purpose: Analysis of pigment composition        #####
##### Author: Luka Seamus Wright                             #####
##################################################################

#### 1.    Data preparation ####
#### 1.1   Load data ####
P <- read.csv("~/PATH/Pigments.csv")

#### 1.2   Reorder levels of species factor ####
P <- within(P,{
  species <- factor(species, levels = c("o","h","d"))
})

#### 1.3   Extract data from West Hoe ####
WH <- P[49:129,]

#### 1.4   Rename variables ####
chla <- WH$Chl.a # Chlorophyll a (μg gDW-1)
chlc1 <- WH$Chl.c1 # Chlorophyll c1 (μg gDW-1)
chlc2 <- WH$Chl.c2 # Chlorophyll c2 (μg gDW-1)
chlc <- chlc1 + chlc2 # Chlorophyll c (μg gDW-1)
fuco <- WH$Fuco # Fucoxanthin (μg gDW-1)
ant.chla <- (fuco + chlc)/chla # antenna:chlorophyll a ratio
bbcar <- WH$bb.Car # β,β-Carotene (μg gDW-1)
zea <- WH$Zea # Zeaxanthin (μg gDW-1)
caro <- bbcar + zea # minor carotenoids (μg gDW-1)
tot <- chla + chlc + fuco + caro # total pigments (μg gDW-1)
age <- WH$age # tissue age in days
sp <- WH$species
bag <- WH$bag

#### 1.5   Extract multivariate dataframes ####
p <- data.frame(WH[,c(1, 6:13)], row.names = 1)

#### 2.    Multivariate data analysis WH ####
#### 2.1   Calculate Euclidian dissimilarity ####
require(vegan)
dist <- vegdist(p, method = "euclidian")

#### 2.2   Test for homogeneity of dispersion ####
betad <- betadisper(dist, sp) # dispersion by species
permutest(betad, permutations = 9999) # homogenous

betad <- betadisper(dist, age) # dispersion by age
permutest(betad, permutations = 9999) # homogenous

#### 2.3   PERMANOVA ####
adonis(p ~ sp * age, permutations = 9999, method = "euclidian")
#           Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# sp         2   5992366 2996183 28.3951 0.42535 0.0001 ***
# age        1    104041  104041  0.9860 0.00738 0.3220    
# sp:age     2     78017   39008  0.3697 0.00554 0.7204    
# Residuals 75   7913819  105518         0.56173           
# Total     80  14088243                 1.00000 


#### 3.    Univariate data analysis ####
#### 3.1   Chlorophyll a ####
#### 3.1.1 Determine random components ####
require(lme4)
m1 <- lm(chla ~ sp * age)
m2 <- lmer(chla ~ sp * age + (age|bag), REML = F) 
# singular fit already indicates that model is too complex
anova(m2, m1) # continue with m1

#### 3.1.2 Determine fixed components ####
drop1(m1, test = "F") # remove interaction
m1 <- lm(chla ~ sp + age)

#### 3.1.3 Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m1, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m1) ~ age) # variance is smaller in old detritus
plot(resid(m1) ~ sp) # variance differs between species
# homogeneity could be improved

hist(resid(m1))
qqnorm(resid(m1))
qqline(resid(m1))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.1.4 Determine random components ####
require(nlme)
m3 <- gls(chla ~ sp * age,
          weights = varIdent(form = ~1|sp))
m4 <- lme(chla ~ sp * age,
          random = ~1|bag,
          weights = varIdent(form = ~1|sp)) # random intercept model
# is run here because the random intercept and slope model is too complex
anova(m4, m3) # continue with m3
m3 <- gls(chla ~ sp * age,
          weights = varIdent(form = ~1|sp),
          method = "ML")

#### 3.1.5 Determine fixed components ####
drop1(m3, test = "Chisq") # remove interaction
m3 <- gls(chla ~ sp + age,
          weights = varIdent(form = ~1|sp),
          method = "REML")

#### 3.1.6 Test model fit ####
plot(m3, col = sp)
plot(resid(m3, type = "normalized") ~ sp) 
# variance no longer differs dramatically between species
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m3))
qqnorm(resid(m3))
qqline(resid(m3))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.1.7 Interpret model ####
require(car)
Anova(m3, type = 2) # Type II sums of squares test
# Response: chla
#     Df   Chisq Pr(>Chisq)    
# sp   2 58.8248  1.684e-13 ***
# age  1  0.0635     0.8011  

summary(m3)
# L. och. vs. L. hyp. intercept, t = 5.610711, p < 0.001 ***
# L. och. vs. L. dig. intercept, t = 7.345089, p < 0.001 ***

sp <- factor(sp, levels = c("d", "h", "o"))
m3 <- gls(chla ~ sp + age,
          weights = varIdent(form = ~1|sp),
          method = "REML")
summary(m3)
# L. dig. vs. L. hyp. intercept, t = -0.568282, p = 0.57

sp <- factor(sp, levels = c("o", "h", "d"))

#### 3.2   Chlorophyll c ####
#### 3.2.1 Determine random components ####
m5 <- lm(chlc ~ sp * age)
m6 <- lmer(chlc ~ sp * age + (age|bag), REML = F)
anova(m6, m5) # continue with m5

#### 3.1.2 Determine fixed components ####
drop1(m5, test = "F") # remove interaction
m5 <- lm(chlc ~ sp + age)

#### 3.1.3 Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m5, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m5) ~ age) # variance decreases with age
plot(resid(m5) ~ sp) # variance differs between species
# homogeneity could be improved

hist(resid(m5))
qqnorm(resid(m5))
qqline(resid(m5))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.2.4 Determine random components ####
m7 <- gls(chlc ~ sp * age,
          weights = varComb(varIdent(form = ~1|sp), 
                            varExp(form = ~age)))
m8 <- lme(chlc ~ sp * age,
          random = ~1|bag,
          weights = varComb(varIdent(form = ~1|sp), 
                            varExp(form = ~age))) # random intercept model
# is run here because the random intercept and slope model is too complex
anova(m8, m7) # continue with m7
m7 <- gls(chlc ~ sp * age,
          weights = varComb(varIdent(form = ~1|sp), 
                            varExp(form = ~age)),
          method = "ML")

#### 3.2.5 Determine fixed components ####
drop1(m7, test = "Chisq") # remove interaction
m7 <- gls(chlc ~ sp + age,
          weights = varComb(varIdent(form = ~1|sp), 
                            varExp(form = ~age)),
          method = "REML")

#### 3.2.6 Test model fit ####
plot(m7, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m7, type = "normalized") ~ age) # variance no longer decreases with age
plot(resid(m7, type = "normalized") ~ sp) # variance no longer differs between species
# homogeneity is improved

hist(resid(m7))
qqnorm(resid(m7))
qqline(resid(m7))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.2.7 Interpret model ####
Anova(m7, type = 2) # Type II sums of squares test
# Response: chlc
#     Df   Chisq Pr(>Chisq)   
# sp   2 11.8018   0.002737 **
# age  1  2.2097   0.137146

summary(m7)
# L. och. vs. L. hyp. intercept, t = 2.819991, p = 0.006 **
# L. och. vs. L. dig. intercept, t = 3.214204, p = 0.002 **

sp <- factor(sp, levels = c("d", "h", "o"))
m7 <- gls(chlc ~ sp + age,
          weights = varComb(varIdent(form = ~1|sp), 
                            varExp(form = ~age)),
          method = "REML")
summary(m7)
# L. dig. vs. L. hyp. intercept, t = 0.430170, p = 0.67

sp <- factor(sp, levels = c("o", "h", "d"))


#### 3.3   Fucoxanthin ####
#### 3.3.1 Determine random components ####
m9 <- lm(fuco ~ sp * age)
m10 <- lmer(fuco ~ sp * age + (age|bag), REML = F)
anova(m10, m9) # continue with m9

#### 3.3.2 Determine fixed components ####
drop1(m9, test = "F") # remove interaction
m9 <- lm(fuco ~ sp + age)

#### 3.3.3 Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m9, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m9) ~ age) # variance is quite similar across age
plot(resid(m9) ~ sp) # variance differs between species
# homogeneity could be improved

hist(resid(m9))
qqnorm(resid(m9))
qqline(resid(m9))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.3.4 Determine random components ####
m11 <- gls(fuco ~ sp * age,
           weights = varIdent(form = ~1|sp))
m12 <- lme(fuco ~ sp * age,
           random = ~age|bag,
           weights = varIdent(form = ~1|sp))
anova(m12, m11) # continue with m11
m11 <- gls(fuco ~ sp * age,
           weights = varIdent(form = ~1|sp),
           method = "ML")

#### 3.3.5 Determine fixed components ####
drop1(m11, test = "Chisq") # remove interaction
m11 <- gls(fuco ~ sp + age,
           weights = varIdent(form = ~1|sp),
           method = "REML")

#### 3.3.6 Test model fit ####
plot(m11, col = sp)
plot(resid(m11, type = "normalized") ~ sp)
# variance no longer differs between species
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m11))
qqnorm(resid(m11))
qqline(resid(m11))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.3.7 Interpret model ####
Anova(m11, type = 2) # Type II sums of squares test
# Response: fuco
#     Df   Chisq Pr(>Chisq)    
# sp   2 63.3250  1.775e-14 ***
# age  1  2.5085     0.1132 

summary(m11)
# L. och. vs. L. hyp. intercept, t = 5.706148, p < 0.001 ***
# L. och. vs. L. dig. intercept, t = 7.824223, p < 0.001 ***

sp <- factor(sp, levels = c("d", "h", "o"))
m11 <- gls(fuco ~ sp + age,
           weights = varIdent(form = ~1|sp),
           method = "REML")
summary(m11)
# L. dig. vs. L. hyp. intercept, t = -1.943501, p = 0.06

sp <- factor(sp, levels = c("o", "h", "d"))

#### 3.4   Antenna/Chl a ####
#### 3.4.1 Determine random components ####
m13 <- lm(ant.chla ~ sp * age)
m14 <- lmer(ant.chla ~ sp * age + (age|bag), REML = F)
anova(m14, m13) # continue with m13

#### 3.4.2 Determine fixed components ####
drop1(m13, test = "F") # retain interaction

#### 3.4.3 Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m13, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m13) ~ age) # variance is quite similar across age
plot(resid(m13) ~ sp) # variance slightly differs between species
# homogeneity could be improved

hist(resid(m13))
qqnorm(resid(m13))
qqline(resid(m13))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.4.4 Determine random components ####
m15 <- gls(ant.chla ~ sp * age,
           weights = varExp(form = ~age|sp))
m16 <- lme(ant.chla ~ sp * age,
           random = ~age|bag,
           weights = varExp(form = ~age|sp))
anova(m16, m15) # continue with m15
m15 <- gls(ant.chla ~ sp * age,
           weights = varExp(form = ~age|sp),
           method = "ML")

#### 3.4.5 Determine fixed components ####
drop1(m15, test = "Chisq") # retain interaction
m15 <- gls(ant.chla ~ sp * age,
           weights = varExp(form = ~age|sp),
           method = "REML")

#### 3.4.6 Test model fit ####
plot(m15, col = sp)
plot(resid(m15, type = "normalized") ~ sp)
# variance no longer differs between species
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m15))
qqnorm(resid(m15))
qqline(resid(m15))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.4.7 Interpret model ####
Anova(m15, type = 3) # Type III sums of squares test
# Response: ant.chla
#             Df    Chisq Pr(>Chisq)    
# (Intercept)  1 213.4088  < 2.2e-16 ***
# sp           2  10.0035  0.0067263 ** 
# age          1   6.2227  0.0126120 *  
# sp:age       2  14.4305  0.0007353 ***

summary(m15)
# L. och. vs. L. hyp. intercept, t = 2.402400, p = 0.02 *
# L. och. vs. L. dig. intercept, t = -0.334836, p = 0.74
# y = 0.005x + 0.62

sp <- factor(sp, levels = c("d", "h", "o"))
m15 <- gls(ant.chla ~ sp * age,
           weights = varExp(form = ~age|sp),
           method = "REML")
Anova(m15, type = 3)
# Response: ant.chla
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 250.171  < 2.2e-16 ***
# sp           2  10.003  0.0067264 ** 
# age          1  26.393  2.785e-07 ***
# sp:age       2  14.431  0.0007353 ***

summary(m15)
# L. dig. vs. L. hyp. intercept, t = 2.308643, p = 0.005 **
# y = 0.008x + 0.6

sp <- factor(sp, levels = c("h", "d", "o"))
m15 <- gls(ant.chla ~ sp * age,
           weights = varExp(form = ~age|sp),
           method = "REML")
Anova(m15, type = 3)
# Response: ant.chla
#             Df    Chisq Pr(>Chisq)    
# (Intercept)  1 432.1003  < 2.2e-16 ***
# sp           2  10.0034  0.0067264 ** 
# age          1   0.0034  0.9536756    
# sp:age       2  14.4305  0.0007353 ***
  
summary(m15)
# y = -0.00009x + 0.75

sp <- factor(sp, levels = c("o", "h", "d"))


#### 3.5   Minor carotenoids ####
#### 3.5.1 Determine random components ####
m17 <- lm(caro ~ sp * age)
m18 <- lmer(caro ~ sp * age + (age|bag), REML = F)
anova(m18, m17) # continue with m17

#### 3.5.2 Determine fixed components ####
drop1(m17, test = "F") # remove interaction
m17 <- lm(caro ~ sp + age)

#### 3.5.3 Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m17, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m17) ~ age) # variance is quite similar across age
plot(resid(m17) ~ sp) # variance differs dramatically between species
# homogeneity could be improved

hist(resid(m17))
qqnorm(resid(m17))
qqline(resid(m17))
# quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.5.4 Determine random components ####
m19 <- gls(caro ~ sp * age,
           weights = varIdent(form = ~1|sp))
m20 <- lme(caro ~ sp * age,
           random = ~age|bag,
           weights = varIdent(form = ~1|sp))
anova(m20, m19) # continue with m19
m19 <- gls(caro ~ sp * age,
           weights = varIdent(form = ~1|sp),
           method = "ML")

#### 3.5.5 Determine fixed components ####
drop1(m19, test = "Chisq") # retain interaction
m19 <- gls(caro ~ sp * age,
           weights = varIdent(form = ~1|sp),
           method = "REML")

#### 3.5.6 Test model fit ####
plot(m19, col = sp)
plot(resid(m19, type = "normalized") ~ sp)
# variance no longer differs dramatically between species
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m19))
qqnorm(resid(m19))
qqline(resid(m19))
# quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.5.7 Interpret model ####
Anova(m19, type = 3) # Type III sums of squares test
# Response: caro
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1  5.4289    0.01981 *  
# sp           2 26.8707  1.462e-06 ***
# age          1  0.5089    0.47564    
# sp:age       2  6.5409    0.03799 * 

summary(m19)
# L. och. vs. L. hyp. intercept, t = 2.827272, p = 0.006 **
# L. och. vs. L. dig. intercept, t = 4.828576, p < 0.001 ***
# y = -0.28x + 22.55

sp <- factor(sp, levels = c("d", "h", "o"))
m19 <- gls(caro ~ sp * age,
           weights = varIdent(form = ~1|sp),
           method = "REML")
Anova(m19, type = 3)
# Response: caro
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 50.2022  1.387e-12 ***
# sp           2 26.8707  1.462e-06 ***
# age          1 10.7145   0.001063 ** 
# sp:age       2  6.5409   0.037990 *

summary(m19)
# L. dig. vs. L. hyp. intercept, t = -0.692860, p = 0.49
# y = —2.1x + 112.13

sp <- factor(sp, levels = c("h", "d", "o"))
m19 <- gls(caro ~ sp * age,
           weights = varIdent(form = ~1|sp),
           method = "REML")
Anova(m19, type = 3)
# Response: caro
#             Df   Chisq Pr(>Chisq)    
# (Intercept)  1 16.4382  5.026e-05 ***
# sp           2 26.8707  1.462e-06 ***
# age          1  0.0019    0.96521    
# sp:age       2  6.5409    0.03799 * 

summary(m19)
# y = 0.04x + 92.84

sp <- factor(sp, levels = c("o", "h", "d"))


#### 3.6   Total pigment ####
#### 3.6.1 Determine random components ####
m21 <- lm(tot ~ sp * age)
m22 <- lmer(tot ~ sp * age + (age|bag), REML = F)
anova(m22, m21) # continue with m21

#### 3.6.2 Determine fixed components ####
drop1(m21, test = "F") # remove interaction
m21 <- lm(tot ~ sp + age)

#### 3.6.3 Test model fit ####
par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m21, col = sp)
par(mfrow = c(1,2), mar = c(2,2,2,1))
plot(resid(m21) ~ age) # variance is smaller for old detritus
plot(resid(m21) ~ sp) # variance differs between species
# homogeneity could be improved

hist(resid(m21))
qqnorm(resid(m21))
qqline(resid(m21))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.6.4 Determine random components ####
m23 <- gls(tot ~ sp * age,
           weights = varIdent(form = ~1|sp))
m24 <- lme(tot ~ sp * age,
           random = ~1|bag,
           weights = varIdent(form = ~1|sp))
anova(m24, m23) # continue with m23
m23 <- gls(tot ~ sp * age,
           weights = varIdent(form = ~1|sp),
           method = "ML")

#### 3.6.5 Determine fixed components ####
drop1(m23, test = "Chisq") # remove interaction
m23 <- gls(tot ~ sp + age,
           weights = varIdent(form = ~1|sp),
           method = "REML")

#### 3.6.6 Test model fit ####
plot(m23, col = sp)
plot(resid(m23, type = "normalized") ~ sp)
# variance no longer differs dramatically between species
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m23))
qqnorm(resid(m23))
qqline(resid(m23))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.6.7 Interpret model ####
Anova(m23, type = 2) # Type II sums of squares test
# Response: tot
#     Df   Chisq Pr(>Chisq)    
# sp   2 64.5623  9.561e-15 ***
# age  1  0.1187     0.7305  

summary(m23)
# L. och. vs. L. hyp. intercept, t = 5.936001, p < 0.001 ***
# L. och. vs. L. dig. intercept, t = 7.794723, p < 0.001 ***

sp <- factor(sp, levels = c("d", "h", "o"))
m23 <- gls(tot ~ sp + age,
           weights = varIdent(form = ~1|sp),
           method = "REML")
summary(m23)
# L. dig. vs. L. hyp. intercept, t = -0.783523, p = 0.44

sp <- factor(sp, levels = c("o", "h", "d"))

#### 3.7   Antenna/Chl a (single explanatory variable) ####
#### 3.7.1 Test model fit ####
m25 <- lm(ant.chla ~ sp)

par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m25, col = sp)
par(mfrow = c(1,1), mar = c(2,2,2,1))
plot(resid(m25) ~ sp) # variance differs between species
# homogeneity could be improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m25))
qqnorm(resid(m25))
qqline(resid(m25))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.7.2 Test model fit ####
m26 <- gls(ant.chla ~ sp,
           weights = varIdent(form = ~1|sp))
plot(m26, col = sp)
plot(resid(m26, type = "normalized") ~ sp)
# variance no longer differs between species
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m26))
qqnorm(resid(m26))
qqline(resid(m26))
# normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.7.3 Interpret model ####
Anova(m26, type = 2) # Type II sums of squares test
# Response: ant.chla
#    Df  Chisq Pr(>Chisq)  
# sp  2 8.0043    0.01828 *

summary(m26)
# L. och. vs. L. hyp. intercept, t = 1.03572, p = 0.3
# L. och. vs. L. dig. intercept, t = 2.76841, p = 0.007 **

sp <- factor(sp, levels = c("d", "h", "o"))
m26 <- gls(ant.chla ~ sp,
           weights = varIdent(form = ~1|sp))
summary(m26)
# L. dig. vs. L. hyp. intercept, t = -2.12481, p = 0.04 *

sp <- factor(sp, levels = c("o", "h", "d"))

#### 3.8   Minor carotenoids (single explanatory variable) ####
#### 3.8.1 Test model fit ####
m27 <- lm(caro ~ sp)

par(mfrow = c(2,2), mar = c(2,2,2,1))
plot(m27, col = sp)
par(mfrow = c(1,1), mar = c(2,2,2,1))
plot(resid(m27) ~ sp) # variance differs between species
# homogeneity could be improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m27))
qqnorm(resid(m27))
qqline(resid(m27))
# quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.8.2 Test model fit ####
m28 <- gls(caro ~ sp,
           weights = varIdent(form = ~1|sp))

plot(m28, col = sp)
plot(resid(m28, type = "normalized") ~ sp)
# variance no longer differs dramatically between species
# homogeneity is improved

par(mfrow = c(1,2), mar = c(2,2,2,1))
hist(resid(m28))
qqnorm(resid(m28))
qqline(resid(m28))
# quite normal
par(mfrow = c(1,1), mar = c(5, 4, 4, 2) + 0.1)

#### 3.8.3 Interpret model ####
Anova(m28, type = 2) # Type II sums of squares test
# Response: caro
#    Df  Chisq Pr(>Chisq)    
# sp  2 128.54  < 2.2e-16 ***

summary(m28)
# L. och. vs. L. hyp. intercept, t = 9.995114, p < 0.001 ***
# L. och. vs. L. dig. intercept, t = 7.065885, p < 0.001 ***

sp <- factor(sp, levels = c("d", "h", "o"))
m28 <- gls(caro ~ sp,
           weights = varIdent(form = ~1|sp))
summary(m28)
# L. dig. vs. L. hyp. intercept, t = 3.313858, p = 0.001 **

sp <- factor(sp, levels = c("o", "h", "d"))

#### 4.    Data visualisation ####
#### 4.1   Calculate descriptive statistics ####
require(psych)
tot.stat <- describeBy(tot, sp, mat = T)
chla.stat <- describeBy(chla, sp, mat = T)
chlc1.stat <- describeBy(chlc1, sp, mat = T)
chlc2.stat <- describeBy(chlc2, sp, mat = T)
chlc.stat <- describeBy(chlc, sp, mat = T)
fuco.stat <- describeBy(fuco, sp, mat = T)
zea.stat <- describeBy(zea, sp, mat = T)
bbcar.stat <- describeBy(bbcar, sp, mat = T)
caro.stat <- describeBy(caro, sp, mat = T)
ant.stat <- describeBy(ant.chla, sp, mat = T)
caro.stat2 <- describeBy(caro, list(sp, age), mat = T)
caro.stat2$group2 <- as.integer(caro.stat2$group2)
ant.stat2 <- describeBy(ant.chla, list(sp, age), mat = T)
ant.stat2$group2 <- as.integer(ant.stat2$group2)

#### 4.2   Calculate model predictions ####
# only models with significant slopes are considered
sp <- factor(sp, levels = c("d", "h", "o"))
m15 <- gls(ant.chla ~ sp * age,
           weights = varExp(form = ~age|sp),
           method = "REML")
m19 <- gls(caro ~ sp * age,
           weights = varIdent(form = ~1|sp),
           method = "REML")

new <- data.frame(age = 13:32,
                  sp = c(rep("d", 20), rep("h", 20), rep("o", 20)))

new$ant.fit <- predict(m15, newdata = new)
modmat <-  model.matrix(formula(m15)[-2], new)
int <- diag(modmat %*% vcov(m15) %*% t(modmat))
new$ant.lo <- with(new, ant.fit - 1.96*sqrt(int))
new$ant.hi <- with(new, ant.fit + 1.96*sqrt(int))

new$caro.fit <- predict(m19, newdata = new)
modmat <-  model.matrix(formula(m19)[-2], new)
int <- diag(modmat %*% vcov(m19) %*% t(modmat))
new$caro.lo <- with(new, caro.fit - 1.96*sqrt(int))
new$caro.hi <- with(new, caro.fit + 1.96*sqrt(int))

#### 4.3  Calculate nMDS data ####
nMDS <- metaMDS(p, distance = "euclidian", autotransform = FALSE)
nMDS$stress # Stress = 0.02

#### 4.4  Extract points and calculate 95% confidence ellipses ####
df <- data.frame(nMDS$points, group1 = sp, group2 = age)
df$group1 <- factor(df$group1, levels = c("d", "h", "o"))
ellipse <- ordiellipse(nMDS, sp, display = "sites", kind = "se",
                       conf = 0.95, label = F)

#### 4.5  Make data compatible with ggplot2 ####
# https://github.com/jarioksa/vegan/blob/master/R/veganCovEllipse.R
veganCovEllipse <- function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  Q <- chol(cov, pivot = TRUE)
  o <- attr(Q, "pivot")
  t(center + scale * t(Circle %*% Q[,o]))
}

dfell <- data.frame()
for(g in levels(df$group1)) {
  dfell <- rbind(dfell, 
                 cbind(as.data.frame(with(df[df$group1==g,],
                                          veganCovEllipse(ellipse[[g]]$cov,
                                                          ellipse[[g]]$center,
                                                          ellipse[[g]]$scale))),
                       group1 = g))
}


#### 4.6  Define theme ####
require(ggplot2)
mytheme <- theme(panel.background = element_blank(),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(fill = NA, size = 1),
                 axis.line = element_blank(),
                 axis.title = element_blank(),
                 axis.text = element_blank(),
                 axis.ticks = element_blank(),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 legend.text = element_text(size = 12),
                 legend.text.align = 0,
                 legend.title = element_blank(),
                 legend.margin = margin(t = -7, r = -7, b = -7, l = -7, 
                                        unit = "pt"),
                 text = element_text(family = "Helvetica Neue"))


#### 4.7  nMDS plot ####
np <- ggplot(data = df, aes(MDS1, MDS2)) +
                geom_point(aes(colour = group1, fill = group1, shape = as.factor(group2)),
                           size = 2.5) +
                geom_polygon(data = dfell, aes(NMDS1, NMDS2, colour = group1,
                                                 fill = group1), size = 0.5, alpha = 0.5) +
                scale_fill_manual(values = c("#333b08","#81a512","#f1c700"),
                                  labels = c(expression(italic("L. digitata")),
                                             expression(italic("L. hyperborea")),
                                             expression(italic("L. ochroleuca"))),
                                  guide = guide_legend(order = 1)) +
                scale_colour_manual(values = c("#333b08","#81a512","#f1c700"),
                                    labels = c(expression(italic("L. digitata")),
                                               expression(italic("L. hyperborea")),
                                               expression(italic("L. ochroleuca"))),
                                    guide = guide_legend(order = 1)) +
                scale_shape_discrete(labels = c("13 d", "25 d", "32 d"),
                                     guide = guide_legend()) +
                annotate("text", label = "Stress = 0.02", x = 800, y = -200,
                         size = 4.2) +
                scale_y_reverse() +
                mytheme +
                theme(legend.position = c(0.19, 0.8))

np # dimensions: 4 x 4 in

#### 4.8   Univariate plots ####
#### 4.8.1 Redefine theme ####
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

#### 4.8.2 Antenna/Chl a ####
ap <- ggplot() +
        geom_line(data = new, aes(age, ant.fit, colour = sp, lty = sp), size = 0.5) +
        geom_ribbon(data = new, aes(age, ymin = ant.lo, ymax = ant.hi, fill = sp),
                    alpha = 0.5) +
        geom_pointrange(data = ant.stat2, aes(group2, mean, ymin = mean - se,
                                              ymax = mean + se, colour = group1),
                        size = 0.5, position = position_dodge(width = 1.5)) +
        scale_colour_manual(values = c("#333b08","#81a512","#f1c700"),
                            labels = c(expression(italic("L. digitata")*"             y = 0.008x + 0.6"),
                                       expression(italic("L. hyperborea")*"       y = —0.00009x + 0.75"),
                                       expression(italic("L. ochroleuca")*"       y = 0.005x + 0.62")),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#333b08","#81a512","#f1c700"),
                          labels = c(expression(italic("L. digitata")*"             y = 0.008x + 0.6"),
                                     expression(italic("L. hyperborea")*"       y = —0.00009x + 0.75"),
                                     expression(italic("L. ochroleuca")*"       y = 0.005x + 0.62")),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5, 1),
                              guide = F) +
        ylab(expression("Antenna pigment : chlorophyll "*italic(a))) +
        xlab("Detrital age (d)") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        coord_cartesian(ylim = c(0.6, 1), xlim = c(10, 35)) +
        theme(legend.position = c(0.27, 0.9)) +
        mytheme

ap # dimsenisons: 4 x 7 in

#### 4.8.3 Minor carotenoids ####
cp <- ggplot() +
        geom_line(data = new, aes(age, caro.fit, colour = sp, lty = sp), size = 0.5) +
        geom_ribbon(data = new, aes(age, ymin = caro.lo, ymax = caro.hi, fill = sp),
                    alpha = 0.5) +
        geom_pointrange(data = caro.stat2, aes(group2, mean, ymin = mean - se,
                                              ymax = mean + se, colour = group1),
                        size = 0.5, position = position_dodge(width = 1.5)) +
        scale_colour_manual(values = c("#333b08","#81a512","#f1c700"),
                            labels = c(expression(italic("L. digitata")*"             y = —2.1x + 112.13"),
                                       expression(italic("L. hyperborea")*"       y = 0.04x + 92.84"),
                                       expression(italic("L. ochroleuca")*"       y = —0.28x + 22.55")),
                            guide = guide_legend()) +
        scale_fill_manual(values = c("#333b08","#81a512","#f1c700"),
                          labels = c(expression(italic("L. digitata")*"             y = —2.1x + 112.13"),
                                     expression(italic("L. hyperborea")*"       y = 0.04x + 92.84"),
                                     expression(italic("L. ochroleuca")*"       y = —0.28x + 22.55")),
                          guide = guide_legend()) +
        scale_linetype_manual(values = c(1, 5, 5),
                              guide = F) +
        ylab(expression("Minor carotenoids ("*mu*"g g"^-1*")")) +
        xlab("Detrital age (d)") +
        scale_y_continuous(expand = c(0,0)) +
        scale_x_continuous(expand = c(0,0)) +
        coord_cartesian(ylim = c(0, 150), xlim = c(10, 35)) +
        theme(legend.position = c(0.27, 0.9)) +
        mytheme

cp # dimsenisons: 4 x 7 in

#### 5.   Clean up ####
detach(package:vegan)
detach(package:car)
detach(package:lme4)
detach(package:nlme)
detach(package:psych)
detach(package:ggplot2)
rm(list = ls())
graphics.off()
cat("\014")

