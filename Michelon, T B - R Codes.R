################################################################################
############################# Data visualization ###############################
################################################################################

# Experiment data
ivg <- dget("https://raw.githubusercontent.com/thomasbrunomichelon/emergence.data/main/original.data.txt")
str(ivg)
#Add = dummy factor for additional treatment (a = additional treatment, b = factor1*factor2)
#Time and Package = factors 4 x 2 (1 to 4 months x 0,10 and 0,20 mm package) + additional treatment (0 months + no package (p0))
#Diasbef/Daysaf = days between evaluation
#rep = repetitions
#germ = number of seedlings emerged between evaluations
#germac = accumulated number of seedlings emerged between evaluations

# Emergence speed index calculated from emergence test
# Made using excel = SUM(days from the beginning of the experiment to the emergence / number of emerged seedlings)
ESI <-  dget("https://raw.githubusercontent.com/thomasbrunomichelon/emergence.data/main/ESI.txt")

# Packages
library(ggplot2)
library(drc)
library(dplyr)
library(survival)
library(flexsurv)
library(rms)
library(ExpDes.pt)
library(emmeans)
library(lmtest)
library(drcSeedGerm)
library(devtools)
library(sandwich)
library(car)
library(multcomp)
# Session info

#R version 4.0.2 (2020-06-22)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 18363)

#Matrix products: default

#locale:
#[1] LC_COLLATE=Portuguese_Brazil.1252  LC_CTYPE=Portuguese_Brazil.1252    LC_MONETARY=Portuguese_Brazil.1252 LC_NUMERIC=C                      
#[5] LC_TIME=Portuguese_Brazil.1252    

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] emmeans_1.4.7   ExpDes.pt_1.2.0 rms_6.0-0       SparseM_1.78    Hmisc_4.4-0     Formula_1.2-3   lattice_0.20-41 flexsurv_1.1.1  survival_3.1-12
#[10] dplyr_0.8.5     drc_3.0-1       MASS_7.3-51.6   ggplot2_3.3.0  

#loaded via a namespace (and not attached):
#[1] nlme_3.1-148        RColorBrewer_1.1-2  tools_4.0.2         backports_1.1.6     R6_2.4.1            rpart_4.1-15        colorspace_1.4-1   
#[8] nnet_7.3-14         withr_2.2.0         tidyselect_1.1.0    gridExtra_2.3       curl_4.3            compiler_4.0.2      quantreg_5.55      
#[15] htmlTable_2.0.0     sandwich_2.5-1      scales_1.1.1        checkmate_2.0.0     polspline_1.1.19    mvtnorm_1.1-0       quadprog_1.5-8     
#[22] multcompView_0.1-8  stringr_1.4.0       digest_0.6.25       foreign_0.8-80      rmarkdown_2.3       rio_0.5.16          base64enc_0.1-3    
#[29] jpeg_0.1-8.1        pkgconfig_2.0.3     htmltools_0.5.0     plotrix_3.7-8       htmlwidgets_1.5.1   rlang_0.4.6         readxl_1.3.1       
#[36] rstudioapi_0.11     farver_2.0.3        zoo_1.8-8           gtools_3.8.2        acepack_1.4.1       zip_2.0.4           car_3.0-7          
#[43] magrittr_1.5        Matrix_1.2-18       Rcpp_1.0.5          munsell_0.5.0       abind_1.4-5         lifecycle_0.2.0     stringi_1.4.6      
#[50] multcomp_1.4-13     yaml_2.2.1          mstate_0.2.12       carData_3.0-3       plyr_1.8.6          grid_4.0.2          forcats_0.5.0      
#[57] crayon_1.3.4        haven_2.3.1         stargazer_5.2.2     splines_4.0.2       hms_0.5.3           knitr_1.29          pillar_1.4.4       
#[64] estimability_1.3    codetools_0.2-16    glue_1.4.0          evaluate_0.14       latticeExtra_0.6-29 data.table_1.12.8   deSolve_1.28       
#[71] png_0.1-7           vctrs_0.3.0         MatrixModels_0.4-1  cellranger_1.1.0    gtable_0.3.0        purrr_0.3.4         muhaz_1.2.6.1      
#[78] tidyr_1.1.0         assertthat_0.2.1    xfun_0.15           openxlsx_4.1.5      xtable_1.8-4        coda_0.19-4         tibble_3.0.1       
#[85] cluster_2.1.0       TH.data_1.0-10      ellipsis_0.3.0     


## FIGURE 1 - cumulative proportion of emerged seeds
# create a cumulative proportion (25 seeds in total)
ivg$prop <- ivg$prop <- (ivg$germac)/(25)
# Factors mean value
Data_mean <- aggregate(prop ~ diasaf:Pack:Time, data = ivg, mean)
Data_mean <- rename(Data_mean, "Packing" = "Pack")
levels(Data_mean$Time) <- c("Without storing", "1 Month", "2 Months", "3 Months", "4 Months")
levels(Data_mean$Pack) <- c("No packing", "0.10 mm", "0.20 mm")

ggplot(Data_mean, aes(x = diasaf, y = prop)) +geom_point(aes(colour = Pack), size = 0.1) + facet_wrap(.~Time) +
  scale_y_continuous(name="Cumulative proportion of emerged seeds", breaks = seq(0.0, 0.75, 0.1)) +
  scale_x_continuous(name="Time(days)", breaks = seq(20, 120, 30)) +
  theme(axis.text.x = element_text(color="black", size=12), axis.title = element_text(size = 12), axis.text.y = element_text(color="black", size=12), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  ggplot(etpred, aes(x=1:120, y=V1))+geom_line()+geom_line(aes(x=1:120, y=V4))+geom_line(aes(x=1:120, y=V3))


################################################################################
######################## Non-linear regression (NLR) ###########################
################################################################################

## Model selection (Cumulative proportion ~ days after beginning test)
ivg_dose_LN.3 <- drm(prop ~ diasaf, curveid = Time:Pack, data = ivg, fct = LN.3())

## Testing alternative distribuctions 
mselect(ivg_dose_LN.3, list(LL.4(), LL.5(), W1.3(), W1.4(), W2.4(),LN.4()))
## Lognormal distribution selected by AIC value.

# Additional treatment effect
ivg_dose_LN.3.null <- drm(prop ~ diasaf, data = ivg, fct = LN.3())
ivg_dose_LN.3.additionaltreatment  <- drm(prop ~ diasaf, curveid = Add, data = ivg, fct = LN.3())
anova(ivg_dose_LN.3.null, ivg_dose_LN.3.additionaltreatment)

## Checking model reduction
ivg_dose_LN.3.Time <- drm(prop ~ diasaf, Time, data = ivg, fct = LN.3())
ivg_dose_LN.3.Pack <- drm(prop ~ diasaf, Pack, data = ivg, fct = LN.3())

anova(ivg_dose_LN.3, ivg_dose_LN.3.null, test = "F")
anova(ivg_dose_LN.3, ivg_dose_LN.3.Time, test = "F")
anova(ivg_dose_LN.3, ivg_dose_LN.3.Pack, test = "F")
# covariates reductions can't be made

# parameters reductions
ivg_dose_LN.3.slope <- drm(prop ~ diasaf, curveid = Time:Pack, data = ivg, pmodels=list(~1, ~ivg$Time:ivg$Pack, ~ivg$Time:ivg$Pack), fct = LN.3())
ivg_dose_LN.3.asymptote <- drm(prop ~ diasaf, curveid = Time:Pack, data = ivg, pmodels=list(~ivg$Time:ivg$Pack, ~1, ~ivg$Time:ivg$Pack), fct = LN.3())
ivg_dose_LN.3.inflection <- drm(prop ~ diasaf, curveid = Time:Pack, data = ivg, pmodels=list(~ivg$Time:ivg$Pack, ~ivg$Time:ivg$Pack, ~1), fct = LN.3())

anova(ivg_dose_LN.3, ivg_dose_LN.3.slope, test = "F")
anova(ivg_dose_LN.3, ivg_dose_LN.3.asymptote, test = "F")
anova(ivg_dose_LN.3, ivg_dose_LN.3.inflection, test = "F")
# parameters reductions can't be made

# checking graphicaly
plot(ivg_dose_LN.3,
     log="",
     ylim = c(0,0.7),
     type = "average",
     lwd = 1,
     xlab = "Time(days)",
     ylab = "Proportion of germinated seeds",
     legend = F)

################################################################################
############################# Time-to-event (TTE) ##############################
################################################################################

## Model selection
#loglogistic, lognormal, weibull - three parameters
ivg_time_LL.3 <- drm(germ ~ diasbef + diasaf, curveid = Time:Pack, data = ivg, fct = LL.3(), type = "event")
ivg_time_LN.3 <- drm(germ ~ diasbef + diasaf, curveid = Time:Pack, data = ivg, fct = LN.3(), type = "event")
ivg_time_W1.3 <- drm(germ ~ diasbef + diasaf, curveid = Time:Pack, data = ivg, fct = W1.3(),  type = "event")

AIC(ivg_time_LL.3, ivg_time_LN.3, ivg_time_W1.3)

# Lognormal distribuction prodives lower AIC.

# Additional treatment effect
ivg_time_LN.3.null <- drm(germ ~ diasbef + diasaf, data = ivg, fct = LN.3(), type = "event")
ivg_time_LN.3.additionaltreatment <- drm(germ ~ diasbef + diasaf, curveid = Add, data = ivg, fct = LN.3(), type = "event")
anova(ivg_time_LN.3.additionaltreatment, ivg_time_LN.3.null, test = "F")

# Test model covariates reduction:
ivg_time_LN.3.Time <- drm(germ ~ diasbef + diasaf, curveid = Time, data = ivg, fct = LN.3(), type = "event")
ivg_time_LN.3.Pack <- drm(germ ~ diasbef + diasaf, curveid = Pack, data = ivg, fct = LN.3(), type = "event")

anova(ivg_time_LN.3, ivg_time_LN.3.null, test = "F")
anova(ivg_time_LN.3, ivg_time_LN.3.Time, test = "F")
anova(ivg_time_LN.3, ivg_time_LN.3.Pack, test = "F")
# covariates reductions can't be made.

# parameters reductions
ivg_time_LN.3.slope <- drm(germ ~ diasbef + diasaf, curveid = Time:Pack, data = ivg,  type = "event", pmodels=list(~1, ~ivg$Time:ivg$Pack, ~ivg$Time:ivg$Pack), fct = LN.3())
ivg_time_LN.3.asymptote <- drm(germ ~ diasbef + diasaf, curveid = Time:Pack, data = ivg,  type = "event", pmodels=list(~ivg$Time:ivg$Pack, ~1, ~ivg$Time:ivg$Pack), fct = LN.3())
ivg_time_LN.3.inflection <- drm(germ ~ diasbef + diasaf, curveid = Time:Pack, data = ivg, type = "event", pmodels=list(~ivg$Time:ivg$Pack, ~ivg$Time:ivg$Pack, ~1), fct = LN.3())

anova(ivg_time_LN.3, ivg_time_LN.3.slope, test = "F")
anova(ivg_time_LN.3, ivg_time_LN.3.asymptote, test = "F")
anova(ivg_time_LN.3, ivg_time_LN.3.inflection, test = "F")
# Slope (parameter b) can be the same for all the curves.

# final model:
ivg_time_LN.3.slope

# checking graphicaly
plot(ivg_time_LN.3.slope,
     log="",
     ylim = c(0,0.7),
     type = "average",
     lwd = 1,
     xlab = "Time(days)",
     ylab = "Proportion of germinated seeds",
     legend = F)

### TABLE 1 - estimated parameters and standard errors for TTE and NLR models
# Robust standard errors by sandwich method (vcovCL) to implement the cluster effect of the petri dishes
estfun.drc <- drc::estfun.drc
bread.drc <- drc::bread.drc
registerS3method("estfun", "drc", drc::estfun.drc)
registerS3method("bread", "drc", drc::bread.drc)

# NLR parameters with robust SE
coeftest(ivg_dose_LN.3, vcov= vcovCL)
# TTE parameters with robust SE
coeftest(ivg_time_LN.3.slope, vcov=vcovCL)

### TABLE 2 - Time to emerge 30% of total seeds in TTE and NRL models
# TTE - T30 - factor1(time):factor2(packing)
# *T30 was chosen because some treatments didn't reach the 50th percentile
ED(ivg_time_LN.3.slope, 0.4, interval = "delta", vcov = vcovCL, type="absolute")
# Mean comparison
TTE.T30 <- ED(ivg_time_LN.3.slope, 0.40, interval = "delta", vcov = vcovCL, multcomp = T, display = F, type="absolute")
summary(glht(TTE.T30[["EDmultcomp"]],
             linfct = contrMat(1:9, "Tukey")))

# TTE - T30 - additional treatment and interaction
ED(ivg_time_LN.3.additionaltreatment, 0.3, interval = "delta", vcov = vcovCL, type = "absolute")
TTE.add.T30 <- ED(ivg_time_LN.3.additionaltreatment, 0.3, interval = "delta", vcov = vcovCL, multcomp = T, display = F, type="absolute")
summary(glht(TTE.add.T30[["EDmultcomp"]],
             linfct = contrMat(1:2, "Tukey")))

# NLR - T30 - factor1(time):factor2(packing)
ED(ivg_dose_LN.3, 0.3, interval = "delta", vcov = vcovCL, type="absolute")
# Mean comparison
NLR.T30 <- ED(ivg_dose_LN.3, 0.30, interval = "delta", vcov = vcovCL, multcomp = T, display = F, type="absolute")
summary(glht(NLR.T30[["EDmultcomp"]],
             linfct = contrMat(1:9, "Tukey")))

# NLR - T30 - additional treatment and interaction
ED(ivg_dose_LN.3.additionaltreatment, 0.3, interval = "delta", vcov = vcovCL, type = "absolute")
NLR.add.T30 <- ED(ivg_dose_LN.3.additionaltreatment, 0.3, interval = "delta", vcov = vcovCL, multcomp = T, display = F, type="absolute")
summary(glht(NLR.add.T30[["EDmultcomp"]],
             linfct = contrMat(1:2, "Tukey")))


### FIGURE 3 - NLR and TTT plot

# CONFIDENCE BANDS FOR EACH CURVE

#### time 0 - additional treatment #### 
## TTE
coefVec_tte <- coef(ivg_time_LN.3.slope)
names(coefVec_tte) <- c("b","d1","d2","d3","d4","d5","d6","d7","d8","d9","e1","e2","e3","e4","e5","e6","e7","e8","e9")

predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec_tte, paste("d1*(pnorm(b*(log(",tival,") - log(e1))))"),
                         vcovCL(ivg_time_LN.3.slope)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_TTE_t0 <- t(predFctv(1:120))

## Predictions NLR
coefVec <- coeftest(ivg_dose_LN.3, vcov= sandwich)[,1]
names(coefVec) <- c("b1","b2","b3","b4","b5","b6","b7","b8","b9","d1","d2","d3","d4","d5","d6","d7","d8","d9","e1","e2","e3","e4","e5","e6","e7","e8","e9")

predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec, paste("d1*(pnorm(b1*(log(",tival,") - log(e1))))"),
                         vcovCL(ivg_dose_LN.3)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_NLR_T0 <- t(predFctv(1:120))

### T1P1 ###
### TTE
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec_tte, paste("d2*(pnorm(b*(log(",tival,") - log(e2))))"),
                         vcovCL(ivg_time_LN.3.slope)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_TTE_T1P1 <- t(predFctv(1:120))

## NLR
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec, paste("d2*(pnorm(b2*(log(",tival,") - log(e2))))"),
                         vcovCL(ivg_dose_LN.3)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_NLR_T1P1 <- t(predFctv(1:120))

#### T1P2  #### 
## TTE
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec_tte, paste("d3*(pnorm(b*(log(",tival,") - log(e3))))"),
                         vcovCL(ivg_time_LN.3.slope)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_TTE_T1P2 <- t(predFctv(1:120))

## NLR
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec, paste("d3*(pnorm(b3*(log(",tival,") - log(e3))))"),
                         vcovCL(ivg_dose_LN.3)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_NLR_T1P2 <- t(predFctv(1:120))

#### T2P1  #### 
## TTE
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec_tte, paste("d4*(pnorm(b*(log(",tival,") - log(e4))))"),
                         vcovCL(ivg_time_LN.3.slope)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_TTE_T2P1 <- t(predFctv(1:120))

## NLR
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec, paste("d4*(pnorm(b4*(log(",tival,") - log(e4))))"),
                         vcovCL(ivg_dose_LN.3)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_NLR_T2P1 <- t(predFctv(1:120))

#### T2P2  #### 
## TTE
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec_tte, paste("d5*(pnorm(b*(log(",tival,") - log(e5))))"),
                         vcovCL(ivg_time_LN.3.slope)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_TTE_T2P2 <- t(predFctv(1:120))

## NLR
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec, paste("d5*(pnorm(b5*(log(",tival,") - log(e5))))"),
                         vcovCL(ivg_dose_LN.3)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_NLR_T2P2 <- t(predFctv(1:120))

#### T3P1  #### 
## TTE

predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec_tte, paste("d6*(pnorm(b*(log(",tival,") - log(e6))))"),
                         vcovCL(ivg_time_LN.3.slope)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_TTE_T3P1 <- t(predFctv(1:120))

## NLR

predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec, paste("d6*(pnorm(b6*(log(",tival,") - log(e6))))"),
                         vcovCL(ivg_dose_LN.3)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_NLR_T3P1 <- t(predFctv(1:120))

#### T3P2  #### 
## TTE

predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec_tte, paste("d7*(pnorm(b*(log(",tival,") - log(e7))))"),
                         vcovCL(ivg_time_LN.3.slope)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_TTE_T3P2 <- t(predFctv(1:120))

## NLR
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec, paste("d7*(pnorm(b7*(log(",tival,") - log(e7))))"),
                         vcovCL(ivg_dose_LN.3)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_NLR_T3P2 <- t(predFctv(1:120))


#### T4P1  #### 
## TTE
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec_tte, paste("d8*(pnorm(b*(log(",tival,") - log(e8))))"),
                         vcovCL(ivg_time_LN.3.slope)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_TTE_T4P1 <- t(predFctv(1:120))

## NLR
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec, paste("d8*(pnorm(b8*(log(",tival,") - log(e8))))"),
                         vcovCL(ivg_dose_LN.3)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_NLR_T4P1 <- t(predFctv(1:120))

#### T4P2  #### 
## TTE
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec_tte, paste("d9*(pnorm(b*(log(",tival,") - log(e9))))"),
                         vcovCL(ivg_time_LN.3.slope)))
}
predFctv <- Vectorize(predFct, "tival")

etpred_TTE_T4P2 <- t(predFctv(1:120))


## NLR
predFct <- function(tival)
{
  as.numeric(deltaMethod(coefVec, paste("d9*(pnorm(b9*(log(",tival,") - log(e9))))"),
                         vcovCL(ivg_dose_LN.3)))
}
predFctv <- Vectorize(predFct, "tival")
etpred_NLR_T4P2 <- t(predFctv(1:120))

# plot all together
par(mfrow=c(3,3), mar=c(3.2,3.2,0.5,0.5), mgp=c(2,.7,0), xpd=T)

plot(ivg_time_LN.3.slope, level=c(1), type="average", lty=2, ylim=c(0,1), bty="L")
polygon(c(1:120, 120:1), c(etpred_TTE_t0[,1]-1.96*etpred_TTE_t0[,2], rev(etpred_TTE_t0[,1]+1.96*etpred_TTE_t0[,2])),
        col=rgb(0,0,1, alpha = 0.2), border = F)
plot(ivg_dose_LN.3, level="0meses:p0",type="average", add = T)
polygon(c(1:120, 120:1), c(etpred_NLR_T0[,1]-1.96*etpred_NLR_T0[,2], rev(etpred_NLR_T0[,1]+1.96*etpred_NLR_T0[,2])),
        col=rgb(0,0,0, alpha = 0.2), border = F)

plot(ivg_time_LN.3.slope, level=c(2), type="average", lty=2, ylim=c(0,1), bty="L")
polygon(c(1:120, 120:1), c(etpred_TTE_T1P1[,1]-1.96*etpred_TTE_T1P1[,2], rev(etpred_TTE_T1P1[,1]+1.96*etpred_TTE_T1P1[,2])),
        col=rgb(0,0,1, alpha = 0.2), border = F)
plot(ivg_dose_LN.3, level="1mes:p0,10",type="average", add = T)
polygon(c(1:120, 120:1), c(etpred_NLR_T1P1[,1]-1.96*etpred_NLR_T1P1[,2], rev(etpred_NLR_T1P1[,1]+1.96*etpred_NLR_T1P1[,2])),
        col=rgb(0,0,0, alpha = 0.2), border = F)

plot(ivg_time_LN.3.slope, level=c(3), type="average", lty=2, ylim=c(0,1), bty="L")
polygon(c(1:120, 120:1), c(etpred_TTE_T1P2[,1]-1.96*etpred_TTE_T1P2[,2], rev(etpred_TTE_T1P2[,1]+1.96*etpred_TTE_T1P2[,2])),
        col=rgb(0,0,1, alpha = 0.2), border = F)
plot(ivg_dose_LN.3, level="1mes:p0,20",type="average", add = T)
polygon(c(1:120, 120:1), c(etpred_NLR_T1P2[,1]-1.96*etpred_NLR_T1P2[,2], rev(etpred_NLR_T1P2[,1]+1.96*etpred_NLR_T1P2[,2])),
        col=rgb(0,0,0, alpha = 0.2), border = F)

plot(ivg_time_LN.3.slope, level=c(4), type="average", lty=2, ylim=c(0,1), bty="L")
polygon(c(1:120, 120:1), c(etpred_TTE_T2P1[,1]-1.96*etpred_TTE_T2P1[,2], rev(etpred_TTE_T2P1[,1]+1.96*etpred_TTE_T2P1[,2])),
        col=rgb(0,0,1, alpha = 0.2), border = F)
plot(ivg_dose_LN.3, level="2mes:p0,10",type="average", add = T)
polygon(c(1:120, 120:1), c(etpred_NLR_T2P1[,1]-1.96*etpred_NLR_T2P1[,2], rev(etpred_NLR_T2P1[,1]+1.96*etpred_NLR_T2P1[,2])),
        col=rgb(0,0,0, alpha = 0.2), border = F)

plot(ivg_time_LN.3.slope, level=c(5), type="average", lty=2, ylim=c(0,1), bty="L")
polygon(c(1:120, 120:1), c(etpred_TTE_T2P2[,1]-1.96*etpred_TTE_T2P2[,2], rev(etpred_TTE_T2P2[,1]+1.96*etpred_TTE_T2P2[,2])),
        col=rgb(0,0,1, alpha = 0.2), border = F)
plot(ivg_dose_LN.3, level="2mes:p0,20",type="average", add = T)
polygon(c(1:120, 120:1), c(etpred_NLR_T2P2[,1]-1.96*etpred_NLR_T2P2[,2], rev(etpred_NLR_T2P2[,1]+1.96*etpred_NLR_T2P2[,2])),
        col=rgb(0,0,0, alpha = 0.2), border = F)

plot(ivg_time_LN.3.slope, level=c(6), type="average", lty=2, ylim=c(0,1), bty="L")
polygon(c(1:120, 120:1), c(etpred_TTE_T3P1[,1]-1.96*etpred_TTE_T3P1[,2], rev(etpred_TTE_T3P1[,1]+1.96*etpred_TTE_T3P1[,2])),
        col=rgb(0,0,1, alpha = 0.2), border = F)
plot(ivg_dose_LN.3, level="3mes:p0,10",type="average", add = T)
polygon(c(1:120, 120:1), c(etpred_NLR_T3P1[,1]-1.96*etpred_NLR_T3P1[,2], rev(etpred_NLR_T3P1[,1]+1.96*etpred_NLR_T3P1[,2])),
        col=rgb(0,0,0, alpha = 0.2), border = F)

plot(ivg_time_LN.3.slope, level=c(7), type="average", lty=2, ylim=c(0,1), bty="L")
polygon(c(1:120, 120:1), c(etpred_TTE_T3P2[,1]-1.96*etpred_TTE_T3P2[,2], rev(etpred_TTE_T3P2[,1]+1.96*etpred_TTE_T3P2[,2])),
        col=rgb(0,0,1, alpha = 0.2), border = F)
plot(ivg_dose_LN.3, level="3mes:p0,20",type="average", add = T)
polygon(c(1:120, 120:1), c(etpred_NLR_T3P2[,1]-1.96*etpred_NLR_T3P2[,2], rev(etpred_NLR_T3P2[,1]+1.96*etpred_NLR_T3P2[,2])),
        col=rgb(0,0,0, alpha = 0.2), border = F)

plot(ivg_time_LN.3.slope, level=c(8), type="average", lty=2, ylim=c(0,1), bty="L")
polygon(c(1:120, 120:1), c(etpred_TTE_T4P1[,1]-1.96*etpred_TTE_T4P1[,2], rev(etpred_TTE_T4P1[,1]+1.96*etpred_TTE_T4P1[,2])),
        col=rgb(0,0,1, alpha = 0.2), border = F)
plot(ivg_dose_LN.3, level="4mes:p0,10",type="average", add = T)
polygon(c(1:120, 120:1), c(etpred_NLR_T4P1[,1]-1.96*etpred_NLR_T4P1[,2], rev(etpred_NLR_T4P1[,1]+1.96*etpred_NLR_T4P1[,2])),
        col=rgb(0,0,0, alpha = 0.2), border = F)

plot(ivg_time_LN.3.slope, level=c(9), type="average", lty=2, ylim=c(0,1), bty="L")
polygon(c(1:120, 120:1), c(etpred_TTE_T4P2[,1]-1.96*etpred_TTE_T4P2[,2], rev(etpred_TTE_T4P2[,1]+1.96*etpred_TTE_T4P2[,2])),
        col=rgb(0,0,1, alpha = 0.2), border = F)
plot(ivg_dose_LN.3, level="4mes:p0,20",type="average", add = T)
polygon(c(1:120, 120:1), c(etpred_NLR_T4P2[,1]-1.96*etpred_NLR_T4P2[,2], rev(etpred_NLR_T4P2[,1]+1.96*etpred_NLR_T4P2[,2])),
        col=rgb(0,0,0, alpha = 0.2), border = F)
legend(100, 0.5, legend=c("NLR"), cex = 1,
       col=c("black"), lty=c(1), box.lty = 0)
legend(50, 1.2, legend=c("TTE"),cex = 1.5,
       col=c("blue"), lty=c(2), box.lty = 4)

################################################################################
############################# Emergence Speed Index (ESI) ######################
################################################################################

# Adjusting different models
ESI.gamma <- glm(GSI ~ Time*Pack, family=Gamma, data=ESI)
ESI.gauss <- glm(GSI ~ Time*Pack, data = ESI)
ESI.invgauss <- glm(GSI ~ Time*Pack, family=inverse.gaussian, data=ESI)

AIC(ESI.gamma, ESI.gauss, ESI.invgauss)
# Gaussian provides lower AIC

#checking normality
plot(ESI.gauss) # ok
shapiro.test(residuals(ESI.gauss)) # ok

#checkin model reductions
ESI.gauss.2 <- glm(GSI ~ Time+Pack, data=ESI) #normal
ESI.gauss.3 <- glm(GSI ~ Time, data=ESI) #normal
ESI.gauss.4 <- glm(GSI ~ Pack, data=ESI) #normal
#reductions can't be made.
# As normal distribution fit better, ANOVA counting the additional treatment vs interaction will be evaluated by the package Exp.Des

# Separating the interaction from the additional treatment and creating repetition column
Add <- ESI[1:8,]
Interaction <- subset(ESI, Time!="0mes")
Interaction$rep <- seq(1,8, by=1)
# Anova
fat2.ad.dic(Interaction$Time, Interaction$Pack, repet = F, Interaction$GSI, Add$GSI, quali = c(TRUE, TRUE), mcomp = "tukey", 
            fac.names = c("tempo", "embalagam"), sigT = 0.05, sigF = 0.05)
# the additional treatment vs interaction non-significance 

#Table 3
# Using emmeans for contrasts between all treatmens
a.emm <-  emmeans(ESI.gauss, transform = "response", specs = ~ Time*Pack, data = ESI)
multcomp::cld (a.emm, reversed = TRUE, Letter = "abcdfghijk") #mostrar as letras do tukey


