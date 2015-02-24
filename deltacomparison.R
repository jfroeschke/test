### FishMod delta lognormal and deltaBRT model comparison
### JT Froeschke
### January 07 2015
### comment 2
install.packages("fishMod")
install.packages("bestglm")
install.packages("pscl")

library(bestglm)
library(fishMod)
library(pscl)
### General approach for deltaLN
### I don't see a way to do model selection so will try
### Exhaustive search on both submodels and fit these models
### into the the deltaLN function.



##test data
data("bioChemists", package = "pscl")
## article production by graduate students in biochemistry Ph.D. programs
## logit-negbin
##pscl hurdle method
bc.hurdle<- hurdle(art ~ ., data = bioChemists, dist = "negbin")

##get fits


## fishMod deltaLN

#run the delta lognormal GLM for Spanish Mackerel (SMK)
require(fishMod)
bc.deltaLN <- deltaLN(ln.form=art~fem + mar + kid5 + phd + ment,
                      binary.form=~fem + mar + kid5 + phd + ment,data = bioChemists)
sanity.check <- subset(bioChemists, bioChemists$art > 0)
lm.sc <- lm(log(art)~ -1 + fem + mar + kid5 + phd + ment, data=sanity.check)
summary(lm.sc)
summary(bc.deltaLN$lnMod)


summary(bc.deltaLN$binMod)
anova(bc.deltaLN$binMod,test='Chisq')
summary(bc.deltaLN$lnMod)
anova(bc.deltaLNlnMod)



## delta BRT methods
library(gbm)
library(seascaper)
source("http://dl.dropbox.com/s/7g0rug82bqpp0b1/brt.functions.R") 
set.seed(123)  # set seed for reproducability

DF2 <- gbm.prep(DATA=bioChemists,GBM.Y=1);head(DF2) 

library(dismo)
## Fit models
Delta.bin <-gbm.step(data=DF2,
                     gbm.x=c(2:6),
                     gbm.y=(ncol(DF2)-1),
                     family = "bernoulli",
                     tree.complexity=5,
                     learning.rate=0.001,
                     bag.fraction=0.5) 

Input.var.length = data.frame(Delta.bin[19])
Input.var.length <- length(Input.var.length[, 1])
model.simp <- gbm.simplify(Delta.bin, n.drops = 3)
Mod.remove.test <- model.simp$deviance.summary[1, 1] + model.simp$deviance.summary[1,2]
if (Mod.remove.test < 0) {## if statement
     Delta.bin <- gbm.step(data = DF2,
                           gbm.x = model.simp$pred.list[[3]], 
                           gbm.y = (ncol(DF2)-1), 
                           family = "bernoulli", 
                           tree.complexity = 5, 
                           learning.rate = .01, 
                           bag.fraction = .5, 
                           plot.main = TRUE)    
} ## end if statement

### Positive submodel
DF3 <- subset(DF2, bin.y==1)
Delta.pos <-gbm.step(data=DF3,
                     gbm.x=c(2:6),
                     gbm.y=(ncol(DF2)),
                     family = "gaussian",
                     tree.complexity=5,
                     learning.rate=0.001,
                     bag.fraction=0.5) 

Input.var.length = data.frame(Delta.pos[19])
Input.var.length <- length(Input.var.length[, 1])
model.simp <- gbm.simplify(Delta.pos, n.drops = 3)
Mod.remove.test <- model.simp$deviance.summary[1, 1] + model.simp$deviance.summary[1,2]
if (Mod.remove.test < 0) {## if statement
     Delta.pos <- gbm.step(data = DF3,
                           gbm.x = model.simp$pred.list[[3]], 
                           gbm.y=(ncol(DF2)), 
                           family = "gaussian", 
                           tree.complexity = 5, 
                           learning.rate = .001, 
                           bag.fraction = .5, 
                           plot.main = TRUE)    
} ## end if statement


## Section 9.1 Predition to binomial data
## bin Predict to train and test data
preds.bin <- gbm.predict.grids(Delta.bin, DF2,
             want.grids = F, sp.name = "preds") 
preds.all.bin <- data.frame(index=DF2$index,
                    binomial=preds.bin)

## Section 9.2 Predition to positive  data
preds.pos <- gbm.predict.grids(Delta.pos, DF3,
                want.grids = F, sp.name = "preds") #predict to training data
preds.pos2 <- gbm.predict.grids(Delta.pos, DF2,
               want.grids = F, sp.name = "preds") #predict to training data
preds.all.pos <- data.frame(index=DF3$index,
                            positive=c(preds.pos))

## Section 9.3  Compute Annual index of relative abundance
## Combine into index
preds.all.index <- merge(preds.all.bin, preds.all.pos, by="index", all.x=TRUE)
preds.all.index$unlog.positive<- exp(preds.all.index$positive)
preds.all.index$unlog.positive[is.na(preds.all.index$unlog.positive)] <- 1

## Section 9.3.1 Create index at points
## merge with original data
DF.output <- merge(DF2, preds.all.index,by="index", all.x=TRUE)
DF.output$delta <- DF.output$binomial * DF.output$unlog.positive

head(DF.output)


### Summary dataset
hurdle.prob <- predict(bc.hurdle, bioChemists, 
                        type=c("prob"))

hurdle.count <- predict(bc.hurdle, bioChemists, 
                       type=c("count"))

hurdle.zero <- predict(bc.hurdle, bioChemists, 
                        type=c("zero"))

hurdle.response <- predict(bc.hurdle, bioChemists, 
                       type=c("response"))

hurdle.df <- data.frame(index=1:nrow(bioChemists),
                        #hurdle.prob=hurdle.prob,
                        hurdle.count=hurdle.count,
                        hurdle.zero=hurdle.zero,
                        hurdle.response=hurdle.response,
                        hurdle.check=hurdle.count*hurdle.zero) ##combined is product of each

###This matches and makes sense.
###Try same for deltaLN

#run the delta lognormal GLM for Spanish Mackerel (SMK)
require(fishMod)
## Fit the model
bc.deltaLN <- deltaLN(ln.form=art~fem + mar + kid5 + phd + ment,
                      binary.form=~fem + mar + kid5 + phd + ment,data = bioChemists)

bc.deltaLN$fitted  ##fitted values from delta lognormal model
deltaLN.bin <- predict(bc.deltaLN$binMod, bioChemists, type=c("response"))
deltaLN.pos <- predict(bc.deltaLN$lnMod, bioChemists, type=c("response"))



f <- bc.deltaLN$lnMod
f2 <- data.frame(f$fitted.values, fexp=exp(f$fitted.values))
copy <- bioChemists

##source of confusion
confusion <- cbind(bioChemists, fitted=bc.deltaLN$fitted,
                   bin.pred=deltaLN.bin, 
                   leftover=bc.deltaLN$fitted-deltaLN.bin)


copy$index <- 1:nrow(copy)
sanity.check <- subset(copy, copy$art > 0)
sc.glm <- glm(as.numeric(copy$art>0) ~ 
                   fem + mar + kid5 + phd + ment,data = bioChemists, family="binomial")
sc.lm <- lm(log(art) ~ -1 + fem + mar + kid5 + phd + ment, data=sanity.check)

sanity.check.df.glm <- data.frame(index=copy$index,
                                  glm.fitted=sc.glm$fitted.values
                                  deltaLN=bc.deltaLN)

sanity.check.df.lm <- data.frame(index=sanity.check$index,
                                 lm.fitted=sc.lm$fitted.values)

sanity.check.merge <- merge(sanity.check.df.glm, sanity.check.df.lm,
                            by="index", all.x=TRUE)
sanity.check.merge$delta <- bc.deltaLN$fitted
sanity.check.merge$unlog <- exp(sanity.check.merge$lm.fitted)

summary(sc.glm)
summary(bc.deltaLN$binMod)
sanity.check <- subset(bioChemists, bioChemists$art > 0)
lm.sc <- lm(log(art)~ -1 + fem + mar + kid5 + phd + ment, data=sanity.check)
summary(lm.sc)
summary(bc.deltaLN$lnMod)
########################################end here








####compare linear model with brt log normal
sanity.check$logart <- log(sanity.check$art)
tmp1 <- lm(logart ~ -1 + fem + mar + kid5 + phd + ment, data=sanity.check)
summary(tmp1)
set.seed(123)
tmp2 <- gbm.step(data = sanity.check,
               gbm.x = c(2:6), 
               gbm.y=7, 
               family = "gaussian", 
               tree.complexity = 1, 
               learning.rate = .001, 
               bag.fraction = 0.5, 
               plot.main = TRUE)
preds.pos2 <- gbm.predict.grids(tmp2, sanity.check,
                                want.grids = F, sp.name = "preds") 

x <- data.frame(lm=exp(tmp1$fitted.values),
                brt=exp(tmp2$fitted),
                brt.pred=exp(preds.pos2),
                art=sanity.check$art)
plot(brt~lm, data=x);abline(0,1)
summary(lm(log(brt)~-1 + log(lm), data=x))
plot(brt~art, data=x)
plot(lm~art, data=x)

### try to manually create delta lognormal model




### Will use this exhaustive serach
#fit binomial model
# best.bin <- cbind(DF2[, c(2,4:7)], cpue=DF2[,10])
# bestAIC <- bestglm(best.bin, family=binomial, IC="AIC", TopModels=20)
# #get best model
# best.bin.mods <- bestAIC$BestModels
# #convert to table, rename and code in LaTeX
# best.with.year <- subset(best.bin.mods, Year==TRUE)
# glm.bin <- glm(cpue ~ Year + Month + Block, data=best.bin, family="binomial")
# AIC(glm.bin)

#table for report
# 
# #fit positive model
# best.pos <- cbind(DF3[, c(2,4:7)], cpue=DF3[,11])
# bestAICpos <- bestglm(best.pos, family=gaussian, IC="AIC", TopModels=20)
# #get best model
# best.bin.mods.pos <- bestAICpos$BestModels
# #convert to table, rename and code in LaTeX
# x.best.pos <- xtable(best.bin.mods.pos)
# 
# 
# best.with.year.pos <- subset(best.bin.mods.pos, Year==TRUE)
# glm.pos <- glm(cpue ~ Year + Month + Block + Longitude, data=best.pos, family="gaussian")
# AIC(glm.pos)
# ##model validation
