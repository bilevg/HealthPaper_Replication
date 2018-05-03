## Replication R code for Cammett, Melani, Julia Lynch, and Gavril Bilev. "The Influence of Private Health Care Financing on Citizen Trust in Government." (2015).
## Code by Gavril Bilev (bilevg@merrimack.edu)
library(lme4)
library(ez)
library(ggplot2)
library(scales)
library(xtable)
library(Amelia)
library(plyr)
library(apsrtable)
library(ltm)
library(sandwich)
library(lmtest)
## be sure to set the working dir to the script dir with setwd('/path/to/script/')
source(file="customapsr.R") ## my mod to do latex for "mer"; new function is "apsr" which takes as its argument a list of "mer" class models objects
options(digits=4, max.print=1000, scipen=5)
load("ESS2008.zip")
attach(working)
## impute missing data for the variables we end up using
dataset <- data.frame(trstgov, stfgov, stfeco , voted.pgc ,
                      GDPpc.log , agea, agea.sq,
                      hincfel.r,  meded.dummy , polintr.r , ppltrst , happy ,
                      dscrgrp.r, gincdif, dfincac, gvhlthc , dcndleq, sick,
                      PocketFin, cntry.r, Postcom.state, Nordic.state, stfhlth,
                      PrivFin, PrivFinWho, PrivInsFin,
                      close.pgc, lrscale, lrscale.sq = lrscale^2, lknhlcn,
                      hlthhmp, hamper.dummy,
                      hinctnta, pro.state,  leftist, Gini, THEcap, CofCoruption,
                      P90P10)
## scaling function for lme4 models
rescale.ess <- function(x) {
    dat <- x@frame
    scaled.frommodel <- numcolwise(scale)(dat)
    scaled.frommodel$cntry.r <- dat$cntry.r
    out <- lmer(eval(formula(x)), data=scaled.frommodel)
    out
}
## robust SEs
robust.se <- function(model, cluster){
    M <- length(unique(cluster))
    N <- length(cluster)
    K <- model$rank
    dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
    uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
    rcse.cov <- dfc * sandwich(model, meat = crossprod(uj)/N)
    rcse.se <- coeftest(model, rcse.cov)
    return(rcse.se) ## modified from above - just return t-tests + p-values
}
## rescale for OLS
rescale.ols <- function(x) {
    dat <- x$model
    require(plyr)
    scaled.frommodel <- numcolwise(scale)(dat)
    scaled.frommodel$cntry.r <- dat$cntry.r
    out <- lm(eval(formula(x)), data=scaled.frommodel)
    out
}
## function to impute by country then combine randomly selecting within country
#
## combine them into 1 dataset by randoml sampling 1/m*n from each BY COUNTRY and combine into one, then run lmer on it. using m=10.
impute.ess <- function(data, cs, m){
    require(Amelia)
    mi <- amelia(data, cs= cs, m=m)
    data.mi <- data.frame(NULL)
    m <- length(mi$imputations)
    for(i in 1:m){
        data.mi <- rbind(data.mi,
                         do.call("rbind", lapply(split(mi$imputations[[i]],
                                                       (mi$imputations[[i]])[ ,paste(cs)]),
                                                 function(x) x[sample(nrow(x), nrow(x)/m, replace=FALSE),])))
    }
    data.mi
}


###################################################################################################################
## desriptives:
## mean by country
means <- aggregate(dataset[ ,-c(grep("cntry.r", names(dataset)))], list(cntry=cntry.r), FUN = "mean", na.rm=T, na.action=na.omit)
## ## health.r descriptive combines "Bad" (4) and "Very Bad" (5) health as percentage for each country
list.sick <- tapply(dataset$sick, cntry.r, table)
sick.r <- unlist(lapply(list.sick, function(x) round((x[4] + x[5])/sum(x[1:5], na.rm=T)*100, 0)))
## table for descriptives
mat <- data.frame(MeanTrust=round(means$trstgov,1), HCS=round(means$stfhlth,1),
                  Risk=round(means$lknhlcn,1), PrivateFinWHO=round(means$PrivFinWho,0),
                  Income=round(means$hincfel.r,1), Sick=sick.r)
rownames(mat) <- means$cntry
mat <- mat[with(mat, order(-MeanTrust)), ]
mat$PrivateFinWHO <- paste(mat$PrivateFinWHO, "%", sep="")
mat$Sick <- paste(mat$Sick, "%", sep="")
mat

# violin plot of the distribution of Trust
## function for taking the mean with confidence interval
stat_sum_df <- function(fun, geom="crossbar", ...) {
    stat_summary(fun.data=fun, colour="black", geom=geom, width=0.2, ...)
}
## reorder according to trust
data.reordered <- within(dataset,
                  cntry.r.2 <- reorder(cntry.r, trstgov, fun=mean, order=T, na.rm=T ))
## plot of means
dev.new(width=5, height=6)
plot1 <- ggplot(data.reordered, aes(y=trstgov, x=cntry.r.2)) + geom_violin(fill="cyan") +
    stat_sum_df(aes(group = "cntry.r.2"),fun="mean_cl_normal", geom="pointrange", shape=16, size=2/3) +
        coord_flip()  + ylab("Trust") + xlab(NULL) + ylim(0,10) + labs(title="Trust")
ggsave(file="Violin plot of mean trust.png", dpi=600, width=7, height=7, units=c("in"))
plot1

## Violion plot of the distribution of government's responsibility to take care of the sick
## GVHLTHC
data.reordered <- within(dataset,
                  cntry.r.2 <- reorder(cntry.r, gvhlthc, fun=mean, order=T, na.rm=T ))
## plot of means
dev.new(width=5, height=6)
plot1 <- ggplot(data.reordered, aes(y=gvhlthc, x=cntry.r.2)) + geom_violin(fill="cyan") +
    stat_sum_df(aes(group = "cntry.r.2"),fun="mean_cl_normal", geom="pointrange", shape=16, size=2/3) +
        coord_flip()  + ylab("Responsibility") + xlab(NULL) + ylim(0,10) + labs(title="Violin plot of gov't\nresponsibility to take care of sick")
ggsave(file="Violin plot of gov't responsibility to take care of sick.png", dpi=600, width=7, height=7, units=c("in"))
plot1

## Stacked bar plot for Income
## data preparation
df.income <- dataset[ ,c("cntry.r", "hincfel.r")]
df.income <- within(df.income,
                  cntry.r.2 <- reorder(cntry.r, hincfel.r, fun=mean, order=T, na.rm=T ))
df.income$type <- factor(df.income$hincfel.r, labels=c("Finding it\nvery difficult", "Difficult", "Coping", "Living\ncomfortably"), levels=c(1:4))
## reorder
df.income <- plyr::arrange(df.income, cntry.r.2, type)
## calculate percentages for each category
df.income <- ddply(df.income, .(cntry.r.2, type), function(x)
                        count=nrow(x))
## vector for labeling the percent in "Very difficult"
poor.2 <- ddply(df.income, .(cntry.r.2), plyr::summarize,
                    poor.2=round((V1)[1]/sum(V1)*100, 0))[ ,2]
## the plot with geom_bar and position="fill"
ggplot(na.omit(df.income), aes(x=cntry.r.2, fill=type)) +
    geom_bar(aes(weight=V1, fill = type), position = 'fill', width=.7) +
        annotate("text", x=1:25, y = -0.05,
                 label = paste(poor.2, "%", sep=""), size=3.5) +
                     coord_flip() + ylab("") + xlab(NULL) +
                         scale_y_continuous(labels=percent) +
    ## guides(fill = guide_legend(label.position = "bottom")) +
    labs(title="Income with Percent of Respondents\n in \"Very Difficult\" Category") + scale_fill_grey(name="Living on \nPresent Income") + theme_minimal()
ggsave(file="Stacked Bar Plot of Income.png", dpi=600, width=7, height=7, units=c("in"))


#####################################################################################
## plot of means for financing and risk by income category
## no missing in the wrapping factor
df.trans <- dataset[!is.na(dataset$hincfel.r),]
df.trans <- ddply(df.trans, .(cntry.r, hincfel.r), transform, Financing = mean(PrivFinWho/100))
df.trans <- ddply(df.trans, .(cntry.r, hincfel.r), transform, Risk = mean(lknhlcn, na.rm=T))
df.trans <- ddply(df.trans, .(cntry.r, hincfel.r), transform, Count = length(Risk))
df.trans$Income <- factor(df.trans$hincfel.r, levels=c(1:4), labels=c("Finding it\n very difficult", "Difficult", "Coping", "Living\n comfortably"))
## Label Switzerland:
df.swiss <- data.frame(Income = levels(df.trans$Income), Risk = unique(df.trans$Risk[df.trans$cntry.r == "Switzerland"]))
## ## and Belgium
## df.belg <- data.frame(Income = levels(df.trans$Income), Risk = unique(df.trans$Risk[df.trans$cntry.r == "Belgium"]))
## and Bulgaria
df.bg <- data.frame(Income = levels(df.trans$Income), Risk = unique(df.trans$Risk[df.trans$cntry.r == "Bulgaria"]))
## and Ukraine
df.ur <- data.frame(Income = levels(df.trans$Income), Risk = unique(df.trans$Risk[df.trans$cntry.r == "Ukraine"]))

## Figure 5
ggplot(data=df.trans, aes(x=Financing, y=Risk)) +
    geom_point(aes(size=Count), shape=1) +
    geom_text(data=df.swiss, size=3.5, aes(y=Risk + .085, x=.4088), label="Switzerland", hjust=.8) +
    geom_point(data=df.swiss, size=5, aes(y=Risk, x=.4088), shape=0) +
    geom_text(data=df.ur, size=3.5, aes(y=Risk + .085, x=.4326), label="Ukraine", hjust=.8) +
    geom_point(data=df.ur, size=5, aes(y=Risk, x=.4326), shape=0) +
    geom_text(data=df.bg, size=3.5, aes(y=Risk + .085, x=.4306), label="Bulgaria", hjust=.9) +
    geom_point(data=df.bg, size=5, aes(y=Risk, x=.4306), shape=0) +
    facet_wrap( ~ Income, ncol=4) +
    xlab("Private Financing") + geom_smooth(method="lm", aes(x=PrivFinWho/100, y=lknhlcn)) +
    ylab("Mean Risk: How likely NOT to receive \nhealthcare in next 12 months if needed")  + ggtitle("Financing and Mean Risk for Different Income Categories")  +
    scale_x_continuous(label=percent, limits=c(.1, .45)) +
    scale_y_continuous(labels=c("Not at\n all likely", "Not very\n likely", "Likely", "Very\n Likely"), breaks=c(1:4),  limits=c(1,3.9)) + theme(legend.position=c(.88,.85)) +
    scale_size(breaks=c(100, 200, 500, 1000), name="Respondents in income\ncategory for country")
ggsave(file="Financing and Mean Risk per Income Category.png", dpi=600, width=7, height=7, units=c("in"))

############################################################################################################################################################################################################################################################################################################

## make a summary descriptive table: mean, min, max, sd
mean.d <- sapply(dataset[ , -c(grep("cntry.r", names(dataset)))], mean, na.rm=T)
min.d <- sapply(dataset[ , -c(grep("cntry.r", names(dataset)))], min, na.rm=T)
max.d <- sapply(dataset[ , -c(grep("cntry.r", names(dataset)))], max, na.rm=T)
sd.d <- sapply(dataset[ , -c(grep("cntry.r", names(dataset)))], sd, na.rm=T)
descr <- data.frame(MEAN=mean.d, MIN=min.d, MAX=max.d, SD=sd.d)
descr <- round(descr, 2)
write.csv(descr, "descr.csv")


################################################################################################################################################################################################################################################################################################################################
# OLS with Robust Standard Errors
## Model 1 with Robust SEs
ols.trst <- trstgov ~  ppltrst + happy + voted.pgc +  stfeco + polintr.r +
    dscrgrp.r +
    hincfel.r  +  agea + agea.sq + sick + meded.dummy +  pro.state
model1.RSE <- lm(ols.trst, data=dataset)
## ## First rescale:
## scaled.data <- numcolwise(scale)(dataset)
## scaled.data$cntry.r <- dataset$cntry.r
## ## run the model
## model1.RSE <- lm(ols.trst, data=scaled.data)
## Then RSE
index <- cntry.r[-(model1.RSE$na.action)] # need this because of missing cases
rse <- robust.se(model1.RSE, cluster=index) # get vector of robust SEs
##
model1.RSE$se <- rse[,2] # 2nd column is robust SEs; pass this to apsrtable w/ se="robust"
textols <- apsrtable(model1.RSE, stars="default", digits=3, se="robust",
          order="longest")
textols <- rbind(textols[3],  textols[5])

############################################################################################################################################################################################################################################################################################################
## Main Models
## Model 2
formula.trstgov2 <- trstgov ~ PrivFinWho +  Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + pro.state + (1 | cntry.r)
model2 <- lmer(formula.trstgov2, data=dataset)
model2.r <-  rescale.ess(model2)
## Model 3
formula.lknhlcn3 <- lknhlcn ~ PrivFinWho + Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + PrivFinWho*hincfel.r +
    (1 | cntry.r)
model3 <- lmer(formula.lknhlcn3, data=dataset)
model3.r <- rescale.ess(model3)
## Model 4
formula.stfhlth4 <- stfhlth ~ PrivFinWho + Nordic.state + Postcom.state +
               ppltrst + happy + voted.pgc + stfeco + polintr.r +  dscrgrp.r + hincfel.r + agea + agea.sq + sick + lknhlcn +  meded.dummy + PrivFinWho*hincfel.r + PrivFinWho*lknhlcn + (1 | cntry.r)
model4 <- lmer(formula.stfhlth4, data=dataset)
model4.r <-  rescale.ess(model4)
## Model 5
formula.trstgov5 <- trstgov ~ PrivFinWho + Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  + agea + agea.sq + sick + meded.dummy  + pro.state + lknhlcn + stfhlth + (1 | cntry.r)
model5 <- lmer(formula.trstgov5, data=dataset)
model5.r <-  rescale.ess(model5)
## models
models <- list(model2.r, model3.r, model4.r, model5.r)
dv.names <- sapply(models,function(x) all.vars(terms(x))[1])
textm <- apsr(models, stars="default", digits=3,
              model.names=paste("model", c(2:5),": ",dv.names, sep=""),
                order="longest", #float="sidewaystable",
                tsize=ifelse(length(models)<3, 1/4, 1), #label=IV.vars[[i]],
                ##caption=paste("Models with", int.cntr)
)


## financing for each country
ddply(dataset, .(cntry.r), summarise, financing = mean(PrivFinWho))

############################################################################################################################################################################################################################################################################################################
## Imputation robustness check
data.mi <- impute.ess(dataset, "cntry.r", 5)

## Model 2
model2 <- lmer(formula.trstgov2, data=data.mi)
model2.r <-  rescale.ess(model2)
## Model 3
model3 <- lmer(formula.lknhlcn3, data=data.mi)
model3.r <- rescale.ess(model3)
## Model 4
model4 <- lmer(formula.stfhlth4, data=data.mi)
model4.r <-  rescale.ess(model4)
## Model 5
model5 <- lmer(formula.trstgov5, data=data.mi)
model5.r <-  rescale.ess(model5)
## models
models <- list(model2.r, model3.r, model4.r, model5.r)
dv.names <- sapply(models,function(x) all.vars(terms(x))[1])
textm <- apsr(models, stars="default", digits=3,
              model.names=paste("model", c(2:5),": ",dv.names, sep=""),
                order="longest", #float="sidewaystable",
                tsize=ifelse(length(models)<3, 1/4, 1), #label=IV.vars[[i]],
                ##caption=paste("Models with", int.cntr)
)


############################################################################################################################################################################################################################################################################################################
## Another robustness check - exclude the Voted for Party in Governing Coalition variable to see the effect of the missingness
## Model 2
formula.trstgov2 <- trstgov ~ PrivFinWho +  Nordic.state + Postcom.state +
    ppltrst + happy +    stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + pro.state + (1 | cntry.r)
model2 <- lmer(formula.trstgov2, data=dataset)
model2.r <-  rescale.ess(model2)
## Model 3
formula.lknhlcn3 <- lknhlcn ~ PrivFinWho + Nordic.state + Postcom.state +
    ppltrst + happy +    stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + PrivFinWho*hincfel.r +
    (1 | cntry.r)
model3 <- lmer(formula.lknhlcn3, data=dataset)
model3.r <- rescale.ess(model3)
## Model 4
formula.stfhlth4 <- stfhlth ~ PrivFinWho + Nordic.state + Postcom.state +
               ppltrst + happy +   stfeco + polintr.r +  dscrgrp.r + hincfel.r + agea + agea.sq + sick + lknhlcn +  meded.dummy + PrivFinWho*hincfel.r + PrivFinWho*lknhlcn + (1 | cntry.r)
model4 <- lmer(formula.stfhlth4, data=dataset)
model4.r <-  rescale.ess(model4)
## Model 5
formula.trstgov5 <- trstgov ~ PrivFinWho + Nordic.state + Postcom.state +
    ppltrst + happy +    stfeco + polintr.r +  dscrgrp.r + hincfel.r  + agea + agea.sq + sick + meded.dummy  + pro.state + lknhlcn + stfhlth + (1 | cntry.r)
model5 <- lmer(formula.trstgov5, data=dataset)
model5.r <-  rescale.ess(model5)
## models
models <- list(model2.r, model3.r, model4.r, model5.r)
dv.names <- sapply(models,function(x) all.vars(terms(x))[1])
textm <- apsr(models, stars="default", digits=3,
              model.names=paste("model", c(2:5),": ",dv.names, sep=""),
                order="longest", #float="sidewaystable",
                tsize=ifelse(length(models)<3, 1/4, 1), #label=IV.vars[[i]],
                ##caption=paste("Models with", int.cntr)
)

## Another robustness check - exclude Switzerland
data.limited <- subset(dataset, cntry.r != "Switzerland")
## Model 2
formula.trstgov2 <- trstgov ~ PrivFinWho +  Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + pro.state + (1 | cntry.r)
model2 <- lmer(formula.trstgov2, data=data.limited)
model2.r <-  rescale.ess(model2)
## Model 3
formula.lknhlcn3 <- lknhlcn ~ PrivFinWho + Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + PrivFinWho*hincfel.r +
    (1 | cntry.r)
model3 <- lmer(formula.lknhlcn3, data=data.limited)
model3.r <- rescale.ess(model3)
## Model 4
formula.stfhlth4 <- stfhlth ~ PrivFinWho + Nordic.state + Postcom.state +
               ppltrst + happy + voted.pgc + stfeco + polintr.r +  dscrgrp.r + hincfel.r + agea + agea.sq + sick + lknhlcn +  meded.dummy + PrivFinWho*hincfel.r + PrivFinWho*lknhlcn + (1 | cntry.r)
model4 <- lmer(formula.stfhlth4, data=data.limited)
model4.r <-  rescale.ess(model4)
## Model 5
formula.trstgov5 <- trstgov ~ PrivFinWho + Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  + agea + agea.sq + sick + meded.dummy  + pro.state + lknhlcn + stfhlth + (1 | cntry.r)
model5 <- lmer(formula.trstgov5, data=data.limited)
model5.r <-  rescale.ess(model5)
## models
models <- list(model2.r, model3.r, model4.r, model5.r)
dv.names <- sapply(models,function(x) all.vars(terms(x))[1])
textm <- apsr(models, stars="default", digits=3,
              model.names=paste("model", c(2:5),": ",dv.names, sep=""),
                order="longest", #float="sidewaystable",
                tsize=ifelse(length(models)<3, 1/4, 1), #label=IV.vars[[i]],
                ##caption=paste("Models with", int.cntr)
)


## Robustness check - use Pocket Financing
## Model 2
formula.trstgov2 <- trstgov ~ PocketFin +  Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + pro.state + (1 | cntry.r)
model2 <- lmer(formula.trstgov2, data=dataset)
model2.r <-  rescale.ess(model2)
## Model 3
formula.lknhlcn3 <- lknhlcn ~ PocketFin + Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + PocketFin*hincfel.r +
    (1 | cntry.r)
model3 <- lmer(formula.lknhlcn3, data=dataset)
model3.r <- rescale.ess(model3)
## Model 4
formula.stfhlth4 <- stfhlth ~ PocketFin + Nordic.state + Postcom.state +
               ppltrst + happy + voted.pgc + stfeco + polintr.r +  dscrgrp.r + hincfel.r + agea + agea.sq + sick + lknhlcn +  meded.dummy + PocketFin*hincfel.r + PocketFin*lknhlcn + (1 | cntry.r)
model4 <- lmer(formula.stfhlth4, data=dataset)
model4.r <-  rescale.ess(model4)
## Model 5
formula.trstgov5 <- trstgov ~ PocketFin + Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  + agea + agea.sq + sick + meded.dummy  + pro.state + lknhlcn + stfhlth + (1 | cntry.r)
model5 <- lmer(formula.trstgov5, data=dataset)
model5.r <-  rescale.ess(model5)
## models
models <- list(model2.r, model3.r, model4.r, model5.r)
dv.names <- sapply(models,function(x) all.vars(terms(x))[1])
textm <- apsr(models, stars="default", digits=3,
              model.names=paste("model", c(2:5),": ",dv.names, sep=""),
                order="longest", #float="sidewaystable",
                tsize=ifelse(length(models)<3, 1/4, 1), #label=IV.vars[[i]],
                ##caption=paste("Models with", int.cntr)
)

## Another robustness check - exclude Bulgaria and Ukraine
data.limited <- subset(data.mi, cntry.r != "Bulgaria" & cntry.r != "Ukraine" )
## Model 2
formula.trstgov2 <- trstgov ~ PrivFinWho +  Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + pro.state + (1 | cntry.r)
model2 <- lmer(formula.trstgov2, data=data.limited)
model2.r <-  rescale.ess(model2)
## Model 3
formula.lknhlcn3 <- lknhlcn ~ PrivFinWho + Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + PrivFinWho*hincfel.r +
    (1 | cntry.r)
model3 <- lmer(formula.lknhlcn3, data=data.limited)
model3.r <- rescale.ess(model3)
## Model 4
formula.stfhlth4 <- stfhlth ~ PrivFinWho + Nordic.state + Postcom.state +
               ppltrst + happy + voted.pgc + stfeco + polintr.r +  dscrgrp.r + hincfel.r + agea + agea.sq + sick + lknhlcn +  meded.dummy + PrivFinWho*hincfel.r + PrivFinWho*lknhlcn + (1 | cntry.r)
model4 <- lmer(formula.stfhlth4, data=data.limited)
model4.r <-  rescale.ess(model4)
## Model 5
formula.trstgov5 <- trstgov ~ PrivFinWho + Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  + agea + agea.sq + sick + meded.dummy  + pro.state + lknhlcn + stfhlth + (1 | cntry.r)
model5 <- lmer(formula.trstgov5, data=data.limited)
model5.r <-  rescale.ess(model5)
## models
models <- list(model2.r, model3.r, model4.r, model5.r)
dv.names <- sapply(models,function(x) all.vars(terms(x))[1])
textm <- apsr(models, stars="default", digits=3,
              model.names=paste("model", c(2:5),": ",dv.names, sep=""),
                order="longest", #float="sidewaystable",
                tsize=ifelse(length(models)<3, 1/4, 1), #label=IV.vars[[i]],
                ##caption=paste("Models with", int.cntr)
)

## Unscaled models
## Main Models
## Model 2
formula.trstgov2 <- trstgov ~ PrivFinWho +  Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + pro.state + (1 | cntry.r)
model2 <- lmer(formula.trstgov2, data=dataset)
## Model 3
formula.lknhlcn3 <- lknhlcn ~ PrivFinWho + Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + PrivFinWho*hincfel.r +
    (1 | cntry.r)
model3 <- lmer(formula.lknhlcn3, data=dataset)
## Model 4
formula.stfhlth4 <- stfhlth ~ PrivFinWho + Nordic.state + Postcom.state +
               ppltrst + happy + voted.pgc + stfeco + polintr.r +  dscrgrp.r + hincfel.r + agea + agea.sq + sick + lknhlcn +  meded.dummy + PrivFinWho*hincfel.r + PrivFinWho*lknhlcn + (1 | cntry.r)
model4 <- lmer(formula.stfhlth4, data=dataset)
## Model 5
formula.trstgov5 <- trstgov ~ PrivFinWho + Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  + agea + agea.sq + sick + meded.dummy  + pro.state + lknhlcn + stfhlth + (1 | cntry.r)
model5 <- lmer(formula.trstgov5, data=dataset)
## models
models <- list(model2, model3, model4, model5)
dv.names <- sapply(models,function(x) all.vars(terms(x))[1])
textm <- apsr(models, stars="default", digits=3,
              model.names=paste("model", c(2:5),": ",dv.names, sep=""),
                order="longest", #float="sidewaystable",
                tsize=ifelse(length(models)<3, 1/4, 1), #label=IV.vars[[i]],
                ##caption=paste("Models with", int.cntr)
)


## short illustration of 2 different regression lines for public and private financing, response variable is risk, based on model 3, unscaled
ggplot() + xlim(0,5) + ylim(0,5) + geom_abline(intercept=2.355, slope=-.127,col="blue") + geom_abline(intercept=2.955, slope=-2.087, col="red") + labs(xlab="Income", ylab="Risk")

## Simulations and figures 3-4
## Model 3; figure 3
formula.lknhlcn3 <- lknhlcn ~ PrivFinWho + Nordic.state + Postcom.state +
    ppltrst + happy + voted.pgc +  stfeco + polintr.r +  dscrgrp.r + hincfel.r  +
    agea + agea.sq +  sick + meded.dummy  + PrivFinWho*hincfel.r +
    (1 | cntry.r)
model3 <- lmer(formula.lknhlcn3, data=dataset)
## simulations with package ez
range <- expand.grid(PrivFinWho=seq(5,45,5), hincfel.r=c(1,4))
## means
data <- model3@frame
cn <- names(data) %in% c("cntry.r", "PrivFinWho", "hincfel.r")
means <-  colMeans(data[!cn])
means <- data.frame(t(means))
means.df <- means[rep(1:nrow(means), nrow(range)), ]
to_predict <- data.frame(cbind(range,means.df))
preds <- ezPredict(fit = model3, to_predict=to_predict, iterations=1000)
##
myplot <- ezPlot2(preds=preds,x=PrivFinWho, split=hincfel.r, do_lines=F,
        levels=list(hincfel.r=list(new_names=c("Low Income","High Income"))),
        y_lab="Risk", split_lab="Income", bar_width=3, x_lab="Private Financing")
myplot + scale_color_manual(values=c( "gray50", "black")) + ylim(1,3) +
        scale_shape_manual(values=c(17, 16))
## ggtitle("Figure 3: Risk by Income and Financing (Model 3)")
ggsave(file="Risk by Income and Financing (Model 3).png", dpi=600, width=7, height=5, units=c("in"))


## Model 4; figure 4
formula.stfhlth4 <- stfhlth ~ PrivFinWho + Nordic.state + Postcom.state +
               ppltrst + happy + voted.pgc + stfeco + polintr.r +  dscrgrp.r + hincfel.r + agea + agea.sq + sick + lknhlcn +  meded.dummy + PrivFinWho*hincfel.r + PrivFinWho*lknhlcn + (1 | cntry.r)
model4 <- lmer(formula.stfhlth4, data=dataset)
## simulations with package ez
range <- expand.grid(PrivFinWho=seq(5,45,5), lknhlcn=c(1,4))
## means
data <- model4@frame
cn <- names(data) %in% c("cntry.r", "stfhlth", "PrivFinWho", "lknhlcn")
means <-  colMeans(data[!cn])
means <- data.frame(t(means))
means.df <- means[rep(1:nrow(means), nrow(range)), ]
to_predict <- data.frame(cbind(range,means.df))
preds <- ezPredict(fit = model4, to_predict=to_predict, iterations=1000)
myplot2 <- ezPlot2(preds=preds,x=PrivFinWho, split=lknhlcn,
                   do_lines=F, levels=list(lknhlcn=list(new_names=c("Low Risk","High Risk"))),
        y_lab="Health System Satisfaction", split_lab="Risk", x_lab="Private Financing", do_plot=T)
myplot2 + scale_color_manual(values=c( "black", "gray50")) + ylim(1,8) +
    scale_shape_manual(values=c(16, 17)) +
    ggtitle("")
ggsave(file="Health Satisfaction (Model 4) by Risk and Financing.png", dpi=600,
       width=7, height=5)

## robustness check, model 4 without Risk
## Model 4
formula.stfhlth4.norisk <- stfhlth ~ PrivFinWho + Nordic.state +
    Postcom.state + ppltrst + happy + voted.pgc + stfeco +
    polintr.r +  dscrgrp.r + hincfel.r + agea + agea.sq +
    sick +  meded.dummy + PrivFinWho*hincfel.r + (1 | cntry.r)
model4.norisk <- lmer(formula.stfhlth4.norisk, data=dataset)
model4.r.norisk <-  rescale.ess(model4.norisk)

## END
