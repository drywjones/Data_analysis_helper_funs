# Variance inflation factor calculations are used to test if 
# multicollinearity is adversely impacting estimates for coefficients
# of multivariate linear models. There is some debate as to what VIF is 
# considered too high. One segestion is to eliminate anything greater than
# 3 (Zuur et al. 2010), while others have suggested values as high as 10 
# (Montgomery & Peck, 1992). A VIF value of 1.9 indicates that the variance 
# of a given coefficient is 90% higher than would be expected if there was 
# no multicollinearity. The below function will calculate VIF values for all 
# independent variables (fixed only) in a linear or mixed effects linear model.
#
# A more complete description of VIF can be found in this reference:
#
# Zuur, A.F.; Ieno, E.N.; Elphick, C.S. A protocol for data exploration to avoid common statistical problems.
# Methods Ecol. Evol. 2010, 1, 3-14.

library(nlme)
library(lme4)

## This function takes a linear model (lm, lmer, lme) and calculates
## the VIF value for the individual fixed effects.
vif.fun <- function (fit) {
  ## adapted from rms::vif
  v <- vcov(fit)
  ## lm models don't have a fixef option so have to use
  ## a different names option for those models
  if(class(fit)[1]=="lm"){
    nam<- names(coef(fit))
  }else{
    nam <- names(fixef(fit))
  }

  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)] }
  ## calculate VIFs
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  ## name output values with fixed effects names
  names(v) <- nam
  v }

## Create data for nested random effects: #####
test.data<-data.frame(indep.var=sample(1:100,5000,replace=T)/10)
test.data$dep.var<-test.data$indep.var*(1+sample(1:100,5000,replace=T)/200)
test.data$indep.var2<-test.data$indep.var+sample(1:100,5000,replace=T)/10
test.data$indep.var3<-test.data$indep.var*sample(1:100,5000,replace=T)/10
test.data$indep.var4<-test.data$indep.var+sample(1:100,5000,replace=T)/60
test.data$SITE<-c(rep("SITE1",1000),rep("SITE2",1000),rep("SITE3",1000),rep("SITE4",1000),rep("SITE5",1000))
test.data$plot<-c(rep(c(rep("plot1",250),rep("plot2",250),rep("plot3",250),rep("plot4",250)),5))
test.data$dep.var<-test.data$dep.var+
  ifelse(test.data$SITE=="SITE1",(sample(1:100,1)/100+4),
         ifelse(test.data$SITE=="SITE2",(sample(1:100,1)/100+3),
                ifelse(test.data$SITE=="SITE3",(sample(1:100,1)/100+9),
                       ifelse(test.data$SITE=="SITE4",(sample(1:100,1)/100+6),
                              (sample(1:100,1)/100+7)))))
test.data$dep.var<-test.data$dep.var+
  ifelse(test.data$plot=="plot1",(sample(1:100,1)/100+4),
         ifelse(test.data$SITE=="plot2",(sample(1:100,1)/100+3),
                ifelse(test.data$SITE=="plot3",(sample(1:100,1)/100+5),
                       (sample(1:100,1)/100+2))))

## Fit three types of linear models to data: ####
lm.fit<-lm(dep.var~indep.var+indep.var2+indep.var3,
           data=test.data)

lme.fit<-lme(dep.var~indep.var+indep.var2+indep.var3,
             random=~1|SITE/plot,
             data=test.data)
lmer.fit<-lmer(dep.var~indep.var+indep.var2+indep.var3+(1|SITE/plot),
               data=test.data)

## Get VIF values for three model fits:
vif.lm<-vif.fun(lm.fit)
vif.lme<-vif.fun(lme.fit)
vif.lmer<-vif.fun(lmer.fit)

## Check the VIF values
vif.lm
vif.lme
vif.lmer
