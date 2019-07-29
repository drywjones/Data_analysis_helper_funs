# Marginal and conditional R2 values for nlme and lme4 objects
# R2 calculations based on suggested approaches in:
# Nakagawa, S., and Holger, S. 2013. A general and simple method for obtaining R2 from generalized
# linear mixed-effects models. Methods in Ecology and Evolution.

# These libraries are needed for fitting and graphing but 
# the function works independently from them.
library(nlme)
library(lme4)
library(ggplot2)

## Create data for nested random effects: #####
test.data<-data.frame(indep.var=sample(1:100,5000,replace=T)/10)
test.data$dep.var<-test.data$indep.var*(1+sample(1:100,5000,replace=T)/200)
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
                              
### functions ####

## This function is used within the r2.for.mm function to reduce the R2 value
## to the desired number of decimal places but I put it here as well in case 
## someone finds it useful on it's own.  

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

## The r2.for.mm function takes a mixed model object created in lme4 or nlme
## and creates labels for conditional (Rc) and marginal (Rm) R2 values. Three 
## formats are created, unformated text, text formated to display formatted 
## text for base R graphic plots, and text formatted for ggplot plots. 
## Additionally R2 labels are limited for lme4 to only full model (with random 
## effects) conditional R2 values, and marginal R2 values (fixed effects only), 
## while Rc2 labels for nlme models will be produced for every level of random
## effects which may be useful for looking at impact of random effects 
## on overall model fit for each level of random effects. 
r2.for.mm<-function(mm.object,round.r2.to=3,rel.dist.fr.top=.025,rel.dist.fr.left=.01,rel.spacing=.1,orientation="vertical"){
  ## mm.object        a mixed model object derived using nlme or lme4 packages
  ## round.r2.        sets how many decimal places you wish to show for your 
  ##                  R2 values.
  ## rel.dist.fr.top  sets the relative distance (relative to the range of y data)
  ##                  from the top that top of the first R2 label will line up with. 
  ##                  Consecutive R2 labels will appear below the first one
  ##                  based on the value of rel.spacing. 
  ## rel.dist.fr.left sets the relative distance (relative to the range in x data)
  ##                  from the left that the left side of the R2 labels will line
  ##                  up with. 
  ## rel.spacing      relative distance between R2 label values
  # set proportion from top here
  p.t<-1-(rel.dist.fr.top)
  p.l<-1-(rel.dist.fr.left)
  # separate function to round R2 values to desired decimal place:
  specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
  ## take fixed effect fitted values (level=0)
  # of random effects levels (for nested data):
  r.levels<-length(ranef(mm.object))
  class<-class(mm.object)[1]
  if(class=="lmerMod"){
    fixs<-(model.matrix(mm.object, type = "fixed") %*% fixef(lme4.mod))
  }else{
    fixs<-fitted(mm.object,level=0)
  }
  
  fixs<-fixs[!is.na(fixs)]
  Fix <- as.numeric(var(fixs))
  resids<-resid(mm.object)
  resi<-resids[!is.na(resids)]
  ## store residuals as res
  res<-var(resi)
  rcs.val<-c()
  rcs.lab<-c()
  gg.rcs.lab<-c()

  if(class=="lmerMod"){
    ## only return R2 for fixed effects model and full model
    ## too much work to get R2 for each level of random effects
  fit.ran<-fitted(mm.object)
  fit.ran1<-fit.ran[!is.na(fit.ran)]
  ## subtract fixed only from full model to estimate impact of random effects
  rans<- abs(var(fit.ran1)-(Fix)) 
  ## calculate the variance explained by fixed effects only - marginal R2
  rcs<-specify_decimal((Fix+rans)/(Fix+rans+res),round.r2.to)
  rcs.val<-c(rcs.val,paste("Rc2 = ",eval(parse(text="rcs"))))
  # These labels work for base plots:
  rcs.lab<-c(rcs.lab,bquote(italic(R)[c]^"2"~' = '~.(rcs)))
  # These labels work for ggplot2 plots:
  gg.rcs.lab<-c(gg.rcs.lab,paste("italic(R)[c]^",2,"~' = '~",rcs,sep=""))
}else{
  for(i in 1:r.levels){

      fit.ran<-fitted(mm.object,level=i)
      fit.ran1<-fit.ran[!is.na(fit.ran)]
      ## subtract fixed only from full model to estimate impact of random effects
      rans<- abs(var(fit.ran1)-(Fix)) 
      ## calculate the variance explained by fixed effects only - marginal R2
      rcs<-specify_decimal((Fix+rans)/(Fix+rans+res),round.r2.to)
      rcs.val<-c(rcs.val,paste("Rc2 = ",eval(parse(text="rcs"))))
      # These labels work for base plots:
      i.val<-as.numeric(i)
      rcs.lab<-c(rcs.lab,bquote(italic(R)[c*.(i.val)]^"2"~' = '~.(rcs)))
      # These labels work for ggplot2 plots:
      gg.rcs.lab<-c(gg.rcs.lab,paste("italic(R)[c",i,"]^",2,"~' = '~",rcs,sep=""))
    
  }
}
  Rm2<-specify_decimal(Fix/(Fix+rans+res),round.r2.to)
  # this orientation assumes that you want to plot how well the 
  # modeled values predict the dependent variable.
  y.max<-max(fitted(mm.object)+residuals(mm.object))
  x.max<-max(fitted(mm.object))
  y.min<-min(fitted(mm.object)+residuals(mm.object))
  x.min<-min(fitted(mm.object))
  y.range<-y.max-y.min
  x.range<-x.max-x.min
  ## now calculate the variance explained by the full model with random effects
  R2s<-list()
  R2s[[1]]<-c((paste("Rm2", "=",eval(parse(text = "Rm2")))),rcs.val)
  R2s[[2]]<-c(bquote(italic(R)[m]^"2"~'='~.(Rm2)),rcs.lab)
  if(class(mm.object)[1]=="lmerMod"){
    if(orientation=="vertical"){
      R2s[[3]]<-data.frame(x=c(rep(x.max-x.range*p.l,2)),y=c(y.max-y.range*(1-(seq(p.t*1000,(p.t-rel.spacing)*1000,length.out=2)/1000))),r2.labels=c(paste("italic(R)[m]^",2,"~' = '~",Rm2,sep=""),gg.rcs.lab))
      }
    else{
      R2s[[3]]<-data.frame(x=c(x.max-x.range*((seq(p.l*1000,(p.l-rel.spacing)*1000,length.out=2)/1000))),y=c(rep(y.min+y.range*p.t,2)),r2.labels=c(paste("italic(R)[m]^",2,"~' = '~",Rm2,sep=""),gg.rcs.lab))
    }
   
  }
  else{
    if(orientation=="vertical"){
      R2s[[3]]<-data.frame(x=c(rep(x.max-x.range*p.l,
                            (r.levels+1))),
                            y=c(y.min-y.range*(1-seq(p.l*1000,(p.l*1000-(r.levels*rel.spacing)*1000),
                            length.out=(r.levels+1)))/1000),
                            r2.labels=c(paste("italic(R)[m]^",2,"~' = '~",Rm2,sep=""),gg.rcs.lab))
    }
    else{
      R2s[[3]]<-data.frame(x=c(x.max-x.range*(seq(p.l*1000,(p.l*1000-(r.levels*rel.spacing)*1000),
                               length.out=(r.levels+1)))/1000),
                               y=c(rep(y.min+y.range*p.t,(r.levels+1))),
                              r2.labels=c(paste("italic(R)[m]^",2,"~' = '~",Rm2,sep=""),gg.rcs.lab))
    }
   
  }
 
  names(R2s)<-c("R2.values","R2.base.plot.labels","R2.ggplot.labels.df")

  return(R2s)
}


## Fit models to the made up data: ####
nlme.mod<-nlme(dep.var~a+b*indep.var,
               random=a~1|SITE/plot,
               fixed=list(a~1,b~1),
               start=c(a=10,b=5),
               data=test.data)
lme.mod<-lme(dep.var~indep.var,
               random=~1|SITE/plot,
               data=test.data)

lme4.mod<-lmer(dep.var~indep.var+(1|SITE/plot),
               data=test.data)


## Pick a model that you want to create R2 values and labels for: ####
model1<-nlme.mod

## Create labels:
r2.labels<-r2.for.mm(model1)

## Add predictions to original dataframe (for graphing predicted versus observed)
test.data$pred<-predict(model1)

# Plots in base plot:
plot(test.data$pred,test.data$dep.var)
text(0.2*max(test.data$pred),.95*max(test.data$dep.var),r2.labels$R2.base.plot.labels[[1]],pos=4)
text(0.2*max(test.data$pred),.85*max(test.data$dep.var),r2.labels$R2.base.plot.labels[[2]],pos=4)
## Additional R2 labels only exist for nlme models with more 
## than one level of random effects so this line is separate from the other two:
text(0.2*max(test.data$pred),.75*max(test.data$dep.var),r2.labels$R2.base.plot.labels[[3]],pos=4)


## For plots in ggplot2:
r2.labels<-r2.for.mm(model1,
                     round.r2.to = 3,
                     rel.dist.fr.top = 0,
                     rel.dist.fr.left = 0,
                     rel.spacing = .13,
                     orientation="horizontal")
lab.data<-as.data.frame(r2.labels$R2.ggplot.labels.df)

ggplot(test.data,aes(pred,dep.var))+geom_point(aes(color=SITE))+
  geom_smooth(color="blue",method="loess")+
  labs(x="Predicted values", y="Observed values")+
  geom_label(data=lab.data,aes(x,y,label=r2.labels),parse=T,hjust=0,color="blue")+
  theme_classic()

## Or try with vertical orientation:
## it is very easy to rerun the graph and just change
## the label function values. You will need to change spacing
## between vertical and horizontal orientation as I did not 
## set spacing from the end of one R label to the beginning
## of another - although that would be great!
r2.labels<-r2.for.mm(model1,
                     round.r2.to = 3,
                     rel.dist.fr.top = 0,
                     rel.dist.fr.left = 0,
                     rel.spacing = .08,
                     orientation="vertical")
lab.data<-as.data.frame(r2.labels$R2.ggplot.labels.df)

ggplot(test.data,aes(pred,dep.var))+geom_point(aes(color=SITE))+
  geom_smooth(color="blue",method="loess")+
  labs(x="Predicted values", y="Observed values")+
  geom_label(data=lab.data,aes(x,y,label=r2.labels),parse=T,hjust=0,color="blue")+
  theme_classic()

