library(lme4)
library(car)
library(gtools)
#source a couple of functions needed:
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/diagnostic_fcns.r")
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/glmm_stability.r") 
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/helpers.r")
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/build_all_models.r")
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/get_conf_set.r")
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/pgls_check_fncs.r")

#read data:
all.data=read.table(file="~/projects/trpm8_proj/temp/tempCRUgrid_1000g_annual_under15_wLAT_w1kgAC.tsv", header=T, sep="\t")
## all.data <- all.data[ all.data[,1] != "CEU" , ]
#aggregate response as average proportion of allelle der per population:
test.data=aggregate(all.data$der/(all.data$der+all.data$anc), all.data[, c("annual", "under15", "latitude", "pop")], mean)
names(test.data)[ncol(test.data)]="allel.freq"#name column with the response
test.data$N=aggregate(!is.na(all.data$der/(all.data$der+all.data$anc)), all.data[, c("annual", "under15", "latitude", "pop")], sum)$x
test.data$N#are all 61, if not I'd consider weighting
rownames(test.data)=test.data$pop#add rownames to test.data (not sure whether this is needed)
#inspect distributions of predictor and response:
par(mfrow=c(2, 2))
hist(test.data$annual)
hist(test.data$under15)
hist(test.data$latitude)
hist(test.data$allel.freq)
#=> all +/-okay [extreme outliers would worry (under15)]
par(mfrow=c(1, 1))
#check for absence of collinarity [under15 too much coll...kick out!]:
xx=vif(lm(allel.freq~annual+latitude, data=test.data))
wt.txt(c.tab(as.data.frame(t(xx)), 3))#after you have called this you can paste the following into the editor (cool stuff, isn't it?)
# annual latitude
# 4.224     4.224
## [ok]

#insepct scatter plots of the three predictors (to assess whether we could still hope for reasonable estimates):
plot(test.data$annual, test.data$under15)#
plot(test.data$annual, test.data$latitude)#somewhat okay
plot(test.data$under15, test.data$latitude)#not too bad

#read the phylogeny:
phylo=read.nexus("/home/felix_schulze/projects/trpm8_proj/glmm/1000G_nexus_tree_fst_dropsButCEU_rooted.nex")
plot(phylo)#just as a check
setdiff(phylo$tip.label, rownames(test.data))#check phylogeny and test.data for population name consistency
setdiff(rownames(test.data), phylo$tip.label)#check phylogeny and test.data for population name consistency
#create object needed as input for the phylogenetic regression:
test.data=comparative.data(phy = phylo, data = test.data, names.col = pop, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
test.data$data$z.annual=as.vector(scale(test.data$data$annual))
test.data$data$z.latitude=as.vector(scale(test.data$data$latitude))
test.data$data$s.lat=(test.data$data$latitude-min(test.data$data$latitude))/diff(range(test.data$data$latitude))
test.data$data$s.ann=(test.data$data$annual-min(test.data$data$annual))/diff(range(test.data$data$annual))

require(caper)
require(geiger)
full=pgls(allel.freq~z.annual+z.latitude, data=test.data, lambda="ML")#fit full model
full=pgls(allel.freq~s.ann+s.lat, data=test.data, lambda="ML")#fit full model
## assumptions:
pdf('/home/felix_schulze/projects/trpm8_proj/pdf/pgls/fitted_hist_qqplot_vsFittedValues.pdf')
diagnostics.plot(mod.res=full)#normality and homogeneity of the residuals (looks +/- okay)
savePlot(file="/home/roger/roger/2016/to_do/felix_dec_pgls_diagnostics.png")
dev.off()
#model stability:
model.stab=influence.pgls(pgls.res=full)
model.stab$taxa.failed#one population excluded lead to a problem withh fitting the model
wt.txt(c.tab(model.stab$summary, 4))
#               orig    min    max
# (Intercept) 0.1645 0.1114 0.1864
# s.ann       0.1378 0.1113 0.1796
# s.lat       0.3647 0.3390 0.4295
m.stab.plot(model.stab$summary)
#=>looks pretty good to me

#full null model comparison:
null=pgls(allel.freq~1, data=test.data, lambda=summary(full)$param["lambda"])#fit full model
xx=as.data.frame(anova(null, full, test="Chisq"))
wt.txt(c.tab(xx, 3))
#   Res.Df    RSS    Df Sum of Sq Pr(>Chi)
# 1 19.000 28.518    NA        NA       NA
# 2 17.000 15.474 2.000    13.044    0.001
## =>full model is clearly significantly better than null model

#extract coefficients etc. from full model:
wt.txt(c.tab(summary(full)$coefficients, 3))
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)    0.164      0.184   0.893    0.384
# s.ann          0.138      0.089   1.541    0.142
# s.lat          0.365      0.127   2.863    0.011

##=>latitude is clearly significant, annual is not
#however, we need to be aware of the collinearity issue; so lets see what a model lacking latitude reveals for the other two:
cor(test.data$data$annual, test.data$data$latitude)
plot(test.data$data$annual, test.data$data$latitude)
##they are quite corelated, so lets have a look at what is revealed for annual when latitude isn't in the model:
red=pgls(allel.freq~s.ann, data=test.data, lambda="ML")#fit full model
wt.txt(c.tab(summary(red)$coefficients, 3))
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)    0.464      0.180   2.583    0.019
# s.ann         -0.093      0.044  -2.098    0.050
##annual is significant now, and its estimate changed is sign
##still, the full model clearly supports latitude

#now the multi model inference:
#construct vector with all models to be fitted:
allmodels=all.models(model="s.ann+s.lat")
allmodels=paste("allel.freq", allmodels, sep="~")
#determine AIC of all models:
all.aic=unlist(lapply(allmodels, function(x){
  x=as.formula(x)
  ifull=try(pgls(x, data=test.data, lambda="ML"), silent=T)
  if(class(ifull)[1]!="try-error"){
    return(ifull["aicc"])
  }else{
    return(NA)
  }
}))
c.set=conf.set(aic=all.aic)#determine Akaike weights etc.
rownames(c.set)=all.models(model="annual+latitude")#add rownames to c.set (for ease of interpretation
c.set=c.set[order(c.set$aic), ]#order according to model rank
wt.txt(c.tab(c.set, 3))
##                       aic d.aic w.aic   cum c.set m.rank
## 1+latitude        -49.430 0.000 0.504 0.504 1.000  1.000
## 1+annual+latitude -49.186 0.244 0.446 0.950 1.000  2.000
## 1+annual          -44.147 5.283 0.036 0.986 0.000  3.000
## 1                 -42.255 7.175 0.014 1.000 0.000  4.000
##                       aic d.aic w.aic   cum c.set m.rank
## laitutde best PV
## determine summed Akaike weights per predictor:
allmodels=all.models(model="annual+latitude")
xterms=unlist(strsplit(allmodels[[length(allmodels)]], split="+", fixed=T))
xx=lapply(xterms, function(x){
  xx=rep(NA, length(allmodels))
  xx[grepl(x=rownames(c.set), pattern=x, fixed=T)]=
    c.set$w.aic[grepl(x=rownames(c.set), pattern=x, fixed=T)]
  return(xx)
})
model.weights=matrix(unlist(xx), nrow=length(allmodels), byrow=F)
colnames(model.weights)=xterms
xx=apply(model.weights, 2, sum, na.rm=T)
wt.txt(c.tab(t(xx), 3))#extract Akaike weights per predictor
## 1     annual latitude
## 1.000  0.482    0.950
##pretty clear: strong support for latitude


