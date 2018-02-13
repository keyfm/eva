library(ape)
library(caper)
library(car)
#source a couple of functions needed (I sent you all of them):
source("/home/roger/r_functions/build_all_models.r")
source("/home/roger/r_functions/get_conf_set.r")
source("/home/roger/r_functions/helpers.r")
source("/home/roger/roger/2013/pgls_assumptions/all_stuff/pgls_check_fncs.r")

#read data:
all.data=read.table(file="/home/roger/roger/2016/to_do/felix/felix_temp_1000g_annual_under26_wLAT_w1kgAC_nCEU.tsv", header=T, sep="\t")
#aggregate response as average proportion of allelle der per population:
test.data=aggregate(all.data$der/(all.data$der+all.data$anc), all.data[, c("annual", "under26", "latitude", "pop")], mean)
names(test.data)[ncol(test.data)]="allel.freq"#name column with the response
test.data$N=aggregate(!is.na(all.data$der/(all.data$der+all.data$anc)), all.data[, c("annual", "under26", "latitude", "pop")], sum)$x
test.data$N#are all 61, if not I'd consider weighting
rownames(test.data)=test.data$pop#add rownames to test.data (not sure whether this is needed)
#inspect distributions of predictor and response:
par(mfrow=c(2, 2))
hist(test.data$annual)
hist(test.data$under26)
hist(test.data$latitude)
hist(test.data$allel.freq)
#=> all +/-okay
par(mfrow=c(1, 1))
#check for absence of collinarity:
xx=vif(lm(allel.freq~annual+under26+latitude, data=test.data))
wt.txt(c.tab(as.data.frame(t(xx)), 3))#after you have called this you can paste the following into the editor (cool stuff, isn't it?)
# annual under26 latitude
# 7.592    4.485    3.049

#some VIF are pretty large which could perhaps lead to false negatives
#insepct scatter plots of the three predictors (to assess whether we could still hope for reasonable estimates):
plot(test.data$annual, test.data$under26)#not good
plot(test.data$annual, test.data$latitude)#somewhat okay
plot(test.data$under26, test.data$latitude)#not too bad

#read the phylogeny:
phylo=read.nexus("/home/roger/roger/2016/to_do/felix/new1000G_nexus_tree_fst_drops_rooted.nex")
plot(phylo)#just as a check
setdiff(phylo$tip.label, rownames(test.data))#check phylogeny and test.data for population name consistency
setdiff(rownames(test.data), phylo$tip.label)#check phylogeny and test.data for population name consistency
#create object needed as input for the phylogenetic regression:
test.data=comparative.data(phy = phylo, data = test.data, names.col = pop, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
full=pgls(allel.freq~annual+under26+latitude, data=test.data, lambda="ML")#fit full model
#assumptions:
diagnostics.plot(mod.res=full)#normality and homogeneity of the residuals (looks +/- okay)
#model stability:
model.stab=influence.pgls(pgls.res=full)
model.stab$taxa.failed#one population excluded lead to a problem withh fitting the model
wt.txt(c.tab(model.stab$summary, 4))#after you have called this you can paste the following into the editor (cool stuff, isn't it?)
#                orig     min    max
# (Intercept) -0.0419 -0.0893 0.1731
# annual       0.0078  0.0036 0.0087
# under26      0.1560  0.0819 0.1804
# latitude     0.0068  0.0036 0.0084
#=> kind of okay (but certainly not perfect)

#full null model comparison:
null=pgls(allel.freq~1, data=test.data, lambda=summary(full)$param["lambda"])#fit full model
xx=as.data.frame(anova(null, full, test="Chisq"))
wt.txt(c.tab(xx, 3))#after you have called this you can paste the following into the editor (cool stuff, isn't it?)
#   Res.Df    RSS    Df Sum of Sq Pr(>Chi)
# 1 18.000 18.607    NA        NA       NA
# 2 15.000  7.988 3.000    10.620    0.000
#=>full model is clearly significantly better than null model

#extract coefficients etc. from full model:
wt.txt(c.tab(summary(full)$coefficients, 3))#after you have called this you can paste the following into the editor (cool stuff, isn't it?)
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)   -0.042      0.204  -0.205    0.840
# annual         0.008      0.004   1.855    0.083
# under26        0.156      0.098   1.585    0.134
# latitude       0.007      0.002   3.797    0.002
#=>latitude is clearly significant, annual reveals a trend
#however, we need to be aware of the collinearity issue; so lets see what a model lacking latitude reveals for the other two:
full2=pgls(allel.freq~annual+under26, data=test.data, lambda="ML")#fit full model
wt.txt(c.tab(summary(full2)$coefficients, 3))
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)    0.365      0.250   1.461    0.163
# annual        -0.001      0.004  -0.360    0.723
# under26        0.068      0.123   0.559    0.584
#=>it doesn't seem that latitude is taking effects away from the others (more like the other way round)
#and what happens if we remove the other two (one at a time)?
full2=pgls(allel.freq~under26+latitude, data=test.data, lambda="ML")#fit full model
wt.txt(c.tab(summary(full2)$coefficients, 3))
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)    0.324      0.163   1.986    0.064
# under26        0.006      0.065   0.094    0.926
# latitude       0.003      0.001   2.517    0.023
full2=full=pgls(allel.freq~annual+latitude, data=test.data, lambda="ML")#fit full model
wt.txt(c.tab(summary(full2)$coefficients, 3))
#             Estimate Std. Error t value Pr(>|t|)
# (Intercept)    0.279      0.171   1.626    0.124
# annual         0.002      0.002   0.704    0.491
# latitude       0.004      0.001   2.697    0.016
#collinearity doesn't seem much of an issue

#now the multi model inference:
#construct vector with all models to be fitted:
allmodels=all.models(model="annual+under26+latitude")
allmodels=paste("allel.freq", allmodels, sep="~")
#determine AIC of all models:
all.aic=unlist(lapply(allmodels, function(x){
	x=as.formula(x)
	ifull=pgls(x, data=test.data, lambda="ML")
	return(ifull["aicc"])
}))
c.set=conf.set(aic=all.aic)#determine Akaike weights etc.
rownames(c.set)=all.models(model="annual+under26+latitude")#add rownames to c.set (for ease of interpretation
c.set=c.set[order(c.set$aic), ]#order according to model rank
wt.txt(c.tab(c.set, 3))#after you have called this you can paste the following into the editor (cool stuff, isn't it?)
#                               aic d.aic w.aic   cum c.set m.rank
# 1+latitude                -44.679 0.000 0.541 0.541 1.000  1.000
# 1+annual+latitude         -42.409 2.270 0.174 0.715 1.000  2.000
# 1+under26+latitude        -41.839 2.839 0.131 0.846 1.000  3.000
# 1+annual+under26+latitude -41.082 3.597 0.090 0.935 1.000  4.000
# 1+under26                 -38.351 6.328 0.023 0.958 1.000  5.000
# 1+annual                  -38.137 6.541 0.021 0.979 0.000  6.000
# 1                         -37.541 7.138 0.015 0.994 0.000  7.000
# 1+annual+under26          -35.654 9.024 0.006 1.000 0.000  8.000

#determine summed Akaike weights per predictor:
allmodels=all.models(model="annual+under26+latitude")
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
#     1 annual under26 latitude
# 1.000  0.290   0.249    0.935
apply(!is.na(model.weights), 2, mean, na.rm=T)#expected Akaike weights per predictor

#and finally a plot:
par(mar=c(3, 3, 0.2, 0.2), mgp=c(1.7, 0.4, 0), tcl=-0.25, las=1)
plot(test.data$data$latitude, test.data$data$allel.freq, xlab="latitude", ylab="porportion alleles with mutation", ylim=c(0, 1))
#get the model in there
#fit model with all predictors except latitude being centered:
plot.model=pgls(allel.freq~as.vector(scale(annual))+as.vector(scale(under26))+latitude, data=test.data, lambda="ML")
plot.model=coefficients(plot.model)#extract coefficients
xvals=seq(from=min(test.data$data$latitude), to=max(test.data$data$latitude), length.out=100)#create vecgtor with latitude values from min to max latitude
yvals=plot.model["(Intercept)"]+plot.model["latitude"]*xvals#extract fitted values
lines(xvals, yvals, lty=2)#add model lines
#the plot quite nicely shows why Gaussian models are not such a good idea for proportions...
savePlot("/home/roger/roger/2016/to_do/felix/pgls_result.pgn", type="png")
save.image("/home/roger/roger/2016/to_do/felix/felix_model_pgls.RData")
