library(lme4)
library(car)
library(gtools)
#source a couple of functions needed:
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/diagnostic_fcns.r")
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/glmm_stability.r") 
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/helpers.r")
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/build_all_models.r")
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/get_conf_set.r")

##read data:
all.data=read.table(file="~/projects/trpm8_proj/sgdp/glmm/rs10166942_sgdp_Cteam_111pop2ind_meanT_diffT.txt", header=T, sep="\t")
#data prep (create random effect for subject ID):
s.id=as.factor(1:nrow(all.data))
#read fille with genetic distances:
gen.dist=read.table(file="/mnt/scratch/josh/simonds_data/fst/wc_sample2_vs_yoruba/chr_combined/summary/Yoruba_vs_other_wc_sample2_averaged_fsts.cee.ceemasked.heng.txt", header=F, sep="\t")
rownames(gen.dist) <- as.character(gen.dist[,2])
dist.to.YRI=as.vector(gen.dist[as.character(all.data$pop), "V7"])#turn gen.dist into vector named dist.to.YRI 
all.data$pop[is.na(dist.to.YRI)]
dist.to.YRI[is.na(dist.to.YRI)]=0
## fit the full model:
z.meanT=as.vector(scale(all.data$meanT))
z.diffT=as.vector(scale(all.data$diffT))
z.abs.latitude=as.vector(scale(abs(all.data$latitude)))
z.dist.to.YRI=as.vector(scale(dist.to.YRI))
aida.full=glmer(cbind(der, anc)~z.meanT+z.abs.latitude+z.dist.to.YRI+(1|s.id)+(1|pop), family=binomial, data=all.data)
vif(lm(der~z.meanT+z.abs.latitude+z.dist.to.YRI, data=all.data))#check for collinearity issue
      ##  z.meanT z.abs.latitude  z.dist.to.YRI 
      ## 3.732761       3.717651       1.093607 
#=> looks ok.
ranef.diagn.plot(aida.full)#check distribution of random effects: looks okay
overdisp.test(aida.full)#check for absence of overdisprsion: great (shouldn't be much larger than 1)
##      chisq  df         P dispersion.parameter
## 1 143.6761 216 0.9999583            0.6651671
## -> disp. param. looks good
#model stability:
full.stab.wac=glmm.model.stab(model.res=aida.full, data=data.frame(all.data,z.meanT,z.abs.latitude,z.dist.to.YRI, s.id), use="pop")
table(full.stab.wac$detailed$warnings)#check whether there are any problems with non-converging models (there are none in this case)
wt.txt(c.tab(full.stab.wac$summary[, -1], 3))#see PGLS script for hints about the functions wt.txt and c.tab. 
#                   orig   min   max
# (Intercept)      0.172 0.124 0.203
# z.meanT          0.859 0.732 1.019
# z.abs.latitude   1.379 1.312 1.574
# z.dist.to.YRI    0.661 0.598 0.752
# s.id@(Intercept) 0.002 0.000 0.070
# pop@(Intercept)  1.109 1.022 1.064
## [range of values obained by full.stab.wac; min max should not deviate too much from orig (same order of magnitude)]. 
pdf('/home/felix_schulze/projects/trpm8_proj/pdf/glmm/glmm_model_stability_sgdp.pdf')
m.stab.plot(full.stab.wac$summary[, -1])#looks pretty good
dev.off()
## fit null model (full model comparison):
aida.null=glmer(cbind(der, anc)~dist.to.YRI+(1|s.id)+(1|pop), family=binomial, data=all.data)
xx=as.data.frame(anova(aida.null, aida.full, test="Chisq"))#full null model comparison
c.tab(xx, 5)
## [1] "#                Df       AIC       BIC     logLik  deviance    Chisq  Chi Df Pr(>Chisq)"
## [2] "# aida.null 4.00000 452.40261 465.97712 -222.20131 444.40261       NA      NA         NA"
## [3] "# aida.full 6.00000 435.01965 455.38142 -211.50983 423.01965 21.38296 2.00000    0.00002"
#=> is clearly significant [compares null (only YRI distance) and full (all predictor)
wt.txt(c.tab(summary(aida.full)$coefficients, 3))#look at model coefficients
#                Estimate Std. Error z value Pr(>|z|)
# (Intercept)       0.172      0.161   1.069    0.285
# z.meanT           0.859      0.312   2.757    0.006
# z.abs.latitude    1.379      0.331   4.170    0.000
# z.dist.to.YRI     0.661      0.189   3.501    0.000
## => [lat pos, dis to YRI too,annual no influence]
## source("/home/roger/r_functions/drop1_para.r")
tests=drop1p(aida.full,data=data.frame(all.data,z.meanT,z.abs.latitude,z.dist.to.YRI, s.id))
wt.txt(c.tab(tests$drop1.res, 3))
#                  logLik     AIC  Chisq Chi.Df Pr..Chisq.       n n.opt.warnings n.fun.warnings
# none           -211.510 435.020     NA     NA         NA 220.000             NA             NA
# z.meanT        -215.365 440.730  7.710  1.000      0.005 220.000          0.000          0.000
# z.abs.latitude -220.794 451.588 18.569  1.000      0.000 220.000          2.000          0.000
# z.dist.to.YRI  -218.407 446.814 13.795  1.000      0.000 220.000          0.000          0.000
as.data.frame(summary(aida.full)$varcor)#inspect contribution of random effects

contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000))
aida.red.only.z.meanT=glmer(cbind(der, anc)~z.meanT+z.abs.latitude+z.dist.to.YRI+(1|s.id)+(1|pop), family=binomial, data=all.data, control=contr)
aida.red.only.z.diffT=glmer(cbind(der, anc)~z.diffT+z.abs.latitude+z.dist.to.YRI+(1|s.id)+(1|pop), family=binomial, data=all.data, control=contr)
wt.txt(c.tab(summary(aida.red.only.z.meanT)$coefficients, 3))
#                Estimate Std. Error z value Pr(>|z|)
# (Intercept)       0.172      0.161   1.069    0.285
# z.meanT           0.859      0.312   2.757    0.006
# z.abs.latitude    1.379      0.331   4.170    0.000
# z.dist.to.YRI     0.661      0.189   3.501    0.000

## [extra variance through inds]

## now the multi model inference without diffT:
all.model.eqs=all.models(model="z.meanT+z.abs.latitude")
all.model.eqs.aida=paste("cbind(der, anc)~", all.model.eqs, "+z.dist.to.YRI+(1|s.id)+(1|pop)", sep="")
all.aic.aida=unlist(lapply(all.model.eqs.aida, function(x){
  xx=glmer(as.formula(x), family=binomial, data=all.data, control=contr)
  summary(xx)$AICtab["AIC"]+aic.c.fac(N=nrow(all.data), k=length(fixef(xx)+nrow(as.data.frame(summary(xx)$varcor))))#with the correction for small samples
}))
c.set.aida=conf.set(aic=all.aic.aida)
rownames(c.set.aida)=gsub(x=gsub(x=all.model.eqs, pattern="cbind(der, anc)~", replacement="", fixed=T), pattern="+dist.to.YRI+(1|s.id)+(1|pop)", replacement="", fixed=T)
c.set.aida=c.set.aida[order(c.set.aida$aic), ]
wt.txt(c.tab(c.set.aida, 3))
#                              aic  d.aic w.aic   cum c.set m.rank
# 1+z.meanT+z.abs.latitude 435.206  0.000 0.943 0.943 1.000  1.000
# 1+z.abs.latitude         440.841  5.635 0.056 1.000 1.000  2.000
# 1+z.meanT                451.699 16.494 0.000 1.000 0.000  3.000
# 1                        452.458 17.252 0.000 1.000 0.000  4.000

## [Null model out of best model confid. set, latitude best predictor but annual adds a bit]
## write.table(c.tab(c.set.aida, 4), file="/home/felix_schulze/projects/trpm8_proj/glmm/res_glmm_LatAnnualDistYRI_modern_sgdp.tsv",sep="\t",quote=F,row.names=F,col.names=F)
##extract Akaike weights per predictor:
xterms=unlist(strsplit(all.model.eqs[[length(all.model.eqs.aida)]], split="+", fixed=T))
xx=lapply(xterms, function(x){
	xx=rep(NA, length(all.models))
	xx[grepl(x=rownames(c.set.aida), pattern=x, fixed=T)]=
		c.set.aida$w.aic[grepl(x=rownames(c.set.aida), pattern=x, fixed=T)]
	return(xx)
})
pv.weights.aida=matrix(unlist(xx), nrow=length(all.model.eqs), byrow=F)
colnames(pv.weights.aida)=xterms
wt.txt(c.tab(t(apply(pv.weights.aida, 2, sum, na.rm=T)), 3))
##      1   z.meanT z.abs.latitude
## [1,] 1 0.9434804      0.9995836


### 3D PLOT
source("/home/roger/r_functions/three_d_plots_new.r")
source("/home/roger/r_functions/int_plots_no_persp.r")
all.data$abs.latitude=abs(all.data$latitude)
contr=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=10000))
plot.full=glmer(cbind(der, anc)~z.meanT+z.abs.latitude+z.dist.to.YRI+(1|s.id)+(1|pop), family=binomial, data=all.data, control=contr)
draw.2.w.int.bw.cov.2(plot.data=data.frame(all.data, z.meanT, z.abs.latitude), vars=c("z.meanT", "z.abs.latitude"), coefs=fixef(plot.full), link="logit",
  grid.resol=8, var.names=c("mean temperature", "absolute latitude"), zlab="derived allele frequency", 
  zlim=c(0, 1), theta=30, phi=10, expand=0.6, r=10, cex.lab=1, response=all.data$der/(all.data$der+all.data$anc), size.fac=1, print.NA.cells=T, quiet=T, col="black")

range(all.data$meanT)
range(all.data$abs.latitude)
x.labels=seq(-5, 25, by=10)
x.at=(x.labels-mean(all.data$meanT))/sd(all.data$meanT)
y.labels=seq(10, 60, by=10)
y.at=(y.labels-mean(all.data$abs.latitude))/sd(all.data$abs.latitude)
z.at=seq(0, 1, by=0.2)
z.labels=as.character(z.at)
z.labels[nchar(z.labels)==1]=paste(z.labels[nchar(z.labels)==1], "0", sep=".")
intplot.2.cov(
	plot.data=data.frame(all.data, z.meanT, z.abs.latitude), vars=c("z.meanT", "z.abs.latitude"), coefs=fixef(plot.full), link="logit",
	response=all.data$der/(all.data$der+all.data$anc),
  grid.resol=8, theta=-30, phi=10, expand=0.6, d=50, 
	c("mean temperature", "absolute latitude"), zlab="derived allele frequency", show.axes="ticks", 
	x.at=x.at, x.labels=x.labels, y.at=y.at, y.labels=y.labels, z.at=z.at, z.labels=z.labels, 
	zlim=c(0, 1), cex.lab=1, cex.axis=1, axis.lab.shift.fac=1, arrow.dist=1, size.fac=1, col="black", background=NA, 
	print.NA.cells=T, quiet=T, bottom.grid=T, prep.for.anim=F, tcl=-0.15)
savePlot(file="/home/roger/roger/2017/to_do/felix_jan/three_d_plot.png", type="png")
dev.copy2pdf(file="/home/roger/roger/2017/to_do/felix_jan/three_d_plot.pdf")
save.image("/home/roger/roger/2017/to_do/felix_jan/glmm_newFst_full.RData")

