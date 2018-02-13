library(lme4)
library(car)
library(gtools)
#source a couple of functions needed:
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/diagnostic_fcns.r")
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/glmm_stability.r") 
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/helpers.r")
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/build_all_models.r")
source("/home/felix_schulze/projects/trpm8_proj/glmm/roger/pgls/get_conf_set.r")

#read data:
all.data=read.table(file="~/projects/trpm8_proj/temp/tempCRUgrid_1000g_annual_under15_wLAT_w1kgAC.tsv", header=T, sep="\t")
#data prep (create random effect for subject ID):
s.id=as.factor(1:nrow(all.data))
#read fille with genetic distances:
gen.dist=read.table(file="/home/felix_schulze/projects/trpm8_proj/glmm/averageFST_pairwise_glmmPOPs_wCEU.tsv", header=T, sep="\t")
## align gen.dist such that it is matching all.data (i.e., one row for each entry in all.data, indicating the resp. individuals distance to the YRI population):
gen.dist=gen.dist[as.character(all.data$pop), "YRI"]
dist.to.YRI=as.vector(gen.dist)#turn gen.dist into vector named dist.to.YRI (not really needed but easier to be integrated in the following code)
tapply(dist.to.YRI, all.data$pop, sd)#just a check (should be 0 throughout)
tapply(dist.to.YRI, all.data$pop, mean)
#fit the full model:
aida.full=glmer(cbind(der, anc)~annual+latitude+dist.to.YRI+(1|s.id)+(1|pop), family=binomial, data=all.data)
vif(lm(der~annual+latitude+dist.to.YRI, data=all.data))#check for collinearity issue
   ##   annual    latitude dist.to.YRI 
   ## 4.638225    6.177039    1.795595 
#=> looks ok
ranef.diagn.plot(aida.full)#check distribution of random effects: looks okay
overdisp.test(aida.full)#check for absence of overdisprsion: great (shouldn't be much larger than 1)
##     chisq   df         P dispersion.parameter
## 1 1181.507 1216 0.7557244            0.9716339
#model stability:
full.stab.wac=glmm.model.stab(model.res=aida.full, data=data.frame(all.data, s.id, dist.to.YRI), use="pop")
table(full.stab.wac$detailed$warnings)#check whether there are any problems with non-converging models (there are none in this case)
## [3 pops cause issues if removed.]
wt.txt(c.tab(full.stab.wac$summary[, -1], 3))
##                     orig    min    max
## (Intercept)      -5.038 -7.896 -4.596
## annual            0.061  0.045  0.140
## latitude          0.086  0.077  0.124
## dist.to.YRI       9.370  7.751 10.579
## s.id@(Intercept)  0.096  0.000  0.453
## pop@(Intercept)   0.535  0.672  0.744
pdf('/home/felix_schulze/projects/trpm8_proj/pdf/glmm/glmm_model_stability.pdf')
m.stab.plot(full.stab.wac$summary[, -1])
dev.off()
## fit null model (full model comparison):
aida.null=glmer(cbind(der, anc)~dist.to.YRI+(1|s.id)+(1|pop), family=binomial, data=all.data)
xx=as.data.frame(anova(aida.null, aida.full, test="Chisq"))#full null model comparison
c.tab(xx, 5)
## [1] "#                Df        AIC        BIC     logLik   deviance    Chisq  Chi Df Pr(>Chisq)"
## [2] "# aida.null 4.00000 1946.16087 1966.58729 -969.08043 1938.16087       NA      NA         NA"
## [3] "# aida.full 6.00000 1929.19651 1959.83615 -958.59826 1917.19651 20.96436 2.00000    0.00003"
#=> is clearly significant 
wt.txt(c.tab(summary(aida.full)$coefficients, 3))#look at model coefficients
## [1] "#             Estimate Std. Error z value Pr(>|z|)"
## [2] "# (Intercept)   -5.038      1.250  -4.032    0.000"
## [3] "# annual         0.061      0.040   1.512    0.131"
## [4] "# latitude       0.086      0.020   4.282    0.000"
## [5] "# dist.to.YRI    9.370      2.954   3.172    0.002"
as.data.frame(summary(aida.full)$varcor)#inspect contribution of random effects
#multi model inference:
all.model.eqs=all.models(model="annual+latitude")
all.model.eqs.aida=paste("cbind(der, anc)~", all.model.eqs, "+dist.to.YRI+(1|s.id)+(1|pop)", sep="")
all.aic.aida=unlist(lapply(all.model.eqs.aida, function(x){
  summary(glmer(as.formula(x), family=binomial, data=all.data))$AICtab["AIC"]
}))
c.set.aida=conf.set(aic=all.aic.aida)
rownames(c.set.aida)=gsub(x=gsub(x=all.model.eqs, pattern="cbind(der, anc)~", replacement="", fixed=T), pattern="+dist.to.YRI+(1|s.id)+(1|pop)", replacement="", fixed=T)
c.set.aida=c.set.aida[order(c.set.aida$aic), ]
c.tab(c.set.aida, 3)
## [1] "#                        aic  d.aic w.aic   cum c.set m.rank"
## [2] "# 1+annual+latitude 1929.197  0.000 0.510 0.510 1.000  1.000"
## [3] "# 1+latitude        1929.286  0.090 0.488 0.997 1.000  2.000"
## [4] "# 1+annual          1939.830 10.633 0.003 1.000 0.000  3.000"
## [5] "# 1                 1946.161 16.964 0.000 1.000 0.000  4.000"
## [Null model out of best model confid. set, latitude best predictor but annual adds a bit]
## write.table(c.tab(c.set.aida, 4), file="/home/felix_schulze/projects/trpm8_proj/glmm/res_glmm_LatAnnualDistYRI_modern.tsv",sep="\t",quote=F,row.names=F,col.names=F)


