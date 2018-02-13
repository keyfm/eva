drop1p<-function(model.res, para=F, data=NULL, contr=NULL, n.cores=c("all-1", "all"), to.del=NULL, return.model.results=F){
	#last changed Mai 3rd 2016
	#print("this function was recently (Mai 3rd 2016) modified - could maybe someone compare its results with those of the older version and tell me what it revealed?")
	##determine terms that can be dropped:
	xcall=as.character(model.res@call)[2]
	#xcall=gsub(x=xcall, pattern="\n", replacement=T)
	model.terms=attr(terms(as.formula(xcall)), "term.labels")
	model.terms=model.terms[!grepl(x=model.terms, pattern="|", fixed=T)]
	model.terms=strsplit(model.terms, split=":", fixed=T)
	if(length(to.del)==0){
		to.del=unlist(lapply(1:length(model.terms), function(m1){
			xx=unlist(lapply((1:length(model.terms))[-m1], function(m2){
				length(intersect(model.terms[[m1]], model.terms[[m2]]))==length(model.terms[[m1]])
			}))
			sum(xx)>0
		}))
		model.terms=model.terms[!to.del]
		model.terms=unlist(lapply(model.terms, paste, collapse=":"))
	}else{
		model.terms=to.del
	}
	if(class(model.res)[[1]]=="lmerMod"){
		xfam=NULL##determine family
		if(length(contr)==0){contr=lmerControl()}##create control object if there isn't any
		if(as.character(model.res@call)[3]=="T"){stop("model wasn't fitted with REML=F; you should do that first")}
	}else{
		xfam=model.res@resp$family$family##determine family
		if(grepl(xfam, pattern="Negative Binomial")){xfam="neg.bin"}
		if(length(contr)==0){contr=glmerControl()}##create control object if there isn't any
	}
	##prepare data:
	#model=paste(model.terms, collapse="+")##get model wrt fixed effects
	if(length(data)==0){##prepare data frame when it is not handed over to data
		ii.data=model.res@frame
		if(any(grepl(x=colnames(ii.data), pattern="offset(", fixed=T))){
			colnames(ii.data)=gsub(x=colnames(ii.data), pattern=")", replacement="", fixed=T)
			colnames(ii.data)=gsub(x=colnames(ii.data), pattern="offset(", replacement="", fixed=T)
		}
	}else{##and if data were handed over to data... (need to check how the model looks like when there were weights handed over (also for the model stability function))
		ii.data=data
	}
	##prepare the data, weights (only supported for LMMs):
	if(length(xfam)==0){
		ii.data$xweights=model.res@resp$weights##weights
	}else{
		ii.data$xweights=1
	}
	#all.models=c(xcall, paste(xcall, model.terms, sep="-"))
	all.models=paste(xcall, model.terms, sep="-")
	##define function doing the model comparison:
	model.fun<-function(i.model, idata, contr, model.res, xfamily, xweights){
		i.model=as.formula(i.model)
		if(length(xfam)==0){
			i.res=try(lmer(i.model, data=idata, weights=xweights, control=contr, REML=F), silent=T)
		}else if(xfam!="neg.bin"){
			i.res=try(glmer(i.model, family=xfamily, data=idata, weights=xweights, control=contr), silent=T)
		}else{
			i.res=try(glmer.nb(i.model, family=xfamily, data=idata, weights=xweights), silent=T)
		}
		return(i.res)
	}
	##run all models:
  if(para){
    require(parallel)
    n.cores=n.cores[1]
    cl <- makeCluster(getOption("cl.cores", detectCores()))
		if(n.cores!="all"){
			if(n.cores!="all-1"){n.cores=length(cl)-1}else{n.cores=length(cl)}
		}
		n.cores=min(c(n.cores, length(all.models)))
		if(n.cores<length(cl)){
			cl=cl[1:n.cores]
		}
		if(length(cl)<length(model.terms)){
			cl=cl[1:length(model.terms)]
		}
    parLapply(cl=cl, 1:length(cl), fun=function(x){
      return(invisible(library(lme4)))
    })
    x.all.res=parLapply(cl=cl, X=all.models, fun=model.fun, idata=ii.data, contr=contr, model.res=model.res, xfamily=xfam, xweights="xweights")
    parLapply(cl=cl, X=1:length(cl), fun=function(x){invisible(rm(list=ls()))})
    stopCluster(cl)
  }else{
		x.all.res=lapply(X=all.models, FUN=model.fun, idata=ii.data, contr=contr, model.res=model.res, xfamily=xfam, xweights="weights")
  }
  xclass=unlist(lapply(x.all.res, function(x){class(x)[[1]]}))
  failed.models=model.terms[xclass=="try-error"]
  x.all.res=x.all.res[xclass!="try-error"]
  all.tests=lapply(1:length(x.all.res), function(x){
		#xx=as.data.frame(anova(x.all.res[[x]]$value, x.all.res[[1]]$value))
		ires=c(Chisq=-2*(as.vector(logLik(x.all.res[[x]]))-as.vector(logLik(model.res))), "Chi Df"=length(fixef(model.res))-length(fixef(x.all.res[[x]])))#xx[2, c("AIC", "Chisq", "Chi Df", "Pr(>Chisq)")]
		ires=c(ires, "Pr(>Chisq)"=as.vector(pchisq(q=ires["Chisq"], df=ires["Chi Df"], lower.tail=F)))
		#ires$AIC=xx[1, "AIC"]
		return(list(logLik=as.vector(logLik(x.all.res[[x]])), AIC=extractAIC(x.all.res[[x]])[2], Chisq=as.vector(ires["Chisq"]), "Chi Df"=as.vector(ires["Chi Df"]), 
			"Pr(>Chisq)"=ires["Pr(>Chisq)"], n=length(residuals(x.all.res[[x]])),
			n.opt.warnings=length(unlist(x.all.res[[x]]@optinfo$conv$lme4)), n.fun.warnings=length(unlist(x.all.res[[x]]@optinfo$warnings))))
	})
  aic=ifelse(class(model.res)[[1]]=="lmerMod", extractAIC(model.res)[2], summary(model.res)$AICtab["AIC"])
	all.tests=data.frame(logLik=c(as.vector(logLik(model.res)), unlist(lapply(all.tests, function(x){x$logLik}))),
		AIC=c(aic, unlist(lapply(all.tests, function(x){x$AIC}))),
		Chisq=c(NA, unlist(lapply(all.tests, function(x){x$Chisq}))),
		"Chi Df"=c(NA, unlist(lapply(all.tests, function(x){x["Chi Df"]}))),
		"Pr(>Chisq)"=c(NA, unlist(lapply(all.tests, function(x){x["Pr(>Chisq)"]}))),
		n=c(length(residuals(x.all.res[[1]])), unlist(lapply(all.tests, function(x){x["n"]}))),
		n.opt.warnings=c(NA, unlist(lapply(all.tests, function(x){x$n.opt.warnings}))), n.fun.warnings=c(NA, unlist(lapply(all.tests, function(x){x$n.fun.warnings}))))
  rownames(all.tests)=c("none", model.terms[xclass!="try-error"])
  if(length(failed.models)>0){
		xx=matrix(NA, ncol=ncol(all.tests), nrow=length(failed.models))
		colnames(xx)=names(all.tests)
		rownames(xx)=failed.models
		all.tests=rbind(all.tests, xx)
	}
	#all.tests[ncol(all.tests)]=as.character(all.tests[, ncol(all.tests)])
	if(return.model.results){
		to.return=list(drop1.res=all.tests, model.results=x.all.res)
	}else{
		to.return=list(drop1.res=all.tests, model.results=NULL)
	}
  return(to.return)
}

