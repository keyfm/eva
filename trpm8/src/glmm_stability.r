require(lme4)
glmm.model.stab<-function(model.res, contr=NULL, ind.cases=F, para=F, data=NULL, use=NULL, n.cores=c("all-1", "all"), save.path=NULL){
	print("please carefully evaluate whether the result makes sense, and if not, please contact me")
	#function determining stability of GLMMs (run using lmer or glmer) by excluding levels of random effects, one at a time;
	#supports
		#weights, offset terms and random slopes;
	#does not support
		#correlations between random intercepts and random slopes
		#any terms calculated in the model formula (e.g., log, function I, etc.); but interactions do work
	#latest additions/modifications:
		#copes also with data sets where a level of a (fixed effects) factor is entirely dropped from the data
		#new way of catching warnings (are contained in the detailed output table)
		#includes sample size in the detailed output table
		#catches errors
		#use: an argument taking the names of the random effects for which model stability should be evaluated
			#(useful for models with a random effect having a unique case for each level of the data)
	#written by Roger Mundry
	#last modified April 2016 (added dealing with glmer.nb; fixed a bug in the output of the random effects of tthe original model)
	n.cores=n.cores[1]
  model.eq=as.formula(as.character(model.res@call)[2])
	if(class(model.res)[1]=="lmerMod"){
    REML=model.res@resp$REML==1
    xfam="gaussian"
  }else{
		xfam=model.res@resp$family$family
  }
  if(grepl(x=xfam, pattern="Negative Binomial", fixed=T)){
		xfam="neg.bin"
	}
  weights=model.res@resp$weights
  if(length(data)==0){
		ii.data=model.res@frame
		offs.col=grep(x=names(ii.data), pattern="offset(", fixed=T)
		if(length(offs.col)>0){
			for(i in offs.col){
				#ii.data[,i]=exp(ii.data[,i])
				names(ii.data)[i]=substr(x=names(ii.data)[i], start=8, stop=nchar(names(ii.data)[i]) -1)
			}
		}
		wght.col=grep(x=names(ii.data), pattern="(weights)", fixed=T)
		if(length(wght.col)>0){
			names(ii.data)[wght.col]=gsub(x=names(ii.data)[wght.col], pattern="(", replacement="", fixed=T)
			names(ii.data)[wght.col]=gsub(x=names(ii.data)[wght.col], pattern=")", replacement="", fixed=T)
		}else{
			ii.data$weights=weights
		}
	}else{
		ii.data=data.frame(data, weights)
	}
	if(substr(as.character(model.eq)[2], start=1, stop=6)=="cbind("){
		ii.data$weights=1
	}
  #if(xfam=="binomial"){
  #   response=as.character(model.eq)[2]
  #   if(sum(names(ii.data)==response)==0)
  #}
  ranefs=names(ranef(model.res))
	if(length(use)==0){use=ranefs}
	ranefs=ranefs[ranefs%in%use]
  xlevels=lapply(ranefs, function(x){return(as.vector(unique(ii.data[ ,x])))})
  ranefs=rep(ranefs, unlist(lapply(xlevels, length)))
  to.do=cbind(ranefs, unlist(xlevels))
  if(ind.cases){
    ii.data=data.frame(ii.data, ic=as.factor(1:nrow(ii.data)))
    to.do=rbind(to.do, cbind("ic", levels(ii.data$ic)))
  }
	keepWarnings <- function(expr) {
		localWarnings <- list()
		value <- withCallingHandlers(expr,
			warning = function(w) {
				localWarnings[[length(localWarnings)+1]] <<- w
				invokeRestart("muffleWarning")
			}
		)
		list(value=value, warnings=localWarnings)
	}
	if(length(contr)==0){
		if(class(model.res)[1]=="lmerMod"){
			contr=lmerControl()
		}else{	
			contr=glmerControl()
		}
	}
  ifun=function(x, model.res, to.do, ii.data, contr=contr){
    sel.ii.data=subset(ii.data, ii.data[,to.do[x, 1]]!=to.do[x, 2])
		if(class(model.res)[1]=="lmerMod"){
      sel.ii.res=try(keepWarnings(lmer(model.eq, data=sel.ii.data, weights=weights, REML=REML, control=contr)), silent=T)
    }else if(xfam=="binomial" | xfam=="poisson"){
      sel.ii.res=try(keepWarnings(glmer(model.eq, data=sel.ii.data, family=xfam, control=contr)), silent=T)
    }else{
			sel.ii.res=try(keepWarnings(glmer.nb(model.eq, data=sel.ii.data, control=contr)), silent=T)
		}
		if(length(save.path)>0){
			est.fixed.effects=fixef(sel.ii.res$value)
			est.random.effects=as.data.frame(summary(sel.ii.res$value)$varcor)
			model.warnings=sel.ii.res$warnings
			n=length(residuals(sel.ii.res$value))
			what=to.do[x, ]
			save(file=paste(c(paste(c(paste(c(save.path, "m"), collapse="/"), x), collapse="_"), ".RData"), collapse=""), list=c("what", "est.fixed.effects", "est.random.effects", "model.warnings", "n"))
		}
    if(class(sel.ii.res)!="try-error"){
      return(list(fere=c(fixef(sel.ii.res$value), sqrt(unlist(lapply(summary(sel.ii.res$value)$varcor, function(x){attr(x, "stddev")})))), N=length(residuals(sel.ii.res$value)), warnings=paste(unlist(sel.ii.res$warnings)$message, collapse="/")))
    }else{
      return(list(fere=rep(NA, length(fixef(model.res))+length(unlist(lapply(summary(model.res)$varcor, function(x){attr(x, "stddev")})))), N=NA, warnings=NA))
    }
  }
  if(para){
    require(parallel)
    cl <- makeCluster(getOption("cl.cores", detectCores()))
		if(n.cores!="all"){
			if(n.cores!="all-1"){n.cores=length(cl)-1}
			if(n.cores<length(cl)){
				cl=cl[1:n.cores]
			}
		}
    parLapply(cl=cl, 1:length(cl), fun=function(x){
      library(lme4)
      return(invisible(""))
    })
    all.coeffs=parLapply(cl=cl, X=1:nrow(to.do), fun=ifun, model.res=model.res, to.do=to.do, ii.data=ii.data, contr=contr)
    parLapply(cl=cl, X=1:length(cl), fun=function(x){rm(list=ls())})
    stopCluster(cl)
  }else{
    all.coeffs=lapply(1:nrow(to.do), ifun, model.res=model.res, to.do=to.do, ii.data=ii.data, contr=contr)
  }

	
	all.n=unlist(lapply(all.coeffs, function(x){x$N}))
	all.warnings=unlist(lapply(all.coeffs, function(x){paste(unlist(x$warnings), collapse=", ")}))
	all.coeffs=lapply(all.coeffs, function(x){
		x=x$fere
		to.ret=rep(NA, length(fixef(model.res)))
		names(to.ret)=names(fixef(model.res))
		to.ret[names(x)]=x
		return(to.ret)
	})
	xx=lapply(summary(model.res)$varcor, function(x){attr(x, "stddev")})
  xnames=c(names(fixef(model.res)), paste(rep(names(xx), unlist(lapply(xx, length))), unlist(lapply(xx, names)), sep="@"))
  # xx=lapply(all.coeffs, function(x){
    # xxr=rep(NA, length(xnames))
    # names(xxr)=xnames
    # xxr[match(names(x), xnames)]=x
    # return(xxr)
  # })
  # to.add=unlist(lapply(summary(model.res)$varcor, function(x){unlist(attr(x, "dimnames"))[1]}))
  # to.add=paste(names(summary(model.res)$varcor), to.add, sep="@")
	all.coeffs=matrix(unlist(all.coeffs), nrow=length(all.coeffs), byrow=T)
  colnames(all.coeffs)=xnames#c(names(fixef(model.res)), to.add)
	xsum=apply(all.coeffs, 2, range, na.rm=T)
  xsum=data.frame(what=colnames(all.coeffs), orig=c(fixef(model.res), 
		unlist(lapply(summary(model.res)$varcor, function(x){attr(x, "stddev")}))), t(xsum))
	rownames(xsum)=as.character(xsum$what)
  colnames(to.do)=c("ranef", "level")
  all.coeffs=data.frame(to.do, N=all.n, all.coeffs, warnings=all.warnings)
  names(xsum)[3:4]=c("min", "max")
  return(list(detailed=all.coeffs, summary=xsum))
}

