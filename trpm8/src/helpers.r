anova2<-function(red, full){#wrapper for anova, provides the same output for glm, glm.nb and lmer
  ll.test=anova(red, full, test="Chisq")
  #browser()
  if(length(grep(attr(ll.test, which="heading"), pattern="Negative Binomial Models"))>0){
    ll.test=as.data.frame(ll.test)[2, c(7, 6, 8)]
  }else if(attr(ll.test, which="heading")[1]=="Analysis of Deviance Table\n"){
    ll.test=as.data.frame(ll.test)[2, c(4, 3, 5)]
  }else if(class(full)[1]=="mer" & class(red)[1]=="mer"){
    ll.test=as.data.frame(ll.test)[2, 5:7]
  }else if(class(full)[1]=="glmerMod" & class(red)[1]=="glmerMod"){
    ll.test=as.data.frame(ll.test)[2, 6:8]
  }else{
    ll.test=as.data.frame(ll.test)[2, 6:8]
  }
  names(ll.test)=c("Chisq", "df", "P")
  return(data.frame(c(ll.test)))
}

c.lrt<-function(lrt, space=F, digits=c(2, 3)){#collapses the results of anova2 in a single entry
  #optional arguments control whether signs like = are surrounded by spaces
  #and the number of digits displayed
  if(space){
    collapse.e=" = "
    collapse.s=" < "
  }else{
    collapse.e="="
    collapse.s="<"
  }
  chi=round(lrt[1], digits=digits[1])
  P=round(lrt[3], digits=digits[2])
  paste(c(
    ifelse(chi>0, paste(c("c2", chi), collapse=collapse.e), paste(c("c2", paste(c("0.", rep(0, digits[1]-1), 1), collapse="")), collapse=collapse.s)),
    paste(c("df", lrt[2]), collapse=collapse.e),
    ifelse(P>0, paste(c("P", P), collapse=collapse.e), paste(c("P", paste(c("0.", rep(0, digits[2]-1), 1), collapse="")), collapse=collapse.s))
  ), collapse=", ")
}

c.est<-function(x, add.names=T, space=F, digits=c(2, 2, 2, 3)){
  #collapses the estimates in a single row of a coefficients object in a single entry
  #optional arguments control whether signs like = are surrounded by spaces
  #and the number of digits displayed
  if(space){
    collapse.e=" = "
    collapse.s=" < "
  }else{
    collapse.e="="
    collapse.s="<"
  }
  tsx.name=substr(names(x)[3], start=1, stop=1)
  est=round(x[1], digits=digits[1])
  se=round(x[2], digits=digits[2])
  tsx=round(x[3], digits=digits[3])
  P=round(x[4], digits=digits[4])
  xx.res=paste(c(
    paste(c(ifelse(abs(est)>0, est, paste(c("0.", rep(0, digits[1]-1), 1), collapse="")),
      ifelse(abs(se)>0, se, paste(c("0.", rep(0, digits[2]-1), 1), collapse=""))), collapse="+"),
      ifelse(abs(tsx)>0, paste(c(tsx.name, tsx), collapse=collapse.e), paste(c(tsx.name, paste(c("0.", rep(0, digits[3]-1), 1), collapse="")), collapse=collapse.s)),
      ifelse(P>0, paste(c("P", P), collapse=collapse.e), paste(c("P", paste(c("0.", rep(0, digits[4]-1), 1), collapse="")), collapse=collapse.s))
    ), collapse=", ")
  if(add.names){
    xx.res=paste(c("estimate+SE", xx.res), collapse="=")
  }
  return(xx.res)
}

c.odt<-function(x){#collapse overdispersion test
  paste(paste(c("dispersion parameter", "c2", "df", "P"), c(round(x[4], 2), round(x[1], 2), round(x[2], 0), round(x[3], 3)), sep="="), collapse=", ")
}

c.drop1<-function(lrt, space=F, digits=c(2, 3)){#collapses the results of drop1 in a single entry
  #optional arguments control whether signs like = are surrounded by spaces
  #and the number of digits displayed
  if(space){
    collapse.e=" = "
    collapse.s=" < "
  }else{
    collapse.e="="
    collapse.s="<"
  }
  chi=round(lrt["LRT"], digits=digits[1])
  P=round(lrt["Pr(Chi)"], digits=digits[2])
  paste(c(
    ifelse(chi>0, paste(c("c2", chi), collapse=collapse.e), paste(c("c2", paste(c("0.", rep(0, digits[1]-1), 1), collapse="")), collapse=collapse.s)),
    paste(c("df", lrt["Df"]), collapse=collapse.e),
    ifelse(P>0, paste(c("P", P), collapse=collapse.e), paste(c("P", paste(c("0.", rep(0, digits[2]-1), 1), collapse="")), collapse=collapse.s))
  ), collapse=", ")
}

c.F<-function(anova.res, space=F, digits=c(2, 3), apa=T){#collapses the results of anova2 in a single entry
  #optional arguments control whether signs like = are surrounded by spaces
  #and the number of digits displayed
  if(space){
    collapse.e=" = "
    collapse.s=" < "
  }else{
    collapse.e="="
    collapse.s="<"
  }
	anova.res=as.data.frame(anova.res)
  f=round(anova.res$F[2], digits=digits[1])
	df=paste(c(anova.res$Df[2], anova.res$Res.Df[2]), collapse=",")
	if(apa){df=paste(c("(", df, ")"), collapse="")}
  P=round(anova.res$"Pr(>F)"[2], digits=digits[2])
  paste(c(
    paste(c(paste(c("F", df), collapse=""), f), collapse=collapse.e), 
    ifelse(P>0, paste(c("P", P), collapse=collapse.e), paste(c("P", paste(c("0.", rep(0, digits[2]-1), 1), collapse="")), collapse=collapse.s))
  ), collapse=", ")
}

c.descr<-function(x, digits=3, sep=" and "){
	#determines and rounds (to digits as specified in digits) mean and sd of x, and then collapses them in a single entry
	xx=as.character(round(c(mean(x, na.rm=T), sd(x, na.rm=T)), digits=digits))
	if(sum(grepl(x=xx, pattern=".", fixed=T))>0){
		xx=matrix(unlist(strsplit(as.character(xx), split=".", fixed=T)), ncol=2, byrow=T)
		xx[, 2]=unlist(lapply(xx[, 2], function(x){
			paste(c(x, paste(c(rep("0", times=digits-nchar(x))), collapse="")), collapse="")
		}))
		xx=apply(xx, 1, paste, collapse=".")
	}
	paste(xx, collapse=sep)
}

c.tab<-function(x, digits=NA, n.spaces=1, add.hash=T, incl.fst=F){
	#last updated: 2015 Nov 25
	#collapses a matrix into a vector (with one entry for for each row of x)
	#rounds entries to digits when digits is not NA and adds spaces in between adjacent entries (including row and column names of x)
		#such that everything is aligned nicely
	#requires x to have row and column names
	if(!is.na(digits)){
		for(i in 1:ncol(x)){
			if(is.numeric(x[, i])){x[, i]=round(x[, i], digits=digits)}
		}
	}
	res=format(as.matrix(x))
	res=gsub(x=res, pattern=" ", replacement="", fixed=T)
	res=cbind(rownames(res), res)
	res=rbind(colnames(res), res)
	rownames(res)=NULL
	colnames(res)=NULL
	n.char.range=apply(res, 2, function(x){range(nchar(x))})
	if(incl.fst){
		res[, 1]=unlist(lapply(res[, 1], function(y){paste(c(rep(" ", n.char.range[2, 1]-nchar(y)), y), collapse="")}))
	}else{
		res[, 1]=unlist(lapply(res[, 1], function(x){paste(c(x, rep(" ", n.char.range[2, 1]-nchar(x))), collapse="")}))
	}
	res[, -1]=matrix(unlist(lapply(2:ncol(res), function(x){
		unlist(lapply(res[, x], function(y){paste(c(rep(" ", n.char.range[2, x]-nchar(y)), y), collapse="")}))
	})), nrow=nrow(res), byrow=F)
	res=apply(res, 1, paste, collapse=paste(rep(" ", n.spaces), collapse=""))
	if(add.hash){res=paste("#", res, sep=" ")}
	return(res)
}

get.ref<-function(package){
  #wrapper for citation() which formats the output in a nicer way
  #is intended to work with books (e.g., MASS), 'normal' packages (e.g., lme4) and R
  #but not really tested...
  #argument is character
  if(nchar(package)>0){xx=citation(package)}else{xx=citation()}
	if(nchar(package)>0){
		all.res=rep(NA, length(xx))
		for(i in 1:length(xx)){
			xres=unlist(lapply(xx[[i]]$author, function(x){
				yy=paste(unlist(lapply(strsplit(x$given, split=" ", fixed=T), function(y){
					unlist(strsplit(y, split=""))[1]
				})), collapse="")
				paste(c(x$family, yy), collapse=" ")
			}))
			if(length(xres)>1){
				xres[-length(xres)]=paste(xres[-length(xres)], c(rep(",", length(xres)-2), " &"), sep="")
			}
			xres=paste(xres, collapse=" ")
			if(xx[[i]]$bibtype!="Book"){
				xres=paste(c(xres, xx[[i]]$year, xx[[i]]$title, xx[[i]]$note), collapse=". ")
			}else{
				b.title=xx[[i]]$title
				if(nchar(xx[[i]]$edition)>0){b.title=paste(c(b.title, paste(c(xx[[i]]$edition, "edition"), collapse=" ")), collapse=". ")}
				xres=paste(c(xres, xx[[i]]$year, b.title, paste(c(xx[[i]]$publisher, xx[[i]]$address), collapse=", ")), collapse=". ")
			}
			all.res[i]=paste(c(xres, "."), collapse="")
		}
  }else{
    all.res=paste(c(xx$author$given, xx$year, xx$title, xx$organization, xx$address), collapse=". ")
		all.res=paste(c(all.res, "."), collapse="")
  }
  return(all.res)
}

lmer.n<-function(model.res, rev=T){
  #determines sample size for result of call or lmer (total and number of levels per random effect
  #argument is an object resulting from call of lmer or glmer
  #value is a character vector of length one
  if(class(model.res)=="lmerMod"){
    xnames=c("total", names(ranef(model.res)))
    xvals=c(length(residuals(model.res)), unlist(lapply(ranef(model.res), nrow)))
  }else{
    ires=summary(model.res)$ngrps
    xnames=c("total", names(ires))
    xvals=c(length(residuals(model.res)), ires)
  }
  if(rev){
    paste(paste(xvals, xnames, sep=" "), collapse="; ")
  }else{
    paste(paste(xnames, xvals, sep=": "), collapse="; ")
  }
}

wt.txt<-function(x){##wrapper for write.table to clipboard
  #writes results of likelihood ratio test in convenient output (no row or column names, no quotes around text
	if(grepl(sessionInfo()$platform, pattern="linux")){
		con <- pipe("xclip -selection clipboard -i", open="w")
		write.table(x, con, , col.names=F, row.names=F, sep="\t", quote=F)
		close(con)
	}else if(grepl(sessionInfo()$platform, pattern="apple")){
		clip=pipe("pbcopy", "w")
		write.table(x=x, file=clip, col.names=F, row.names=F, sep="\t", quote=F)
		close(clip)
	}else{
		write.table(x=x, file="clipboard", col.names=F, row.names=F, sep="\t", quote=F)
	}
}

wt2<-function(x, row.names=T){#wrapper for write.table to clipboard, row.names optional
	if(grepl(sessionInfo()$platform, pattern="linux")){
		con <- pipe("xclip -selection clipboard -i", open="w")
		write.table(x, con, col.names=T, row.names=row.names, sep="\t", quote=F)
		close(con)
	}else if(grepl(sessionInfo()$platform, pattern="apple")){
		clip=pipe("pbcopy", "w")
		write.table(x=x, file=clip, col.names=F, row.names=F, sep="\t", quote=F)
		close(clip)
	}else{
		write.table(x=x, file="clipboard", col.names=T, row.names=row.names, sep="\t", quote=F)
	}
}

rt2<-function(header=T){#wrapper for read.table from clipboard, header optional
	if(grepl(sessionInfo()$platform, pattern="linux")){
		
	}else if(grepl(sessionInfo()$platform, pattern="apple")){
		read.table(file=pipe("pbpaste"), header=header, sep="\t", fill=T)
	}else{
		read.table(file="clipboard", header=header, sep="\t", fill=T)
	}
}

wp<-function(type="wmf"){
  #wrapper for savePlot with file = "clipboard" and type = "wmf" (default; others are possible)
  savePlot(file="clipboard", type=type)
}

overdisp.correction<-function(coeffs, disp.param){
	coeffs[, "Std. Error"]=coeffs[, "Std. Error"]*sqrt(disp.param)
	coeffs[, "z value"]=coeffs[, "Estimate"]/coeffs[, "Std. Error"]
	coeffs[, "Pr(>|z|)"]=2*pnorm(q=-abs(coeffs[, "z value"]), mean=0, sd=1)
	return(coeffs)
}
