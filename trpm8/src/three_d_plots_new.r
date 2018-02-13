#functions for displaying interacting effects of 2 to 4 predictors on a single response using 3D plots
#written by Roger Mundry
#last updated Apr 22, 2015

#All these functions are wrappers for persp and show 3D-plots.
#In addition to the surface modeled, they show the response, generally averaged per cell of the surface depicted.

#draw.2.w.int.bw.cov.2 and draw.3.w.int.bw.cov depict interactions between two or three covariates, respectively,
#draw.2.w.int.bw.cov.and.facs depicts interactions between two covariates and two factors, and
#draw.2.w.int.bw.cov.circ depicts an interaction between a covariate and a circular variable.

#There is an aditional function named 'use.mat'; you don't need to worry about it, but it might be needed by the other functions.
#For some of the functions to work the package gtools needs to be installed.

#All four functions accept models up to the maximum complexity, but also much simpler models
#(e.g., draw.2.w.int.bw.cov.and.facs would draw the plot based on the simplest model
#including only the main effects, for the most complex model including all interactions
#up the forth-way interaction, and also all models with intermediate complexity).
#All can handle squared covariates and also interactions involving squared covariates.
  #If squared terms are involved these have to enter the model as 'I(covariate^2)'.
  #Similarly, if the function draw.2.w.int.bw.cov.circ is used the model must have been fitted
  #with the circular variable being entered in the model as sin(var)+cos(var).

#All can handle models including other terms (which would simply be ignored).

#All functions display the various figures with the same limits of all axes.
#Points depicting the average response per cell of the fitted model/surface can be scaled
#according to the number of data in the respective cell and are depicted as filled when the
#average is above the fitted model and as open points when they are below. For better visual
#orientation points are connected to the bottom of the figure by a line which is dashed between
#the bottom of the figure and the fitted surface, and solid above it.

#The use of the functions is very similar since most of their arguments are shared among them
#the arguments are the following:

#plot.data: data frame with the predictors (names must be the same as indicated in vars);
#vars: character string with the names of the covariate (quantitative predictors) as used in the model and existing in plot.data;
#var.3.breaks: character. Argument of the function draw.3.w.int.bw.cov, determining how the third covariate is binned.
  #Currently only one value ("equalnumbers") is permitted which means that roughly equal numbers of values
  #are put into each of the bins
#coefs: named numeric vector with the coefficients revealed by the model.
#link=c("identity","logit","log"): character. Link fuction used in the model fitted.
#grid.resol=11: integer. Deteremines the number of lines used to depict the fitted surface
  #(usually the default looks quite good)
#var.names: character. names of the covariates as they should appear on the axis-labels
  #The sequence in which the variables are indicated must be the same as in vars.
  #The function draw.2.w.int.bw.cov.and.facs uses their names as indicated in covariates as default)
#zlab: character. z-axis label.
#zlim: numeric vector of length 2. Determines the limits of the z-axis. The functions draw.2.w.int.bw.cov.and.facs
  #and draw.2.w.int.bw.cov.circ have a default of NA in which case the limits of the z-axis are chosen such that
  #the fitted values and the averaged response can be entirely shown.
#theta: numeric. Horizontal turning angle of the plot(s).
#phi: numeric. Vertical turning angle of the plot(s). Default is 10.
#expand: numeric. Vertical expansion of the plot(s). Has default of 0.6.
#response: numeric. Vector with the response variable.
#size.fac: numeric or NA (the default). Determines whether the volume of the points should be proportional
  #to the number of data points per cell in the depicted surface. If NA all points will be shown with the same size.
  #If not, determines the size of the dots (decrease size.fac when the points are to large, increase when they are too small).
#print.NA.cells: logical. Determines whether cells of the surface where no data do exist should be depicted
  #(T, the default in some of the functions) or not
#covariates: character (only used by function draw.2.w.int.bw.cov.and.facs). Names of the two covariates.
#factors: character (only used by function draw.2.w.int.bw.cov.and.facs). Names of the two factors.
#connect.lines.to: character. Applies only to function draw.2.w.int.bw.cov.and.facs where it
  #determines whether the vertical lines drawn vor visual orientation connect the points to the
  #fitted surface (value "surface") or the bottom of the figure (value: "bottom", the default).
#circular: character. Name of the circular variable as used in the model and in plot.data (applies only for the function draw.2.w.int.bw.cov.circ).

#note that the functions currently are ridiculously slow for larger data sets (will be improved at some point in the future)
#also, currently no detailed ticks at the axes are provided (does someone need that?)
#don't judge the code for elegance (it was written quickly to solve the task)

###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

draw.3.w.int.bw.cov<-function(plot.data, vars, var.3.breaks="equalnumbers", coefs, link=c("identity", "logit", "log"),
  grid.resol=11, var.names, zlab, zlim=NULL, theta, phi=10, r=sqrt(3), expand=0.6, response, size.fac=NA, print.NA.cells, quiet=T){
  #version nov 25, 2013
  old.par = par(no.readonly = TRUE)
  link=link[1]
	if(length(grid.resol)==1){grid.resol=rep(grid.resol, 2)}
  #extract coefficients needed:
	tk=unlist(lapply(strsplit(names(coefs), split=":", fixed=T), function(cf){
		length(cf)==sum(unlist(lapply(vars, function(v){
			grepl(x=cf, pattern=v)
		})))
	}))
	if(names(coefs)[1]=="(Intercept)"){tk[1]=T}
	coefs=coefs[tk]
	
  xvar1=seq(min(plot.data[,vars[1]]), max(plot.data[,vars[1]]), length.out=grid.resol[1])
  yvar2=seq(min(plot.data[,vars[2]]), max(plot.data[,vars[2]]), length.out=grid.resol[2])
	bin.x=cut(x=plot.data[,vars[1]], breaks=xvar1, labels=F, include.lowest=T)
	bin.y=cut(x=plot.data[,vars[2]], breaks=yvar2, labels=F, include.lowest=T)
	bin.x=min(xvar1)+diff(xvar1)[1]/2+((bin.x-1)*diff(xvar1)[1])
	bin.y=min(yvar2)+diff(yvar2)[1]/2+((bin.y-1)*diff(yvar2)[1])

  if(var.3.breaks[1]=="equalnumbers"){
    var.3.breaks=sort(plot.data[,vars[3]])[round(nrow(plot.data)*c(1/3, 2/3), digits=0)]
  }
  
  z.centers=rep(3, nrow(plot.data))
  z.centers[plot.data[, vars[3]]<=var.3.breaks[2]]=2
  z.centers[plot.data[, vars[3]]<=var.3.breaks[1]]=1
  var3.mid.pts=round(tapply(plot.data[, vars[3]], z.centers, mean, na.rm=T), 10)
  bin.z=var3.mid.pts[match(z.centers, as.numeric(names(var3.mid.pts)))]

  obs.mean=tapply(response, list(bin.x, bin.y, bin.z), mean, na.rm=T)
  obs.mean=as.data.frame(as.table(obs.mean))
  names(obs.mean)[1:3]=c("x", "y", "z")
  obs.mean$x=as.numeric(as.character(obs.mean$x))
  obs.mean$y=as.numeric(as.character(obs.mean$y))
  obs.mean$z=as.numeric(as.character(obs.mean$z))
  obs.mean=na.omit(obs.mean)
  #get the predicted means at the centers of the cells:
	rr=runif(nrow(obs.mean))
	xx=obs.mean[, c("x", "y", "z")]
	colnames(xx)=vars 
  pvs=model.matrix(object=as.formula(paste(c("rr", paste(c(1, names(coefs)[-1]), collapse="+")), collapse="~")), data=xx)
	pred.mean=apply(t(coefs*t(pvs[, names(coefs)])), 1, sum)
  
	if(link=="log"){pred.mean=exp(pred.mean)}
  if(link=="logit"){pred.mean=exp(pred.mean)/(1+exp(pred.mean))}
  #determine complete observations:
	complete.obs=apply(is.na(data.frame(plot.data[, vars], response)), 1, sum)==0
  symbol=rep(1, nrow(obs.mean))
  symbol[obs.mean$Freq>pred.mean]=19
  if(!is.na(size.fac)){
    N=as.data.frame(table(bin.x[complete.obs], bin.y[complete.obs], bin.z[complete.obs]))
    N=subset(N, Freq>0)[,4]
  }else{
    N=rep(1, nrow(pvs))
    size.fac=1
  }
    
  layout(cbind(4, 3:1), widths=c(1,5))
  par(mar=rep(1,4))
	pv.mat=cbind(data.frame(expand.grid(xvar1, yvar2)), 0)
	names(pv.mat)=vars
	rr=runif(nrow(pv.mat))
  #browser()
	all.z=lapply(as.vector(var3.mid.pts), function(zz){
    pv.mat[, vars[3]]=zz
		z=model.matrix(object=as.formula(paste(c("rr", paste(c(1, names(coefs)[-1]), collapse="+")), collapse="~")), data=pv.mat)
		z=apply(t(coefs*t(z[, names(coefs)])), 1, sum)
		return(tapply(z, pv.mat[, c(vars[1], vars[2])], mean))
	})
	if(length(zlim)==0){zlim=range(c(unlist(all.z), obs.mean$z), na.rm=T)}
  for(zz in 1:length(as.vector(var3.mid.pts))){
    z=all.z[[zz]]
    if(link=="log"){z=exp(z)}
    if(link=="logit"){z=exp(z)/(1+exp(z))}
    if(!print.NA.cells){
      obs.mat=tapply(response[bin.z==var3.mid.pts[zz]], list(bin.x[bin.z==var3.mid.pts[zz]], bin.y[bin.z==var3.mid.pts[zz]]), mean, na.rm=T)
      rownames(z)=xvar1
      colnames(z)=yvar2
      use=use.mat(obs.mat=obs.mat, z.mat=z)
      z[!use]=NA
    }
		#browser()
    xplot=persp(x=xvar1, y=yvar2, z=z, theta=theta, phi=phi, expand=expand, r=r, xlab=var.names[1], ylab=var.names[2], zlab=zlab, zlim=zlim)
    for(i in 2:(grid.resol[1]-1)){
      lines(trans3d(x=xvar1[i], y=range(yvar2), z=zlim[1], pmat=xplot), lty=3, col="grey")
    }
    for(i in 2:(grid.resol[2]-1)){
      lines(trans3d(x=range(xvar1), y=yvar2[i], z=zlim[1], pmat=xplot), lty=3, col="grey")
    }
    sel.obs=subset(data.frame(obs.mean, pred.mean), obs.mean$z==var3.mid.pts[zz] & Freq<=pred.mean)
    if(nrow(sel.obs)>0){
      for(i in 1:nrow(sel.obs)){
        lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(sel.obs$Freq[i], zlim[1]), pmat=xplot), lty=2)
      }
    }
    sel.obs=subset(data.frame(obs.mean, pred.mean), obs.mean$z==var3.mid.pts[zz] & Freq>pred.mean)
    #browser()
    if(nrow(sel.obs)>0){
      for(i in 1:nrow(sel.obs)){
        lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(zlim[1], sel.obs$pred.mean[i]), pmat=xplot), lty=2)
        lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(sel.obs$pred.mean[i], sel.obs$Freq[i]), pmat=xplot), lty=1)
      }
    }
    points(trans3d(x=obs.mean$x[obs.mean$z==var3.mid.pts[zz]], y=obs.mean$y[obs.mean$z==var3.mid.pts[zz]], 
      z=obs.mean$Freq[obs.mean$z==var3.mid.pts[zz]], pmat=xplot), pch=symbol[obs.mean$z==var3.mid.pts[zz]],
      cex=size.fac*N[obs.mean$z==var3.mid.pts[zz]]^(1/3))
  }
  plot(1,1, axes=F, bty="n", xlab="", ylab="", type="n", xlim=c(0,2), ylim=c(0,2))
  text(var.names[3], x=1.5, y=1, srt=90)
  arrows(x0=1.8, x1=1.8, y0=0.2, y1=1.8, length=0.1)
  par(old.par)
  if(!quiet){
		names(obs.mean)=c(vars, "mean.response")
		obs.mean=data.frame(obs.mean, N=N)
		return(list(zlim=zlim, plotted.data=obs.mean))
	}
}
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

draw.2.w.int.bw.cov.2<-function(plot.data, vars, coefs, link=c("identity", "logit", "log"),
  grid.resol=11, var.names, zlab, zlim=NULL, theta, phi=10, expand=0.6, r=sqrt(3), cex.lab=1, response, size.fac=NA, print.NA.cells, quiet=T, col="black"){
  #version nov 25, 2013
  #plots also when there is no interaction and also squared terms (potentially involved in interactions)
  old.par = par(no.readonly = TRUE)
  link=link[1]
	if(length(grid.resol)==1){grid.resol=rep(grid.resol, 2)}
	#extract coefficients needed:
	tk=unlist(lapply(strsplit(names(coefs), split=":", fixed=T), function(cf){
		length(cf)==sum(unlist(lapply(vars, function(v){
			grepl(x=cf, pattern=v)
		})))
	}))
	if(names(coefs)[1]=="(Intercept)"){tk[1]=T}
	coefs=coefs[tk]

	xvar1=seq(min(plot.data[,vars[1]]), max(plot.data[,vars[1]]), length.out=grid.resol[1])
  yvar2=seq(min(plot.data[,vars[2]]), max(plot.data[,vars[2]]), length.out=grid.resol[2])
	bin.x=cut(x=plot.data[,vars[1]], breaks=xvar1, labels=F, include.lowest=T)
	bin.y=cut(x=plot.data[,vars[2]], breaks=yvar2, labels=F, include.lowest=T)
	bin.x=min(xvar1)+diff(xvar1)[1]/2+((bin.x-1)*diff(xvar1)[1])
	bin.y=min(yvar2)+diff(yvar2)[1]/2+((bin.y-1)*diff(yvar2)[1])
	if(length(dim(response))==2){
		response=response[1,]/apply(response, 1, sum)
	}
  obs.mean=tapply(response, list(bin.x, bin.y), mean, na.rm=T)
  obs.mean=as.data.frame(as.table(obs.mean))
  names(obs.mean)[1:2]=c("x", "y")
  obs.mean$x=as.numeric(as.character(obs.mean$x))
  obs.mean$y=as.numeric(as.character(obs.mean$y))
  obs.mean=na.omit(obs.mean)
	#determine complete observation:
	complete.obs=apply(is.na(data.frame(plot.data[, vars], response)), 1, sum)==0
  #and the sample zize per cell
  obs.n=tapply(complete.obs, list(bin.x, bin.y), sum)
  obs.n=as.data.frame(as.table(obs.n))
  names(obs.n)[1:2]=c("x", "y")
  obs.n$x=as.numeric(as.character(obs.n$x))
  obs.n$y=as.numeric(as.character(obs.n$y))
  obs.n=na.omit(obs.n)
  if(is.na(size.fac)){obs.n$Freq=1; size.fac=1}
  #get the predicted means at the centers of the cells:
	rr=runif(nrow(obs.mean))
	xx=obs.mean[, c("x", "y")]
	colnames(xx)=vars 
  pvs=model.matrix(object=as.formula(paste(c("rr", paste(c(1, names(coefs)[-1]), collapse="+")), collapse="~")), data=xx)
	pred.mean=apply(t(coefs*t(pvs[, names(coefs)])), 1, sum)
  if(link=="log"){pred.mean=exp(pred.mean)}
  if(link=="logit"){pred.mean=exp(pred.mean)/(1+exp(pred.mean))}
  symbol=rep(1, nrow(obs.mean))
  symbol[obs.mean$Freq>pred.mean]=19
	xcol=rep("black", nrow(obs.mean))
	xcol[obs.mean$Freq>pred.mean]=col
	pv.mat=data.frame(expand.grid(xvar1, yvar2))
	names(pv.mat)=vars
	rr=runif(nrow(pv.mat))
	z=model.matrix(object=as.formula(paste(c("rr", paste(c(1, names(coefs)[-1]), collapse="+")), collapse="~")), data=pv.mat)
	z=apply(t(coefs*t(z[, names(coefs)])), 1, sum)
	z=tapply(z, pv.mat[, c(vars[1], vars[2])], mean)

  if(link=="log"){z=exp(z)}
  if(link=="logit"){z=exp(z)/(1+exp(z))}
  if(!print.NA.cells){
    obs.mat=tapply(response, list(bin.x, bin.y), mean, na.rm=T)
    rownames(z)=xvar1
    colnames(z)=yvar2
    use=use.mat(obs.mat=obs.mat, z.mat=z)
    z[!use]=NA
  }
	if(length(zlim)==0){zlim=range(c(pred.mean, obs.mean$Freq), na.rm=T)}
  xplot=persp(x=xvar1, y=yvar2, z=z, theta=theta, phi=phi, expand=expand, r=r, xlab=var.names[1], ylab=var.names[2], zlab=zlab, zlim=zlim, cex.lab=cex.lab)
  for(i in 2:(grid.resol[1]-1)){
    lines(trans3d(x=xvar1[i], y=range(yvar2), z=min(zlim), pmat=xplot), lty=3, col="grey")
	}
  for(i in 2:(grid.resol[2]-1)){
    lines(trans3d(x=range(xvar1), y=yvar2[i], z=min(zlim), pmat=xplot), lty=3, col="grey")
  }
  points(trans3d(x=obs.mean$x, y=obs.mean$y, z=obs.mean$Freq, pmat=xplot), pch=symbol, cex=size.fac*obs.n$Freq^(1/3), col=xcol)
  sel.obs=data.frame(obs.mean, pred.mean)
  for(i in 1:nrow(sel.obs)){
    if(obs.mean$Freq[i]<=pred.mean[i]){
      lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(min(zlim), sel.obs$Freq[i]), pmat=xplot), lty=3)
    }else{
      lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(pred.mean[i], sel.obs$Freq[i]), pmat=xplot), lty=1)
      lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(min(zlim), pred.mean[i]), pmat=xplot), lty=3)
    }
  }
  par(old.par)
  if(!quiet){
		names(obs.mean)=c(vars, "mean.response")
		obs.mean=data.frame(obs.mean, N=obs.n$Freq)
		return(list(zlim=zlim, plotted.data=obs.mean))
	}
}


###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
draw.2.w.int.bw.cov.and.facs<-function(
  plot.data, covariates, factors, coefs, link=c("identity", "logit", "log"), 
  grid.resol=11, var.names=covariates, zlab, zlim=NULL, theta, phi=10, expand=0.6, r=sqrt(3), response,
  connect.lines.to=c("bottom", "surface"), size.fac=NA, print.NA.cells=T, quiet=T, besides=F){
  #version nov 25, 2013
  old.par = par(no.readonly = TRUE)
  if(!is.factor(plot.data[,factors[1]])){
    stop(paste(c(factors[1], "is not a factor"), collapse=" "))
  }
  if(!is.factor(plot.data[,factors[2]])){
    stop(paste(c(factors[2], "is not a factor"), collapse=" "))
  }
  if(min(table(plot.data[,factors[1]]))==0){
    stop(paste(c(factors[1], "has levels not present in the data frame"), collapse=" "))
  }
  if(min(table(plot.data[,factors[2]]))==0){
    stop(paste(c(factors[2], "has levels not present in the data frame"), collapse=" "))
  }
	#extract coefficients needed:
	tk=unlist(lapply(strsplit(names(coefs), split=":", fixed=T), function(cf){
		length(cf)==sum(unlist(lapply(c(covariates, factors), function(v){
			grepl(x=cf, pattern=v)
		})))
	}))
	if(names(coefs)[1]=="(Intercept)"){tk[1]=T}
	coefs=coefs[tk]

	if(length(grid.resol)==1){grid.resol=rep(grid.resol, 2)}
  link=link[1]
  connect.lines.to=connect.lines.to[1]
  
  #reverse the order of the factors in case that with more levels occurs second
  if(length(levels(plot.data[, factors[2]]))>length(levels(plot.data[, factors[1]]))){  
    factors=rev(factors)#rearragne factors such that the one with more levels occurs first
  }

  #create data frame storing the observed mean values per combination of the values 
  #of the grided covariates and factor levels in the data
  #create vectors denoting the grid cell center values for the two covariates:
  xvar1=seq(min(plot.data[,covariates[1]]), max(plot.data[,covariates[1]]), length.out=grid.resol[1])
  yvar2=seq(min(plot.data[,covariates[2]]), max(plot.data[,covariates[2]]), length.out=grid.resol[2])
	bin.x=cut(x=plot.data[,covariates[1]], breaks=xvar1, labels=F, include.lowest=T)
	bin.y=cut(x=plot.data[,covariates[2]], breaks=yvar2, labels=F, include.lowest=T)
	bin.x=min(xvar1)+diff(xvar1)[1]/2+((bin.x-1)*diff(xvar1)[1])
	bin.y=min(yvar2)+diff(yvar2)[1]/2+((bin.y-1)*diff(yvar2)[1])

  #derive the average response per cell of the grid:
  obs.mean=tapply(response, list(bin.x, bin.y, plot.data[, factors[1]], plot.data[, factors[2]]), mean, na.rm=T)
  #... and turn the result into a data frame:
  obs.mean=as.data.frame(as.table(obs.mean))
  names(obs.mean)[1:4]=c(covariates, factors)
  obs.mean[,covariates[1]]=as.numeric(as.character(obs.mean[,covariates[1]]))
  obs.mean[,covariates[2]]=as.numeric(as.character(obs.mean[,covariates[2]]))
  obs.mean=data.frame(na.omit(obs.mean))
  #add columns for the squared terms:
  obs.mean=cbind(obs.mean, obs.mean[, covariates[1]]^2, obs.mean[, covariates[2]]^2)
  colnames(obs.mean)[(ncol(obs.mean)-1):ncol(obs.mean)]=paste("I(", covariates, "^2)", sep="")
	#determine complete observation:
	complete.obs=apply(is.na(data.frame(plot.data[, c(covariates, factors)], response)), 1, sum)==0
  #get the sample size per combination of the factor levels and binned covariates
  if(!is.na(size.fac)){
    N=table(bin.x[complete.obs], bin.y[complete.obs], plot.data[complete.obs, factors[1]], plot.data[complete.obs, factors[2]])
    N=as.data.frame(N)
    N=as.vector(subset(N, Freq>0)$Freq)
  }else{
    N=rep(1, nrow(obs.mean))
    size.fac=1
  }

  #... and columns for the factor levels
  for(ii in 1:2){
    xlev=levels(plot.data[, factors[ii]])
    for(jj in 2:length(xlev)){
      obs.mean=cbind(obs.mean, as.numeric(obs.mean[, factors[ii]]==xlev[jj]))
      colnames(obs.mean)[ncol(obs.mean)]=paste(c(factors[ii], xlev[jj]), collapse="")
    }
  }
  #and finally a column for the intercept:
  obs.mean=cbind(obs.mean, 1)
  colnames(obs.mean)[ncol(obs.mean)]="(Intercept)"
  #determine fitted values per row in obs.mean:
  pred.mean=lapply((1:length(coefs)), function(i.coef){
    names.i.coef=unlist(strsplit(names(coefs[i.coef]), split=":", fixed=T))
    apply(obs.mean[, names.i.coef, drop=F], 1, prod)*coefs[i.coef]
  })
  pred.mean=matrix(unlist(pred.mean), ncol=length(coefs), byrow=F)
  pred.mean=apply(pred.mean, 1, sum)
  if(link=="log"){pred.mean=exp(pred.mean)}
  if(link=="logit"){pred.mean=exp(pred.mean)/(1+exp(pred.mean))}
  #define z-axis limits if they were not indicated by the user:
  if(length(zlim)==0){
    zlim=range(c(pred.mean, obs.mean$Freq, response, na.rm=T))
  }
  #determine symbol for plotting the points:
  symbol=rep(1, nrow(obs.mean))
  symbol[obs.mean$Freq>pred.mean]=19
  
  #define design matrix for plot:
	layout.mat=matrix(1:(length(levels(plot.data[, factors[1]]))*length(levels(plot.data[, factors[2]]))), ncol=length(levels(plot.data[, factors[2]])), byrow=T)
	layout.mat=cbind((max(layout.mat)+1):(max(layout.mat)+length(levels(plot.data[, factors[1]]))), layout.mat)
	layout.mat=rbind(c(max(layout.mat)+length(levels(plot.data[, factors[2]]))+1, (max(layout.mat)+1):(max(layout.mat)+length(levels(plot.data[, factors[2]])))), layout.mat)
  if(!besides){
		layout(layout.mat, heights=c(1, rep(4, length(levels(plot.data[, factors[1]])))), widths=c(1, rep(4, length(levels(plot.data[, factors[2]])))))
	}else{
		layout(t(layout.mat), widths=c(1, rep(4, length(levels(plot.data[, factors[1]])))), heights=c(1, rep(4, length(levels(plot.data[, factors[2]])))))
	}
  #layout.show(max(layout.mat))

  #plot the stuff:
  par(mar=c(0.5,1.5, 0.5, 0.5))
  #loop through combinations of levels of factors:
  for (ii in 1:length(levels(plot.data[, factors[1]]))){
    for(jj in 1:length(levels(plot.data[, factors[2]]))){
      #construct data.frame comprising the intercept (just 1) and the factor terms (dummy vars)
      #with the appropriate values for the particular combination of levels:
      x.pred.data=data.frame(1, t(as.numeric(levels(plot.data[, factors[1]])[-1]==levels(plot.data[, factors[1]])[ii])))
      x.pred.data=cbind(x.pred.data, data.frame(t(as.numeric(levels(plot.data[, factors[2]])[-1]==levels(plot.data[, factors[2]])[jj]))))
      names(x.pred.data)=c("(Intercept)", paste(factors[1], levels(plot.data[, factors[1]])[-1], sep=""), paste(factors[2], levels(plot.data[, factors[2]])[-1], sep=""))
      #get the predicted values considering also the covariates:
      z=outer(xvar1, yvar2, Vectorize(function(x,y){
        #add values of the covariates
        data.to.add=data.frame(x, y, x^2, y^2)
        names(data.to.add)=c(covariates, paste("I(", covariates, "^2)", sep=""))
        data.to.add=cbind(data.to.add, x.pred.data)
        #(data.to.add is a matrix with a single row)
        return(sum(unlist(lapply((1:length(coefs)), function(i.coef){
          names.i.coef=unlist(strsplit(names(coefs[i.coef]), split=":", fixed=T))
          apply(data.to.add[, names.i.coef, drop=F], 1, prod)*coefs[i.coef]
        }))))
      }))
      #transform according to link function
      if(link=="log"){z=exp(z)}
      if(link=="logit"){z=exp(z)/(1+exp(z))}
      #extract data to plot:
      i.plot.data=subset(data.frame(obs.mean, N), obs.mean[, factors[1]]==levels(plot.data[, factors[1]])[ii] & obs.mean[, factors[2]]==levels(plot.data[, factors[2]])[jj])
      if(!print.NA.cells){#still needs to be done
        obs.mat=tapply(i.plot.data$Freq, list(i.plot.data[, covariates[1]], i.plot.data[, covariates[2]]), mean, na.rm=T)
        rownames(z)=xvar1
        colnames(z)=yvar2
        use=use.mat(obs.mat=obs.mat, z.mat=z)
        z[!use]=NA
      }
      #plot the surface:
      xplot=persp(x=xvar1, y=yvar2, z=z, theta=theta, phi=phi, expand=expand, r=r, xlab=var.names[1], ylab=var.names[2], zlab=zlab, zlim=zlim)
      #extract fitted values to plot:
      plot.pred=pred.mean[obs.mean[, factors[1]]==levels(plot.data[, factors[1]])[ii] & obs.mean[, factors[2]]==levels(plot.data[, factors[2]])[jj]]
      plot.symbol=symbol[obs.mean[, factors[1]]==levels(plot.data[, factors[1]])[ii] & obs.mean[, factors[2]]==levels(plot.data[, factors[2]])[jj]]
      #draw the points:
      points(trans3d(x=i.plot.data[, covariates[1]], y=i.plot.data[, covariates[2]], z=i.plot.data$Freq, pmat=xplot), pch=plot.symbol, cex=size.fac*i.plot.data$N^(1/3))
      #and the lines:
      for(kk in 1:nrow(i.plot.data)){
        if(i.plot.data$Freq[kk]<=pred.mean[kk]){
          if(connect.lines.to=="bottom"){
            lines(trans3d(x=i.plot.data[kk, covariates[1]], y=i.plot.data[kk, covariates[2]], z=c(min(zlim), i.plot.data$Freq[kk]), pmat=xplot), lty=3)
          }else{
            lines(trans3d(x=i.plot.data[kk, covariates[1]], y=i.plot.data[kk, covariates[2]], z=c(plot.pred[kk], i.plot.data$Freq[kk]), pmat=xplot), lty=3)
          }
        }else{
          if(connect.lines.to=="bottom"){
            lines(trans3d(x=i.plot.data[kk, covariates[1]], y=i.plot.data[kk, covariates[2]], z=c(plot.pred[kk], i.plot.data$Freq[kk]), pmat=xplot), lty=1)
            lines(trans3d(x=i.plot.data[kk, covariates[1]], y=i.plot.data[kk, covariates[2]], z=c(min(zlim), plot.pred[kk]), pmat=xplot), lty=3)
          }else{
            lines(trans3d(x=i.plot.data[kk, covariates[1]], y=i.plot.data[kk, covariates[2]], z=c(plot.pred[kk], i.plot.data$Freq[kk]), pmat=xplot), lty=1)
          }
        }
      }
      if(connect.lines.to=="bottom"){
        for(kk in 2:(grid.resol[1]-1)){
          lines(trans3d(x=xvar1[kk], y=range(yvar2), z=min(zlim), pmat=xplot), lty=3, col="grey")
          #lines(trans3d(x=range(xvar1), y=yvar2[kk], z=min(zlim), pmat=xplot), lty=3, col="grey")
        }
        for(kk in 2:(grid.resol[2]-1)){
          #lines(trans3d(x=xvar1[kk], y=range(yvar2), z=min(zlim), pmat=xplot), lty=3, col="grey")
          lines(trans3d(x=range(xvar1), y=yvar2[kk], z=min(zlim), pmat=xplot), lty=3, col="grey")
        }
      }
    }
  }
  par(mar=c(0.5,0.5, 0.5, 0.5))
  for(ii in 1:length(levels(plot.data[, factors[1]]))){
    plot(1,1, type="n", axes=F, xlab="", ylab="")    
    text(labels=levels(plot.data[, factors[1]])[ii], x=1, y=1, srt=90, cex=2)
  }
  for(kk in 1:length(levels(plot.data[, factors[2]]))){
    plot(1,1, type="n", axes=F, xlab="", ylab="")    
    text(labels=levels(plot.data[, factors[2]])[kk], x=1, y=1, cex=2)
  }
  par(old.par)
  if(!quiet){
		obs.mean=obs.mean[, 1:5]
		names(obs.mean)=c(covariates, factors, "mean.response")
		obs.mean=data.frame(obs.mean, N=N)
		return(list(zlim=zlim, plotted.data=obs.mean))
	}
}
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

draw.2.w.int.bw.cov.circ<-function(plot.data, vars, circular, coefs, link=c("identity", "logit", "log"),
  grid.resol=11, var.names, zlab, zlim=NULL, theta, phi=10, expand=0.6, r=sqrt(3), response, size.fac=NA, print.NA.cells=T, quiet=T){
  #version nov 25, 2013
  old.par = par(no.readonly = TRUE)
  link=link[1]
  vars=vars[order(as.numeric(vars%in%circular))]
	#extract coefficients needed:
	tk=unlist(lapply(strsplit(names(coefs), split=":", fixed=T), function(cf){
		length(cf)==sum(unlist(lapply(vars, function(v){
			grepl(x=cf, pattern=v)
		})))
	}))
	if(names(coefs)[1]=="(Intercept)"){tk[1]=T}
	coefs=coefs[tk]
	if(length(grid.resol)==1){grid.resol=rep(grid.resol,2)}
  xvar1=seq(min(plot.data[,vars[1]]), max(plot.data[,vars[1]]), length.out=grid.resol[1])
  yvar2=seq(min(plot.data[,vars[2]]), max(plot.data[,vars[2]]), length.out=grid.resol[2])
	bin.x=cut(x=plot.data[,vars[1]], breaks=xvar1, labels=F, include.lowest=T)
	bin.x=min(xvar1)+diff(xvar1)[1]/2+((bin.x-1)*diff(xvar1)[1])
	bin.y=cut(x=plot.data[,vars[2]], breaks=yvar2, labels=F, include.lowest=T)
	bin.y=min(yvar2)+diff(yvar2)[1]/2+((bin.y-1)*diff(yvar2)[1])

  obs.mean=tapply(response, list(bin.x, bin.y), mean, na.rm=T)
  obs.mean=as.data.frame(as.table(obs.mean))
  names(obs.mean)[1:2]=c("x", "y")
  obs.mean$x=as.numeric(as.character(obs.mean$x))
  obs.mean$y=as.numeric(as.character(obs.mean$y))
  obs.mean=na.omit(obs.mean)
	#determine complete observation:
	complete.obs=apply(is.na(data.frame(plot.data[, vars], response)), 1, sum)==0
  if(!is.na(size.fac)){
    sample.size=table(bin.x[complete.obs], bin.y[complete.obs])
    sample.size=as.data.frame(sample.size)
    sample.size=subset(sample.size, Freq>0)
    sample.size=as.vector(sample.size$Freq)
  }else{
    sample.size=rep(1, nrow(obs.mean))
    size.fac=1
  }

  pred.mean=coefs[names(coefs)=="(Intercept)"]+
        coefs[names(coefs)%in%vars[1]]*obs.mean$x+
        coefs[names(coefs)%in%paste(c("sin(", vars[2], ")"), collapse="")]*sin(obs.mean$y)+
        coefs[names(coefs)%in%paste(c("cos(", vars[2], ")"), collapse="")]*cos(obs.mean$y)+
        coefs[names(coefs)%in%paste(c(vars[1], paste(c("sin(", vars[2], ")"), collapse="")), collapse=":") | names(coefs)%in%paste(c(paste(c("sin(", vars[2], ")"), collapse=""), vars[1]), collapse=":")]*obs.mean$x*sin(obs.mean$y)+
        coefs[names(coefs)%in%paste(c(vars[1], paste(c("cos(", vars[2], ")"), collapse="")), collapse=":") | names(coefs)%in%paste(c(paste(c("cos(", vars[2], ")"), collapse=""), vars[1]), collapse=":")]*obs.mean$x*cos(obs.mean$y)
  if(sum(names(coefs)%in%paste(c("I(", vars[1], "^2)"), collapse=""))>0){
        pred.mean=pred.mean+coefs[names(coefs)%in%paste(c("I(", vars[1], "^2)"), collapse="")]*obs.mean$x
  }
  if(sum(names(coefs)%in%paste(c(paste(c("I(", vars[1], "^2)"), collapse=""), paste(c("sin(", vars[2], ")"), collapse="")), collapse=":"))>0 |
     sum(names(coefs)%in%paste(c(paste(c("sin(", vars[2], ")"), collapse=""), paste(c("I(", vars[1], "^2)"), collapse="")), collapse=":"))>0){
        pred.mean=pred.mean+coefs[names(coefs)%in%paste(c(paste(c("I(", vars[1], "^2)"), collapse=""), paste(c("sin(", vars[2], ")"), collapse="")), collapse=":") | names(coefs)%in%paste(c(paste(c("sin(", vars[2], ")"), collapse=""), paste(c("I(", vars[1], "^2)"), collapse="")), collapse=":")]*(obs.mean$x^2)*sin(obs.mean$y)
        pred.mean=pred.mean+coefs[names(coefs)%in%paste(c(paste(c("I(", vars[1], "^2)"), collapse=""), paste(c("cos(", vars[2], ")"), collapse="")), collapse=":") | names(coefs)%in%paste(c(paste(c("cos(", vars[2], ")"), collapse=""), paste(c("I(", vars[1], "^2)"), collapse="")), collapse=":")]*(obs.mean$x^2)*cos(obs.mean$y)
  }

  if(link=="log"){pred.mean=exp(pred.mean)}
  if(link=="logit"){pred.mean=exp(pred.mean)/(1+exp(pred.mean))}
  symbol=rep(1, nrow(obs.mean))
  symbol[obs.mean$Freq>pred.mean]=19
  z=outer(xvar1, yvar2, Vectorize(function(x,y){
    LP=coefs[names(coefs)=="(Intercept)"]+
      coefs[names(coefs)%in%vars[1]]*x+
      coefs[names(coefs)%in%paste(c("sin(", vars[2], ")"), collapse="")]*sin(y)+
      coefs[names(coefs)%in%paste(c("cos(", vars[2], ")"), collapse="")]*cos(y)+
      coefs[names(coefs)%in%paste(c(vars[1], paste(c("sin(", vars[2], ")"), collapse="")), collapse=":") | names(coefs)%in%paste(c(paste(c("sin(", vars[2], ")"), collapse=""), vars[1]), collapse=":")]*x*sin(y)+
      coefs[names(coefs)%in%paste(c(vars[1], paste(c("cos(", vars[2], ")"), collapse="")), collapse=":") | names(coefs)%in%paste(c(paste(c("cos(", vars[2], ")"), collapse=""), vars[1]), collapse=":")]*x*cos(y)

      if(sum(names(coefs)%in%paste(c("I(", vars[1], "^2)"), collapse=""))>0){
        LP=LP+coefs[names(coefs)%in%paste(c("I(", vars[1], "^2)"), collapse="")]*x
      }
      if(sum(names(coefs)%in%paste(c(paste(c("I(", vars[1], "^2)"), collapse=""), paste(c("sin(", vars[2], ")"), collapse="")), collapse=":"))>0 |
         sum(names(coefs)%in%paste(c(paste(c("sin(", vars[2], ")"), collapse=""), paste(c("I(", vars[1], "^2)"), collapse="")), collapse=":"))>0){
            LP=LP+coefs[names(coefs)%in%paste(c(paste(c("I(", vars[1], "^2)"), collapse=""), paste(c("sin(", vars[2], ")"), collapse="")), collapse=":") | names(coefs)%in%paste(c(paste(c("sin(", vars[2], ")"), collapse=""), paste(c("I(", vars[1], "^2)"), collapse="")), collapse=":")]*(x^2)*sin(y)
            LP=LP+coefs[names(coefs)%in%paste(c(paste(c("I(", vars[1], "^2)"), collapse=""), paste(c("cos(", vars[2], ")"), collapse="")), collapse=":") | names(coefs)%in%paste(c(paste(c("cos(", vars[2], ")"), collapse=""), paste(c("I(", vars[1], "^2)"), collapse="")), collapse=":")]*(x^2)*cos(y)
      }

    return(LP)
  }))
  if(link=="log"){z=exp(z)}
  if(link=="logit"){z=exp(z)/(1+exp(z))}
  if(!print.NA.cells){
    obs.mat=tapply(response, list(bin.x, bin.y), mean, na.rm=T)
    rownames(z)=xvar1
    colnames(z)=yvar2
    use=use.mat(obs.mat=obs.mat, z.mat=z)
    z[!use]=NA
  }
  if(length(zlim)==0){
    zlim=range(c(z, obs.mean$Freq), na.rm=T)
  }
  xplot=persp(x=xvar1, y=yvar2, z=z, theta=theta, phi=phi, expand=expand, r=r, xlab=var.names[1],
    ylab=var.names[2], zlab=zlab, zlim=zlim)
  for(i in 2:(grid.resol[1]-1)){
    lines(trans3d(x=xvar1[i], y=range(yvar2), z=min(zlim), pmat=xplot), lty=3, col="grey")
    lines(trans3d(x=range(xvar1), y=yvar2[i], z=min(zlim), pmat=xplot), lty=3, col="grey")
  }
  for(i in 2:(grid.resol[2]-1)){
    lines(trans3d(x=xvar1[i], y=range(yvar2), z=min(zlim), pmat=xplot), lty=3, col="grey")
    lines(trans3d(x=range(xvar1), y=yvar2[i], z=min(zlim), pmat=xplot), lty=3, col="grey")
  }
  sel.obs=data.frame(obs.mean, pred.mean)
  for(i in 1:nrow(sel.obs)){
    if(obs.mean$Freq[i]<=pred.mean[i]){
      lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(min(zlim), sel.obs$Freq[i]), pmat=xplot), lty=3)
    }else{
      lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(pred.mean[i], sel.obs$Freq[i]), pmat=xplot), lty=1)
      lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(min(zlim), pred.mean[i]), pmat=xplot), lty=3)
    }
  }
  points(trans3d(x=obs.mean$x, y=obs.mean$y, z=obs.mean$Freq, pmat=xplot), pch=symbol, cex=size.fac*(sample.size^(1/3)))
  par(old.par)
  if(!quiet){
		names(obs.mean)=c(vars, "mean.response")
		obs.mean=data.frame(obs.mean, N=sample.size)
		return(list(zlim=zlim, plotted.data=obs.mean))
	}
}
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

draw.2.w.int.bw.cov.and.factor<-function(
  plot.data, covariates, factor, coefs, link=c("identity", "logit", "log"),
  grid.resol=11, var.names=covariates, zlab, zlim=NULL, theta, phi=10, expand=0.6, r = sqrt(3), response,
  connect.lines.to=c("bottom", "surface"), besides=T, size.fac=NA, print.NA.cells=T, quiet=T, show.data=T){
  #version dec 30, 2013
  old.par = par(no.readonly = TRUE)
  if(!is.factor(plot.data[,factor])){
    stop(paste(c(factor, "is not a factor"), collapse=" "))
  }
  if(min(table(plot.data[,factor]))==0){
    stop(paste(c(factor, "has levels not present in the data frame"), collapse=" "))
  }
  link=link[1]
  connect.lines.to=connect.lines.to[1]
  #plots also when there is no interaction and also squared terms (potentially involved in interactions)
	if(length(grid.resol)==1){grid.resol=rep(grid.resol, 2)}
	#extract coefficients needed:
	tk=unlist(lapply(strsplit(names(coefs), split=":", fixed=T), function(cf){
		length(cf)==sum(unlist(lapply(c(covariates, factor), function(v){
			grepl(x=cf, pattern=v)
		})))
	}))
	if(names(coefs)[1]=="(Intercept)"){tk[1]=T}
	coefs=coefs[tk]

  #create data frame storing the observed mean values per combination of the values
  #of the grided covariates and factor levels in the data
  #create vectors denoting the grid cell center values for the two covariates:
  xvar1=seq(min(plot.data[,covariates[1]]), max(plot.data[,covariates[1]]), length.out=grid.resol[1])
  yvar2=seq(min(plot.data[,covariates[2]]), max(plot.data[,covariates[2]]), length.out=grid.resol[2])
	bin.x=cut(x=plot.data[,covariates[1]], breaks=xvar1, labels=F, include.lowest=T)
	bin.x=min(xvar1)+diff(xvar1)[1]/2+((bin.x-1)*diff(xvar1)[1])
	bin.y=cut(x=plot.data[,covariates[2]], breaks=yvar2, labels=F, include.lowest=T)
	bin.y=min(yvar2)+diff(yvar2)[1]/2+((bin.y-1)*diff(yvar2)[1])
  #derive the average response per cell of the grid:
  obs.mean=tapply(response, list(bin.x, bin.y, plot.data[, factor]), mean, na.rm=T)
  #... and turn the result into a data frame:
  obs.mean=as.data.frame(as.table(obs.mean))
  names(obs.mean)[1:3]=c(covariates, factor)
  obs.mean[,covariates[1]]=as.numeric(as.character(obs.mean[,covariates[1]]))
  obs.mean[,covariates[2]]=as.numeric(as.character(obs.mean[,covariates[2]]))
  obs.mean=data.frame(na.omit(obs.mean))
  #add columns for the squared terms:
  obs.mean=cbind(obs.mean, obs.mean[, covariates[1]]^2, obs.mean[, covariates[2]]^2)
  colnames(obs.mean)[(ncol(obs.mean)-1):ncol(obs.mean)]=paste("I(", covariates, "^2)", sep="")
	#determine complete observation:
	complete.obs=apply(is.na(data.frame(plot.data[, c(covariates, factor)], response)), 1, sum)==0
  #get the sample size per combination of the factor levels and binned covariates
  if(!is.na(size.fac)){
    N=table(bin.x[complete.obs], bin.y[complete.obs], plot.data[complete.obs, factor])
    N=as.data.frame(N)
    N=as.vector(subset(N, Freq>0)$Freq)
  }else{
    N=rep(1, nrow(obs.mean))
    size.fac=1
  }

  #... and columns for the factor levels
  xlev=levels(plot.data[, factor])
  for(jj in 2:length(xlev)){
    obs.mean=cbind(obs.mean, as.numeric(obs.mean[, factor]==xlev[jj]))
    colnames(obs.mean)[ncol(obs.mean)]=paste(c(factor, xlev[jj]), collapse="")
  }
  #and finally a column for the intercept:
  obs.mean=cbind(obs.mean, 1)
  colnames(obs.mean)[ncol(obs.mean)]="(Intercept)"
  #determine fitted values per row in obs.mean:
  pred.mean=lapply((1:length(coefs)), function(i.coef){
    names.i.coef=unlist(strsplit(names(coefs[i.coef]), split=":", fixed=T))
    apply(obs.mean[, names.i.coef, drop=F], 1, prod)*coefs[i.coef]
  })
  pred.mean=matrix(unlist(pred.mean), ncol=length(coefs), byrow=F)
  pred.mean=apply(pred.mean, 1, sum)
  if(link=="log"){pred.mean=exp(pred.mean)}
  if(link=="logit"){pred.mean=exp(pred.mean)/(1+exp(pred.mean))}
  #define z-axis limits if they were not indicated by the user:
  if(length(zlim)==0){
    zlim=range(c(pred.mean, obs.mean$Freq, response, na.rm=T))
  }

  #determine symbol for plotting the points:
  symbol=rep(1, nrow(obs.mean))
  symbol[obs.mean$Freq>pred.mean]=19
  #define design matrix for plot:
  layout.mat=matrix(1:(length(levels(plot.data[, factor]))), nrow=1)
  layout.mat=rbind((max(layout.mat)+1):(max(layout.mat)+length(levels(plot.data[, factor]))), layout.mat)
  #layout.mat=rbind(c(max(layout.mat)+length(levels(plot.data[, factors[2]]))+1, (max(layout.mat)+1):(max(layout.mat)+length(levels(plot.data[, factors[2]])))), layout.mat)
  if(besides){
    layout(layout.mat, heights=c(1, 6), widths=rep(4, length(levels(plot.data[, factor]))))
  }else{
    layout.mat=t(layout.mat)
    layout(layout.mat, widths=c(1, 6), heights=rep(4, length(levels(plot.data[, factor]))))
  }
  #layout.show(max(layout.mat))

  #plot the stuff:
  par(mar=c(0.5,1.5, 0.5, 0.5))

  #loop through combinations of levels of factors:
  for (ii in 1:length(levels(plot.data[, factor]))){
    #for(jj in 1:length(levels(plot.data[, factors[2]]))){
      #construct data.frame comprising the intercept (just 1) and the factor terms (dummy vars)
      #with the appropriate values for the particular combination of levels:
      x.pred.data=data.frame(1, t(as.numeric(levels(plot.data[, factor])[-1]==levels(plot.data[, factor])[ii])))
      #x.pred.data=cbind(x.pred.data, data.frame(t(as.numeric(levels(plot.data[, factors[2]])[-1]==levels(plot.data[, factors[2]])[jj]))))
      names(x.pred.data)=c("(Intercept)", paste(factor, levels(plot.data[, factor])[-1], sep=""))
      #get the predicted values considering also the covariates:
      z=outer(xvar1, yvar2, Vectorize(function(x,y){
        #add values of the covariates
        data.to.add=data.frame(x, y, x^2, y^2)
        names(data.to.add)=c(covariates, paste("I(", covariates, "^2)", sep=""))
        data.to.add=cbind(data.to.add, x.pred.data)
        #(data.to.add is a matrix with a single row)
        return(sum(unlist(lapply((1:length(coefs)), function(i.coef){
          names.i.coef=unlist(strsplit(names(coefs[i.coef]), split=":", fixed=T))
          apply(data.to.add[, names.i.coef, drop=F], 1, prod)*coefs[i.coef]
        }))))
      }))
      #transform according to link function
      if(link=="log"){z=exp(z)}
      if(link=="logit"){z=exp(z)/(1+exp(z))}
      #extract data to plot:
      i.plot.data=subset(data.frame(obs.mean, N), obs.mean[, factor]==levels(plot.data[, factor])[ii])
      if(!print.NA.cells){#still needs to be done
        obs.mat=tapply(i.plot.data$Freq, list(i.plot.data[, covariates[1]], i.plot.data[, covariates[2]]), mean, na.rm=T)
        rownames(z)=xvar1
        colnames(z)=yvar2
        use=use.mat(obs.mat=obs.mat, z.mat=z)
        z[!use]=NA
      }
      #plot the surface:
      xplot=persp(x=xvar1, y=yvar2, z=z, theta=theta, phi=phi, expand=expand, r=r, xlab=var.names[1], ylab=var.names[2], zlab=zlab, zlim=zlim)
      if(show.data){
				#extract fitted values to plot:
				plot.pred=pred.mean[obs.mean[, factor]==levels(plot.data[, factor])[ii]]
				plot.symbol=symbol[obs.mean[, factor]==levels(plot.data[, factor])[ii]]
				#draw the points:
				points(trans3d(x=i.plot.data[, covariates[1]], y=i.plot.data[, covariates[2]], z=i.plot.data$Freq, pmat=xplot), pch=plot.symbol, cex=size.fac*i.plot.data$N^(1/3))
				#and the lines:
				for(kk in 1:nrow(i.plot.data)){
					if(i.plot.data$Freq[kk]<=pred.mean[kk]){
						if(connect.lines.to=="bottom"){
							lines(trans3d(x=i.plot.data[kk, covariates[1]], y=i.plot.data[kk, covariates[2]], z=c(min(zlim), i.plot.data$Freq[kk]), pmat=xplot), lty=3)
						}else{
							lines(trans3d(x=i.plot.data[kk, covariates[1]], y=i.plot.data[kk, covariates[2]], z=c(plot.pred[kk], i.plot.data$Freq[kk]), pmat=xplot), lty=3)
						}
					}else{
						if(connect.lines.to=="bottom"){
							lines(trans3d(x=i.plot.data[kk, covariates[1]], y=i.plot.data[kk, covariates[2]], z=c(plot.pred[kk], i.plot.data$Freq[kk]), pmat=xplot), lty=1)
							lines(trans3d(x=i.plot.data[kk, covariates[1]], y=i.plot.data[kk, covariates[2]], z=c(min(zlim), plot.pred[kk]), pmat=xplot), lty=3)
						}else{
							lines(trans3d(x=i.plot.data[kk, covariates[1]], y=i.plot.data[kk, covariates[2]], z=c(plot.pred[kk], i.plot.data$Freq[kk]), pmat=xplot), lty=1)
						}
					}
				}
				if(connect.lines.to=="bottom"){
					for(kk in 2:(grid.resol[1]-1)){
						lines(trans3d(x=xvar1[kk], y=range(yvar2), z=min(zlim), pmat=xplot), lty=3, col="grey")
						#lines(trans3d(x=range(xvar1), y=yvar2[kk], z=min(zlim), pmat=xplot), lty=3, col="grey")
					}
				for(kk in 2:(grid.resol[2]-1)){
					#lines(trans3d(x=xvar1[kk], y=range(yvar2), z=min(zlim), pmat=xplot), lty=3, col="grey")
					lines(trans3d(x=range(xvar1), y=yvar2[kk], z=min(zlim), pmat=xplot), lty=3, col="grey")
				}
			}
    }
  }

  par(mar=c(0.5,0.5, 0.5, 0.5))
  for(ii in 1:length(levels(plot.data[, factor]))){
    plot(1,1, type="n", axes=F, xlab="", ylab="")
    text(labels=levels(plot.data[, factor])[ii], x=1, y=1, srt=90*(1-besides), cex=2)
  }
  par(old.par)
  if(!quiet){
		obs.mean=obs.mean[, 1:4]
		names(obs.mean)=c(c(covariates, factor), "mean.response")
		obs.mean=data.frame(obs.mean, N=N)
		return(list(zlim=zlim, plotted.data=obs.mean))
	}
}
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

intplot.3.cov.1.fac<-function(plot.data, covariates, factor, response, coefs, link=c("identity", "logit", "log"),
	var.3.breaks="equalnumbers", grid.resol=11, 
	var.names, zlab, zlim=NULL, level.seq=NULL, level.labels=NULL,
	theta, phi=10, r=sqrt(3), expand=0.6, 
	size.fac=NA, print.NA.cells, quiet=T, show.data=T){
  #version April 18, 2016
  old.par = par(no.readonly = TRUE)
  link=link[1]
  if(length(level.seq)==0){level.seq=levels(droplevels(na.omit(data.frame(plot.data, response)))[, factor])}
  if(length(level.labels)==0){level.labels=level.seq}
	if(length(grid.resol)==1){grid.resol=rep(grid.resol, 2)}
  #extract coefficients needed:
	tk=unlist(lapply(strsplit(names(coefs), split=":", fixed=T), function(cf){
		length(cf)==sum(unlist(lapply(c(covariates, factor), function(v){
			grepl(x=cf, pattern=v)
		})))
	}))
	if(length(intersect(names(coefs)[1], "(Intercept)"))){tk[1]=T}
	coefs=coefs[tk]
	plot.data=droplevels(plot.data)
	
  xvar1=seq(min(plot.data[,covariates[1]]), max(plot.data[,covariates[1]]), length.out=grid.resol[1])
  yvar2=seq(min(plot.data[,covariates[2]]), max(plot.data[,covariates[2]]), length.out=grid.resol[2])
	bin.x=cut(x=plot.data[,covariates[1]], breaks=xvar1, labels=F, include.lowest=T)
	bin.y=cut(x=plot.data[,covariates[2]], breaks=yvar2, labels=F, include.lowest=T)
	bin.x=min(xvar1)+diff(xvar1)[1]/2+((bin.x-1)*diff(xvar1)[1])
	bin.y=min(yvar2)+diff(yvar2)[1]/2+((bin.y-1)*diff(yvar2)[1])

  if(var.3.breaks[1]=="equalnumbers"){
    var.3.breaks=sort(plot.data[,covariates[3]])[round(nrow(plot.data)*c(1/3, 2/3), digits=0)]
  }
  
  z.centers=rep(3, nrow(plot.data))
  z.centers[plot.data[, covariates[3]]<=var.3.breaks[2]]=2
  z.centers[plot.data[, covariates[3]]<=var.3.breaks[1]]=1
  var3.mid.pts=round(tapply(plot.data[, covariates[3]], z.centers, mean, na.rm=T), 10)
  bin.z=var3.mid.pts[match(z.centers, as.numeric(names(var3.mid.pts)))]
  obs.mean=aggregate(response, list(bin.x, bin.y, bin.z, factor=plot.data[, factor]), mean, na.rm=T)
  #obs.mean=as.data.frame(as.table(obs.mean))
  names(obs.mean)=c("x", "y", "z", factor, "Freq")
  #obs.mean$x=as.numeric(as.character(obs.mean$x))
  #obs.mean$y=as.numeric(as.character(obs.mean$y))
  #obs.mean$z=as.numeric(as.character(obs.mean$z))
  #obs.mean=na.omit(obs.mean)
  #get the predicted means at the centers of the cells:
	model=names(coefs)
	model[model%in%"(Intercept)"]=1
	model=strsplit(model, split=":", fixed=T)
	for(ls in levels(plot.data[, factor])[-1]){
		model=lapply(model, function(x){
			x[x==paste(c(factor, ls), collapse="")]=factor
			return(x)
		})
	}
	model=unlist(lapply(model, paste, collapse=":"))
	model=paste(c("rr", paste(model, collapse="+")), collapse="~")
	rr=runif(nrow(obs.mean))
	xx=obs.mean[, c("x", "y", "z", factor)]
	colnames(xx)=c(covariates, factor)
  pvs=model.matrix(object=as.formula(model), data=xx)
	pred.mean=apply(t(coefs*t(pvs[, names(coefs)])), 1, sum)
	if(link=="log"){pred.mean=exp(pred.mean)}
  if(link=="logit"){pred.mean=exp(pred.mean)/(1+exp(pred.mean))}
  #determine complete observations:
	complete.obs=apply(is.na(data.frame(plot.data[, c(covariates, factor)], response)), 1, sum)==0

	symbol=rep(1, nrow(obs.mean))
  symbol[obs.mean$Freq>pred.mean]=19
  if(!is.na(size.fac)){
    N=aggregate(!is.na(response), list(bin.x[complete.obs], bin.y[complete.obs], bin.z[complete.obs], factor=plot.data[complete.obs, factor]), sum)
    N=subset(N, x>0)$x
  }else{
    N=rep(1, nrow(pvs))
    size.fac=1
  }
  #layout  
  xmat=matrix(1:(length(var3.mid.pts)*length(level.seq)), ncol=length(level.seq), byrow=F)
  xmat=cbind(max(xmat)+1, xmat)
  xmat=rbind(max(xmat)+c(length(level.seq)+1, 1:length(level.seq)), xmat)
  layout(xmat, widths=c(1, rep(5, length(level.seq))), heights=c(1, rep(5, length(var3.mid.pts))))
  par(mar=rep(1, 4))

	pv.mat=data.frame(expand.grid(xvar1, yvar2, as.vector(var3.mid.pts), levels(plot.data[, factor])))
	names(pv.mat)=c(covariates, factor)
	rr=runif(nrow(pv.mat))
	z=model.matrix(object=as.formula(model), data=pv.mat)
	pv.mat$z=apply(t(coefs*t(z[, names(coefs)])), 1, sum)
	
	if(length(zlim)==0){zlim=range(c(unlist(pv.mat$z), obs.mean$Freq), na.rm=T)}
	var3.mid.pts=rev(var3.mid.pts)
  for(ls in 1:length(level.seq)){
		for(zz in 1:length(as.vector(var3.mid.pts))){
			iz=subset(pv.mat, pv.mat[, factor]==level.seq[ls] & pv.mat[, covariates[3]]==var3.mid.pts[zz])
			iz=tapply(iz$z, iz[, c(covariates[1], covariates[2])], mean)
			if(link=="log"){z=exp(z)}
			if(link=="logit"){z=exp(z)/(1+exp(z))}
			if(!print.NA.cells){
				obs.mat=tapply(X=response[bin.z==var3.mid.pts[zz] & plot.data[, factor]==level.seq[ls]], 
					INDEX=list(bin.x[bin.z==var3.mid.pts[zz] & plot.data[, factor]==level.seq[ls]], bin.y[bin.z==var3.mid.pts[zz] & plot.data[, factor]==level.seq[ls]]), 
					FUN=mean, na.rm=T)
				use=use.mat(obs.mat=obs.mat, z.mat=iz)
				iz[!use]=NA
			}
			xplot=persp(x=xvar1, y=yvar2, z=iz, theta=theta, phi=phi, expand=expand, r=r, xlab=var.names[1], ylab=var.names[2], zlab=zlab, zlim=zlim)
			if(show.data){
				for(i in 2:(grid.resol[1]-1)){
					lines(trans3d(x=xvar1[i], y=range(yvar2), z=zlim[1], pmat=xplot), lty=3, col="grey")
				}
				for(i in 2:(grid.resol[2]-1)){
					lines(trans3d(x=range(xvar1), y=yvar2[i], z=zlim[1], pmat=xplot), lty=3, col="grey")
				}
				sel.obs=subset(data.frame(obs.mean, pred.mean), obs.mean$z==var3.mid.pts[zz] & obs.mean[, factor]==level.seq[ls] & Freq<=pred.mean)
				if(nrow(sel.obs)>0){
					for(i in 1:nrow(sel.obs)){
						lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(sel.obs$Freq[i], zlim[1]), pmat=xplot), lty=2)
					}
				}
				sel.obs=subset(data.frame(obs.mean, pred.mean), obs.mean$z==var3.mid.pts[zz] & obs.mean[, factor]==level.seq[ls] & Freq>pred.mean)
				if(nrow(sel.obs)>0){
					for(i in 1:nrow(sel.obs)){
						lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(zlim[1], sel.obs$pred.mean[i]), pmat=xplot), lty=2)
						lines(trans3d(x=sel.obs$x[i], y=sel.obs$y[i], z=c(sel.obs$pred.mean[i], sel.obs$Freq[i]), pmat=xplot), lty=1)
					}
				}
				points(trans3d(x=obs.mean$x[obs.mean$z==var3.mid.pts[zz] & obs.mean[, factor]==level.seq[ls]], 
					y=obs.mean$y[obs.mean$z==var3.mid.pts[zz] & obs.mean[, factor]==level.seq[ls]], 
					z=obs.mean$Freq[obs.mean$z==var3.mid.pts[zz] & obs.mean[, factor]==level.seq[ls]], 
					pmat=xplot), pch=symbol[obs.mean$z==var3.mid.pts[zz]],
					cex=size.fac*N[obs.mean$z==var3.mid.pts[zz] & obs.mean[, factor]==level.seq[ls]]^(1/3))
			}
		}
	}
	plot(1,1, axes=F, bty="n", xlab="", ylab="", type="n", xlim=c(0,2), ylim=c(0,2))
  text(var.names[3], x=1.4, y=1, srt=90)
  arrows(x0=1.8, x1=1.8, y0=0.2, y1=1.8, length=0.1)
  for(ls in 1:length(level.seq)){
		plot(1,1, axes=F, bty="n", xlab="", ylab="", type="n", xlim=c(0, 1), ylim=c(0, 1))
		text(level.labels[ls], x=0.5, y=0.5, cex=1.5)
  }
  
  par(old.par)
  if(!quiet){
		names(obs.mean)=c(covariates, "mean.response")
		obs.mean=data.frame(obs.mean, N=N)
		return(list(zlim=zlim, plotted.data=obs.mean))
	}
}
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################
###################################################################################################

use.mat<-function(obs.mat, z.mat){
  #version Mar. 13, 2012
  #create matrix with fake NAs such that each entry in z has four observations around it:
  xobs=cbind(NA, rbind(NA,obs.mat, NA),NA)
  zr.vals=as.numeric(rownames(z.mat))
  zc.vals=as.numeric(colnames(z.mat))
  rownames(xobs)[1]=min(zr.vals)-mean(diff(zr.vals))/2
  rownames(xobs)[nrow(xobs)]=max(zr.vals)+mean(diff(zr.vals))/2
  colnames(xobs)[1]=min(zc.vals)-mean(diff(zc.vals))/2
  colnames(xobs)[ncol(xobs)]=max(zc.vals)+mean(diff(zc.vals))/2
  #potentially include missing rows and columns:
  if(ncol(xobs)<ncol(z.mat)+1){
    cnames=round(seq(min(as.numeric(colnames(xobs))), max(as.numeric(colnames(xobs))), length.out=ncol(z.mat)+1),digits=4)
    xx=cnames%in%round(as.numeric(colnames(xobs)), digits=4)
    new.xobs=matrix(NA, ncol=ncol(z.mat)+1, nrow=nrow(xobs))
    new.xobs[,xx]=xobs
    colnames(new.xobs)=cnames
    rownames(new.xobs)=rownames(xobs)
    xobs=new.xobs
  }
  if(nrow(xobs)<nrow(z.mat)+1){
    cnames=round(seq(min(as.numeric(rownames(xobs))), max(as.numeric(rownames(xobs))), length.out=nrow(z.mat)+1),digits=4)
    xx=cnames%in%round(as.numeric(rownames(xobs)), digits=4)
    new.xobs=matrix(NA, nrow=nrow(z.mat)+1, ncol=ncol(xobs))
    new.xobs[xx,]=xobs
    rownames(new.xobs)=cnames
    colnames(new.xobs)=colnames(xobs)
    xobs=new.xobs
  }
  x.use=matrix(NA, ncol=ncol(z.mat), nrow=nrow(z.mat))
  for(rz in 1:nrow(x.use)){
    for(cz in 1:ncol(x.use)){
      x.use[rz, cz]=sum(!is.na(xobs[rz:(rz+1), cz:(cz+1)]))>0
    }
  }
  return(x.use)
}

