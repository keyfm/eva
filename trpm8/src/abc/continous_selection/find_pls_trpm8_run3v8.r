args <- commandArgs(trailingOnly = TRUE)
c <- args[1]
rid <- args[2]
v <- args[3]
mod <- args[4]

## open File
numComp <- 5; # this flag specifies the components...we stick to 5!

directory <- paste("/mnt/scratch/felix/trpm8/abc/run",rid,"_out/",sep="")
filename1 <- paste("all_sumstat_5ParamBT_",mod,"_",c,"_sdn.tsv.gz",sep="");

filename2 <- paste("all_sumstat_5ParamBT_",mod,"_",c,"_ssv.tsv.gz",sep="");

filename3 <- paste("all_sumstat_5ParamBT_",mod,"_",c,"_ntr.tsv.gz",sep="");

output.dir <- paste("/mnt/scratch/felix/trpm8/abc/analysis/run",rid,"/version",v,"/",sep="")
output.plsda <- paste(output.dir,"find_pls_",c,"_allMod_v",v,".plsda",sep="")

nr=100000 # upscaled this from 20k

model <- as.factor(rep(0:2,each=nr))

# get columns: zcat /mnt/scratch/felix/trpm8/abc/run3_out/all_sumstat_5ParamBT_eur_sdn.tsv.gz |sed '1!d' |tr '\t' '\n' |nl
parCols <- c(1:5)
statCols <- c(8,12,13,17,18,21,26,27,31,32,33,34,35,36)

    ##  8	XPEHH
    ## 12	wdwNonTrg_p1_TD
    ## 13	wdwNonTrg_p1_FWH
    ## 17	wdwNonTrg_p2_TD
    ## 18	wdwNonTrg_p2_FWH
    ## 21	wdwNonTrg_FstAll
    ## 26	wdw2_p1_TD
    ## 27	wdw2_p1_FWH
    ## 31	wdw2_p2_TD
    ## 32	wdw2_p2_FWH
    ## 33	wdw2_freq_selsite_P1
    ## 34	wdw2_freq_selsite_P2
    ## 35	wdw2_FstAll
    ## 36	wdw2_Fst_selsite


###############################################################################################
# for BoxCoxLambda Rplots being in the correct folder!
setwd(output.dir)

#no changes should be required past here
#read file
sdn <- read.table(gzfile(paste(directory, filename1, sep="")), header=T, nrows=nr, skip=0); # aka 'a'
ssv <- read.table(gzfile(paste(directory, filename2, sep="")), header=T, nrows=nr, skip=0); # aka 'b'
ntr <- read.table(gzfile(paste(directory, filename3, sep="")), header=T, nrows=nr, skip=0); # aka 'b2'
dt <- rbind(sdn,ssv,ntr)
stats <- dt[,statCols]
params <- dt[,parCols]
 
#force stats in [1,2]
myMax<-c(); myMin<-c(); lambda<-c(); myGM<-c();
for(i in 1:length(stats)){
	myMax<-c(myMax, max(stats[,i]));
	myMin<-c(myMin, min(stats[,i]));
	stats[,i]<-1+(stats[,i]-myMin[i])/(myMax[i]-myMin[i]);
}

#transform statistics via boxcox  
library("MASS");	
for(i in 1:length(stats)){		
	d<-cbind(stats[,i], params);
	mylm<-lm(as.formula(d), data=d)			
	myboxcox<-boxcox(mylm, lambda=seq(-20, 20, 1), plotit=T, interp=T, eps=1/50);	
	lambda<-c(lambda, myboxcox$x[myboxcox$y==max(myboxcox$y)]);			
	print(paste(names(stats)[i], myboxcox$x[myboxcox$y==max(myboxcox$y)]));
	myGM<-c(myGM, exp(mean(log(stats[,i]))));			
}

#standardize the BC-stats
myBCMeans<-c(); myBCSDs<-c();
for(i in 1:length(stats)){
	stats[,i]<-(stats[,i]^lambda[i] - 1)/(lambda[i]*myGM[i]^(lambda[i]-1));	
	myBCSDs<-c(myBCSDs, sd(stats[,i]));
	myBCMeans<-c(myBCMeans, mean(stats[,i]));		
	stats[,i]<-(stats[,i]-myBCMeans[i])/myBCSDs[i];
}

#perform pls
library("pls");
require("mixOmics")

myPlsr<-plsr(as.matrix(params) ~ as.matrix(stats), scale=F, ncomp=numComp);
k<-plsda(stats,model,ncomp=numComp)
myPlsrDataFrame<-data.frame(k$loadings$X)
#write pls to a file
#myPlsrDataFrame<-data.frame(comp1=myPlsr$loadings[,1]);
#for(i in 2:numComp) { myPlsrDataFrame<-cbind(myPlsrDataFrame, myPlsr$loadings[,i]); } 
kp = names(stats) %in% dimnames(myPlsrDataFrame)[[1]]
write.table(cbind(names(stats)[kp], myMax[kp], myMin[kp], lambda[kp], myGM[kp], myBCMeans[kp], myBCSDs[kp], myPlsrDataFrame), file=output.plsda, col.names=F, row.names=F, sep="\t", quote=F);


### make RMSE plot
pdf(paste(output.dir, "RMSE_plot_",c,".pdf",sep=""));
plot(RMSEP(myPlsr));
dev.off();



#obsa<-read.table("/mnt/uni/ABC/arvalis/arvalis_both.obs", header=T);
#n<-data.frame(a=1:length(names(obsa)), n=names(obsa));
#pdf(paste(directory, "stats_", filename, ".pdf", sep=""), width=9, height=12);
#par(mfrow=c(5,4), cex=0.5)		
#	for(i in c(1:13,25,26,49:51,63,64,76:80,183:227)){
#	plot(density(stats[,i]), xlim=c(min(stats[,i])-max(stats[,i])+min(stats[,i]),max(stats[,i])+max(stats[,i])-min(stats[,i])), main=names(stats)[i]);	
#	print(paste(n[n[,2]==names(stats)[i],1], obsa[n[n[,2]==names(stats)[i],1]]));
#	lines(c(obsa[,n[n[,2]==names(stats)[i],1]], obsa[,n[n[,2]==names(stats)[i],1]]), c(0,1000), col="red")
#}

#dev.off();

