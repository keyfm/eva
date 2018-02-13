args <- commandArgs(trailingOnly = TRUE)
rid <- args[1]
v <- args[2]
pps <- args[3]
## c <- args[2]
library(parallel)

if (pps == 'eur'){
  pops=c('CEU','GBR','TSI','FIN','IBS')
  mc=colorRampPalette(c("lightblue", "darkblue"))(5)
} else if(pps == 'asi'){
  pops=c('CDX','CHB','CHS','KHV','JPT','BEB','GIH','ITU','PJL','STU')
  mc=c(colorRampPalette(c("lightgreen", "darkgreen"))(5),colorRampPalette(c("pink", "purple"))(5))
} else if(pps == 'all'){
  pops=c('CEU','GBR','TSI','FIN','IBS','CDX','CHB','CHS','KHV','JPT','BEB','GIH','ITU','PJL','STU')
  mc=c(colorRampPalette(c("lightblue", "darkblue"))(5),colorRampPalette(c("lightgreen", "darkgreen"))(5),colorRampPalette(c("pink", "purple"))(5))
}

########################
###### Plot Posterior Prob for Model Choice
########################
pdf(paste('/home/felix_schulze/projects/trpm8_proj/pdf/abc/run',rid,'/version',v,'_modelPosterior_',pps,'.pdf',sep=''))
md <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run',rid,'/version',v,'/marginalDensities_ntr_sdn_ssv.txt',sep=""),header=F)
d <- t(md[,2:4])
colnames(d) <- md[,1]
for (i in 1:ncol(d)){
  d[,i] <- d[,i]/sum(d[,i])
}
barplot(d, col=c('white','lightblue','darkgoldenrod1'),main=paste('run',rid,';ntr=white;sdn=blue,ssv=red',sep=''),las=3)
dev.off()

########################
###### Plot Posterior Prob for Model Choice
########################
pdf(paste('/home/felix_schulze/projects/trpm8_proj/pdf/abc/run',rid,'/version',v,'_paramPosterior_',pps,'.pdf',sep=''),height=10,width=15)
par(mfrow=c(2,3))
cont = c('eur') ## note for now cont only used for prior plotting...ASI and EUR no difference! -> Thats why hard coded!
for (mod in c('ssv','sdn')){
  if(mod=='sdn'){prior_cols <- c(1,2,3);post_cols <- c(3,4,5)} else {prior_cols <- c(1,3,5);post_cols <- c(3,4,5)}
  prior <- read.table(paste('/mnt/scratch/felix/trpm8/abc/run',rid,'_out/transf/all_sumstat_5ParamBT_trans_',cont,'_',mod,'_run',rid,'v',v,'_5PLS.tsv',sep=""),header=T,nrows = 100000)
  for (i in 1:3){
    l=density(prior[,prior_cols[i]])
    mymax <- max(l$y)+((max(l$y))*3)
    hist(prior[,prior_cols[i]],freq=FALSE,ylim=c(0,mymax),main=paste('run',rid,'v',v,'; ',mod,' ',colnames(prior)[ prior_cols[i] ],sep=''),xlab=colnames(prior)[i])
    for (pid in 1:length(pops)){
      post <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run',rid,'/version',v,'/YRI_',pops[pid],'/out_',mod,'.BestSimsParamStats_Obs0.txt',sep=''),header=T)
      lines(density(post[,post_cols[i]]),col=mc[pid])
    }
  }
}
dev.off()

########################
###### Plot PLS clouds and real obs
########################
## We read the entire data into one table (sims,obs)
## NOTE: only selected columsn (assuming 5 PLS!)
mycols <- rep('NULL', 10); mycols[c(6:10)] <- NA
## for transparency (alpha flag in plot() )
library(scales)

if(pps=='all'){pps <- c('eur','asi')}
for (c in pps){
  pdf(paste('/home/felix_schulze/projects/trpm8_proj/pdf/abc/run',rid,'/version',v,'_PLSclouds_',c,'.pdf',sep=''),height=10,width=10)
  d <- data.frame()
  ntr <- read.table(paste('/mnt/scratch/felix/trpm8/abc/run',rid,'_out/transf/all_sumstat_5ParamBT_trans_',c,'_ntr_run',rid,'v',v,'_5PLS.tsv',sep=""),header=T,colClasses=mycols,nrows=10000)
  ssv <- read.table(paste('/mnt/scratch/felix/trpm8/abc/run',rid,'_out/transf/all_sumstat_5ParamBT_trans_',c,'_ssv_run',rid,'v',v,'_5PLS.tsv',sep=""),header=T,colClasses=mycols,nrows=10000)
  sdn <- read.table(paste('/mnt/scratch/felix/trpm8/abc/run',rid,'_out/transf/all_sumstat_5ParamBT_trans_',c,'_sdn_run',rid,'v',v,'_5PLS.tsv',sep=""),header=T,colClasses=mycols,nrows=10000)
  d <- rbind(ntr,ssv,sdn) # order persists
  for (p in pops){
    obs <- read.table(paste('/mnt/scratch/felix/trpm8/abc/obs/transf/run',rid,'/obs_trans_abcRun',rid,'v',v,'scp_YRI_',p,'_5PLS.tsv',sep=""),header=T)
    d <- rbind(d,obs)
  }
  myColors <- c( rep('grey',nrow(ntr)) , rep('darkgoldenrod1',nrow(ssv)) , rep('lightblue',nrow(sdn)) , mc )
  myPCH <- c( rep('.',nrow(ntr)) , rep('.',nrow(ssv)) , rep('.',nrow(sdn)) , rep('*',length(pops)) )
  weighting <- c( rep(1, (nrow(ntr)+nrow(sdn)+nrow(ssv)) ), rep(4,length(pops)) )
  plot(d,pch=myPCH,col=alpha(myColors,0.5),cex=weighting)
  dev.off()
}




