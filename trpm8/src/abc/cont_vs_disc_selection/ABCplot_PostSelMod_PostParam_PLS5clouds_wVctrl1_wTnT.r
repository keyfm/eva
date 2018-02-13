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
pdf(paste('/home/felix_schulze/projects/trpm8_proj/pdf/abc/run',rid,'/version',v,'_modelPosterior_',pps,'_wTnT.pdf',sep=''))
md <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run',rid,'/version',v,'/marginalDensities_ntr_sdn_ssv_nTwT.txt',sep=""),header=F)
d <- t(md[,2:7])
colnames(d) <- md[,1]
for (i in 1:ncol(d)){
  d[,i] <- d[,i]/sum(d[,i])
}
barplot(d, col=c('white','grey','lightblue','lightblue4','darkgoldenrod1','darkgoldenrod4'),main=paste('run',rid,';ntr=white;sdn=blue,ssv=red;nTrunc=light,wTrunc=dark',sep=''),las=3)
dev.off()

########################
###### Plot Posterior Prob for Model Choice
########################
pdf(paste('/home/felix_schulze/projects/trpm8_proj/pdf/abc/run',rid,'/version',v,'_paramPosterior_',pps,'_wTnT.pdf',sep=''),height=20,width=15)
par(mfrow=c(4,3))
cont = c('eur') ## note for now cont only used for prior plotting...ASI and EUR no difference! -> Thats why hard coded!
for (mod in c('ssv','sdn')){
 if(mod=='sdn'){prior_cols <- c(1,2,3);post_cols <- c(3,4,5)} else {prior_cols <- c(1,3,5);post_cols <- c(3,4,5)}
 for (trunc in c('nTrunc','wTrunc')){
  prior <- read.table(paste('/mnt/scratch/felix/trpm8/abc/run',rid,'_out/transf/all_sumstat_5ParamBT_trans_',cont,'_',mod,'_',trunc,'_run',rid,'v',v,'_5PLS.tsv',sep=""),header=T,nrows = 100000)
  for (i in 1:3){
    l=density(prior[,prior_cols[i]])
    mymax <- max(l$y)+((max(l$y))*3)
    hist(prior[,prior_cols[i]],freq=FALSE,ylim=c(0,mymax),main=paste('run',rid,'v',v,'; ',mod,trunc,' ',colnames(prior)[ prior_cols[i] ],sep=''),xlab=colnames(prior)[i])
    for (pid in 1:length(pops)){
      post <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run',rid,'/version',v,'/YRI_',pops[pid],'/out_',mod,'_',trunc,'.BestSimsParamStats_Obs0.txt',sep=''),header=T)
      lines(density(post[,post_cols[i]]),col=mc[pid])
    }
  }
}}
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
  pdf(paste('/home/felix_schulze/projects/trpm8_proj/pdf/abc/run',rid,'/version',v,'_PLSclouds_',c,'_wTnT.pdf',sep=''),height=10,width=10)
  d <- data.frame()
  ntr_nT <- read.table(paste('/mnt/scratch/felix/trpm8/abc/run',rid,'_out/transf/all_sumstat_5ParamBT_trans_',c,'_ntr_nTrunc_run',rid,'v',v,'_5PLS.tsv',sep=""),header=T,colClasses=mycols,nrows=10000)
  ssv_nT <- read.table(paste('/mnt/scratch/felix/trpm8/abc/run',rid,'_out/transf/all_sumstat_5ParamBT_trans_',c,'_ssv_nTrunc_run',rid,'v',v,'_5PLS.tsv',sep=""),header=T,colClasses=mycols,nrows=10000)
  sdn_nT <- read.table(paste('/mnt/scratch/felix/trpm8/abc/run',rid,'_out/transf/all_sumstat_5ParamBT_trans_',c,'_sdn_nTrunc_run',rid,'v',v,'_5PLS.tsv',sep=""),header=T,colClasses=mycols,nrows=10000)
  ntr_wT <- read.table(paste('/mnt/scratch/felix/trpm8/abc/run',rid,'_out/transf/all_sumstat_5ParamBT_trans_',c,'_ntr_wTrunc_run',rid,'v',v,'_5PLS.tsv',sep=""),header=T,colClasses=mycols,nrows=10000)
  ssv_wT <- read.table(paste('/mnt/scratch/felix/trpm8/abc/run',rid,'_out/transf/all_sumstat_5ParamBT_trans_',c,'_ssv_wTrunc_run',rid,'v',v,'_5PLS.tsv',sep=""),header=T,colClasses=mycols,nrows=10000)
  sdn_wT <- read.table(paste('/mnt/scratch/felix/trpm8/abc/run',rid,'_out/transf/all_sumstat_5ParamBT_trans_',c,'_sdn_wTrunc_run',rid,'v',v,'_5PLS.tsv',sep=""),header=T,colClasses=mycols,nrows=10000)
  d <- rbind(ntr_nT,ntr_wT,ssv_nT,ssv_wT,sdn_nT,sdn_wT) # order persists
  for (p in pops){
    obs <- read.table(paste('/mnt/scratch/felix/trpm8/abc/obs/transf/run',rid,'/obs_trans_abcRun',rid,'v',v,'scp_YRI_',p,'_5PLS.tsv',sep=""),header=T)
    d <- rbind(d,obs)
  }
  myColors <- c( rep('grey',nrow(ntr_nT)), rep('darkgrey',nrow(ntr_wT)) , rep('darkgoldenrod1',nrow(ssv_nT)), rep('darkgoldenrod4',nrow(ssv_wT)) , rep('lightblue',nrow(sdn_nT)), rep('lightblue4',nrow(sdn_wT)) , mc )
  myPCH <- c( rep('.',nrow(ntr_nT)+nrow(ntr_wT)+nrow(sdn_nT)+nrow(sdn_wT)+nrow(ssv_nT)+nrow(ssv_wT)) , rep('*',length(pops)) )
  weighting <- c( rep(1, (nrow(ntr_nT)+nrow(ntr_wT)+nrow(sdn_nT)+nrow(sdn_wT)+nrow(ssv_nT)+nrow(ssv_wT)) ), rep(4,length(pops)) )
  plot(d,pch=myPCH,col=alpha(myColors,0.5),cex=weighting)
  dev.off()
}




