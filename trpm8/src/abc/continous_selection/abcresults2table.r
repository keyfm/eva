## bayes factor and SSV posterior prob
rid <- 3
v <- 8
md <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run',rid,'/version',v,'/marginalDensities_ntr_sdn_ssv.txt',sep=""),stringsAsFactors=F)
stime <- matrix(ncol=3,nrow=length(md[,1]));rownames(stime) <- md[,1]
sstr <- matrix(ncol=3,nrow=length(md[,1]));rownames(sstr) <- md[,1]
sfrq <- matrix(ncol=3,nrow=length(md[,1]));rownames(sfrq) <- md[,1]
for (p in md[,1]){
    pp <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run',rid,'/version',v,'/YRI_',p,'/out_ssv.BestSimsParamStats_Obs0.txt',sep=''),header=T)

    stime[p,1] <- round(median(pp[,'seltime_bt']),0)
    stime[p,2:3] <- round(quantile(pp[,'seltime_bt'], probs = c(0.025,0.975)),0)
    sstr[p,1] <- median(pp[,'selstrP2_bt'])
    sstr[p,2:3] <- quantile(pp[,'selstrP2_bt'], probs = c(0.025,0.975))
    sfrq[p,1] <- median(pp[,'freqssP2'])
    sfrq[p,2:3] <- quantile(pp[,'freqssP2'], probs = c(0.025,0.975))
    
}
out.mat <- cbind(md[,1], round((md[,4]/(md[,2]+md[,3])),1),round((md[,4]/(md[,2]+md[,3]+md[,4])),4),stime,sstr,sfrq )
colnames(out.mat) <- c('Pop','Bayes_Factor','Posterior_P_SSV','SelTime_Med','2.5%','97.5%','SelStr_pop2_Med','2.5%','97.5%','FreqencySelection_Med','2.5%','97.5%')
write.table(out.mat, file=paste("/home/felix_schulze/projects/trpm8_proj/table/bayFac_Pssv_STimeMed_SStrMed_SFrwMed_run",rid,"_",v,".txt",sep=""),quote=F,row.names=F)
