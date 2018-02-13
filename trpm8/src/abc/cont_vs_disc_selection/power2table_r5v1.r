## run5v1 ABC (truncated vs. continous selection)
res <- list()
for (c in c('eur','asi')){
    data <- list()
    data[['NTR']] <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run5/power/',c,'_run5v1_ntr_nTrunc_cap10000/marginal_densities_ntr_sdn_ssv_nTwT.tsv',sep=""));colnames(data[['NTR']]) <- c('NTR','NTR2','SDN_continous','SDN_stop','SSV_continous','SSV_stop')
    data[['SDN_continous']] <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run5/power/',c,'_run5v1_sdn_nTrunc_cap10000/marginal_densities_ntr_sdn_ssv_nTwT.tsv',sep=""));colnames(data[['SDN_continous']]) <- c('NTR','NTR2','SDN_continous','SDN_stop','SSV_continous','SSV_stop')
    data[['SSV_continous']] <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run5/power/',c,'_run5v1_ssv_nTrunc_cap10000/marginal_densities_ntr_sdn_ssv_nTwT.tsv',sep=""));colnames(data[['SSV_continous']]) <- c('NTR','NTR2','SDN_continous','SDN_stop','SSV_continous','SSV_stop')
    data[['NTR2']] <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run5/power/',c,'_run5v1_ntr_wTrunc_cap10000/marginal_densities_ntr_sdn_ssv_nTwT.tsv',sep=""));colnames(data[['NTR2']]) <- c('NTR','NTR2','SDN_continous','SDN_stop','SSV_continous','SSV_stop')
    data[['SDN_stop']] <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run5/power/',c,'_run5v1_sdn_wTrunc_cap10000/marginal_densities_ntr_sdn_ssv_nTwT.tsv',sep=""));colnames(data[['SDN_stop']]) <- c('NTR','NTR2','SDN_continous','SDN_stop','SSV_continous','SSV_stop')
    data[['SSV_stop']] <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run5/power/',c,'_run5v1_ssv_wTrunc_cap10000/marginal_densities_ntr_sdn_ssv_nTwT.tsv',sep=""));colnames(data[['SSV_stop']]) <- c('NTR','NTR2','SDN_continous','SDN_stop','SSV_continous','SSV_stop')
    res[[c]] <- matrix(ncol=3,nrow=6);colnames(res[[c]]) <- c('TP','FP','FN');rownames(res[[c]]) <- c('NTR','NTR2','SDN_continous','SDN_stop','SSV_continous','SSV_stop')
    all.m <- c('NTR','NTR2','SDN_continous','SDN_stop','SSV_continous','SSV_stop')
    for (m in c('NTR','NTR2','SDN_continous','SDN_stop','SSV_continous','SSV_stop')){
        other.m <- all.m[ -which(all.m==m) ]
        ## TP
        res[[c]][m,'TP'] <- sum(data[[ m ]][,m] > data[[ m ]][,other.m[1]] & data[[ m ]][,m] > data[[ m ]][,other.m[2]] & data[[ m ]][,m] > data[[ m ]][,other.m[3]] & data[[ m ]][,m] > data[[ m ]][,other.m[4]] & data[[ m ]][,m] > data[[ m ]][,other.m[5]]) / nrow(data[[m]])
        ## FP
        vec.fp <- vector(length=length(other.m));names(vec.fp) <- other.m
        vec.n <- vector(length=length(other.m));names(vec.fp) <- other.m
        for (o.m in other.m){
          a <- table(apply(data[[o.m]],1,which.max));names(a) <- c('NTR','NTR2','SDN_continous','SDN_stop','SSV_continous','SSV_stop')
          vec.fp[o.m] <- a[m]
          vec.n[o.m] <- sum(a)
        }
        res[[c]][m,'FP'] <- sum(vec.fp)/sum(vec.n)
        ## FN
        res[[c]][m,'FN'] <- 1-res[[c]][m,'TP']
    }
    ## write output
    write.table(res[[c]],file=paste('/home/felix_schulze/projects/trpm8_proj/table/abc/abc_TP_FP_FN_cap10000_run5_v1_',c,'.tsv',sep=''),quote=F,sep="\t")
}

