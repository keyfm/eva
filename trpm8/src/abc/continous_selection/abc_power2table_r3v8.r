res <- list()
for (c in c('eur','asi')){
    data <- list()
    data[['ntr']] <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run3/power/',c,'_run3v8_ntr_cap10000/marginal_densities_ntr_sdn_ssv.tsv',sep=""));colnames(data[['ntr']]) <- c('ntr','sdn','ssv')
    data[['sdn']] <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run3/power/',c,'_run3v8_sdn_cap10000/marginal_densities_ntr_sdn_ssv.tsv',sep=""));colnames(data[['sdn']]) <- c('ntr','sdn','ssv')
    data[['ssv']] <- read.table(paste('/mnt/scratch/felix/trpm8/abc/analysis/run3/power/',c,'_run3v8_ssv_cap10000/marginal_densities_ntr_sdn_ssv.tsv',sep=""));colnames(data[['ssv']]) <- c('ntr','sdn','ssv')
    res[[c]] <- matrix(ncol=3,nrow=3);colnames(res[[c]]) <- c('TP','FP','FN');rownames(res[[c]]) <- c('ntr','sdn','ssv')
    all.m <- c('ntr','sdn','ssv')
    for (m in c('ntr','sdn','ssv')){
        other.m <- all.m[ -which(all.m==m) ]
        ## TP
        res[[c]][m,'TP'] <- sum(data[[ m ]][,m] > data[[ m ]][,other.m[1]] & data[[ m ]][,m] > data[[ m ]][,other.m[2]]) / nrow(data[[m]])
        ## FP
        vec.fp <- vector(length=length(other.m));names(vec.fp) <- other.m
        vec.n <- vector(length=length(other.m));names(vec.fp) <- other.m
        for (o.m in other.m){
          a <- table(apply(data[[o.m]],1,which.max));names(a) <- c('ntr','sdn','ssv')
          vec.fp[o.m] <- a[m]
          vec.n[o.m] <- sum(a)
        }
        res[[c]][m,'FP'] <- sum(vec.fp)/sum(vec.n)
        ## FN
        res[[c]][m,'FN'] <- 1-res[[c]][m,'TP']
    }
    ## write output
    write.table(res[[c]],file=paste('/home/felix_schulze/projects/trpm8_proj/table/abc/abc_TP_FP_FN_cap10000_run3_v8_',c,'.tsv',sep=''),quote=F,sep="\t")
}
