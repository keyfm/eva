###################################
####### Plot boxplots of distribution of pairwise counts for each haplotype
###################################
# input data generated w/ pairwiseDiff.py
library(ggplot2)

myymax=80
bx=1.5
ytxt=5
## AFR
for (p in c('YRI','LWK','ESN','GWD','MSL','ACB','ASW')){
d=read.table(paste('/tmp/fmk_',p,'_pw.txt',sep=""),stringsAsFactors=F)
ids=sort(unique(unlist(strsplit(d[,1],split="-",fixed=T))))
pdf(paste("/home/felix_schulze/projects/trpm8_proj/pdf/pairwiseDiff/boxplots/pw_target_allSNPs_condMigraine_",p,"_bp.pdf",sep=""),width=length(ids)/1.5,height=8)
boxplot(rnorm(200), xlim = c(1,length(ids)+1),ylim=c(0,myymax), ylab='Pairwise Differences', boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1),cex.axis=2) # empty plot
text(1:length(ids), par("usr")[3]-ytxt, labels = ids, srt = 45, pos = 1, xpd = TRUE,cex=1.5)
for (i in 1:length(ids)){
 idx=grep(ids[i],d[,1],fixed=T)
 # set color
 if(median(d[idx,2])> 10){
 	mycol='orange'
 } else {
 	mycol='grey'
 }
 boxplot(d[idx,2], xaxt = "n", yaxt = "n", add=T,at=i+0.5,boxwex=bx,col=mycol)
 }
dev.off()
}


## EUR et al.
for (p in c('CEU','GBR','TSI','FIN','IBS')){
d=read.table(paste('/tmp/fmk_',p,'_pw.txt',sep=""),stringsAsFactors=F)
ids=sort(unique(unlist(strsplit(d[,1],split="-",fixed=T))))
pdf(paste("/home/felix_schulze/projects/trpm8_proj/pdf/pairwiseDiff/boxplots/pw_target_allSNPs_condMigraine_",p,"_bp.pdf",sep=""),width=length(ids)/1.5,height=8)
boxplot(rnorm(200), xlim = c(1,length(ids)+1),ylim=c(0,myymax), ylab='Pairwise Differences', boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1),cex.axis=2) # empty plot
text(1:length(ids), par("usr")[3]-ytxt, labels = ids, srt = 45, pos = 1, xpd = TRUE,cex=1.5)
for (i in 1:length(ids)){
 idx=grep(ids[i],d[,1],fixed=T)
 # set color
 if(median(d[idx,2])> 10){
 	mycol='orange'
 } else {
 	mycol='grey'
 }
 boxplot(d[idx,2], xaxt = "n", yaxt = "n", add=T,at=i+0.5,boxwex=bx,col=mycol)
 }
dev.off()
}

for (p in c('CDX','CHB','CHS','KHV','JPT')){
d=read.table(paste('/tmp/fmk_',p,'_pw.txt',sep=""),stringsAsFactors=F)
ids=sort(unique(unlist(strsplit(d[,1],split="-",fixed=T))))
pdf(paste("/home/felix_schulze/projects/trpm8_proj/pdf/pairwiseDiff/boxplots/pw_target_allSNPs_condMigraine_",p,"_bp.pdf",sep=""),width=length(ids)/1.5,height=8)
boxplot(rnorm(200), xlim = c(1,length(ids)+1),ylim=c(0,myymax), ylab='Pairwise Differences', boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1),cex.axis=2) # empty plot
text(1:length(ids), par("usr")[3]-ytxt, labels = ids, srt = 45, pos = 1, xpd = TRUE,cex=1.5)
for (i in 1:length(ids)){
 idx=grep(ids[i],d[,1],fixed=T)
 # set color
 if(median(d[idx,2])> 10){
 	mycol='orange'
 } else {
 	mycol='grey'
 }
 boxplot(d[idx,2], xaxt = "n", yaxt = "n", add=T,at=i+0.5,boxwex=bx,col=mycol)
 }
dev.off()
}

for (p in c('BEB','GIH','ITU','PJL','STU')){
d=read.table(paste('/tmp/fmk_',p,'_pw.txt',sep=""),stringsAsFactors=F)
ids=sort(unique(unlist(strsplit(d[,1],split="-",fixed=T))))
pdf(paste("/home/felix_schulze/projects/trpm8_proj/pdf/pairwiseDiff/boxplots/pw_target_allSNPs_condMigraine_",p,"_bp.pdf",sep=""),width=length(ids)/1.5,height=8)
boxplot(rnorm(200), xlim = c(1,length(ids)+1),ylim=c(0,myymax), ylab='Pairwise Differences', boxfill=rgb(1, 1, 1, alpha=1), border=rgb(1, 1, 1, alpha=1),cex.axis=2) # empty plot
text(1:length(ids), par("usr")[3]-ytxt, labels = ids, srt = 45, pos = 1, xpd = TRUE,cex=1.5)
for (i in 1:length(ids)){
 idx=grep(ids[i],d[,1],fixed=T)
 # set color
 if(median(d[idx,2])> 10){
 	mycol='orange'
 } else {
 	mycol='grey'
 }
 boxplot(d[idx,2], xaxt = "n", yaxt = "n", add=T,at=i+0.5,boxwex=bx,col=mycol)
 }
dev.off()
}

###################################
####### Plot proportion of SNP exclusively present on median pi > 10 haplotypes, being also present on ancestral migraine linked haplotype
###################################
# input table gerneated in bash code
pdf('/Users/fmk/Dropbox/work/trpm8/pdf/pairwise_ana/proportion_highPWder_onAncHap_perPop.pdf')
mc <- c(colorRampPalette(c("orange", "yellow"))(5),colorRampPalette(c("lightblue", "darkblue"))(5),colorRampPalette(c("lightgreen", "darkgreen"))(5),colorRampPalette(c("pink", "purple"))(5))
names(mc) <- c('YRI','ESN','MSL','GWD','LWK','CEU','GBR','TSI','FIN','IBS','CDX','CHB','CHS','KHV','JPT','BEB','GIH','ITU','PJL','STU')
d <- read.table('/Users/fmk/Dropbox/work/trpm8/tmp_from_server/highPWder_on_ancestral_perPop.txt',row.names=1)
d <- d[d[,1] != 0 ,] # remove the once that do not have any Median-PW > 10
plot("",xlim=c(1,nrow(d)),ylim=c(0,1),ylab="Proportion excess pairwise diversity alleles on ancestral haplotype",xlab='',xaxt='n')
ctr=1
for (p in c('YRI','ESN','MSL','GWD','LWK','CEU','GBR','TSI','FIN','IBS','CDX','CHB','CHS','KHV','JPT','BEB','GIH','ITU','PJL','STU')){
    if( p %in% rownames(d)){
        points(x=ctr,y=(d[p,2]/d[p,1]),pch=21 , cex=2.5 , bg=mc[p],col='black')
        ctr <- ctr+1
    }
}
legend('bottom',legend=paste(rownames(d)," (",d[,1],")",sep=""),fill=mc[rownames(d)],ncol=4,bty='n')    
dev.off()

