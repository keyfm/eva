## load data made in trpm8_fst_ihs_xp.r l.231 
load('/home/felix_schulze/projects/trpm8_proj/fst/allPop_vs_YRI_trpm8_refSeq_pm20kb_Fst_FstEmpP.RData')
load("/home/felix_schulze/projects/trpm8_proj/xp_res/allPop_vs_YRI_trpm8_refSeq_pm20kb_XPnormP.RData")
meanFstXP <- read.table('/home/felix_schulze/projects/trpm8_proj/table/fst_xp_meanEmpP_chr2genes.txt',header=T,row.names=1) # /Users/fmk/Documents/eva/trpm8/readme/trpm8_fst_ihs_xp.r l.555 for rev#1
ymax <- max(unlist(lapply(fst.p, "[", c(4))))
Varcex <- 1
TrgCex <- 1.5
AxisCex <- 1.65
varlwd <- 3
Varcex <- 1.1
pdf('/home/felix_schulze/projects/trpm8_proj/pdf/fst_xpehh/trpm8_refseq_pm20kb_allpops_fst_xp.pdf',width=30,height=18)
par(mfrow=c(4,5))
pops=c('CEU','GBR','TSI','FIN','IBS','CDX','CHB','CHS','KHV','JPT','BEB','GIH','ITU','PJL','STU','LWK','ESN','GWD','MSL')
for (p in pops){
  ## plotting empty plot
  par(xpd=F)
  plot("",xlim=c(234806042,234948166),ylim=c(0,ymax),las=1,cex.axis=1.2,xaxt='n',xlab="Position on chromosome 2",ylab="-log(Emprical P-value)",yaxt='n',cex.lab=AxisCex,mgp=c(2,1,.5))
  lim=par("usr")
  ## rect(234928166,lim[3],234826042,lim[4],col=colors()[353],border=colors()[353])
  abline(v=c(234826042,234928166),col='grey',lty=3,lwd=3) # trpm8 refseq
  abline(v=c(234810003,234875139),col='red',lty=3,lwd=3) # target roi
  # abline(h=-log10(0.05),col='grey',lty=5,lwd=3)
  abline(h=c(meanFstXP[p,'Fst_mean'],meanFstXP[p,'XP_mean']),col=c('cornflowerblue','grey'),lty=5,lwd=3)
  axis(side=1,at=c(234928166,234826042),labels=c(234928166,234826042),cex.axis=AxisCex)
  axis(side=2,at=seq(0,5,1),labels=seq(0,5,1),cex.axis=AxisCex,las=1) # yaxis=0,1,2,3
  rect(234928166,ymax+0.2, 234826042, ymax+0.4,col="burlywood", border="black", xpd=T,cex=AxisCex)
  text(mean(234928166:234826042),ymax+0.4,"TRPM8",font=4,pos=3,xpd=T,cex=AxisCex)
  text(mean(234806042:234826042),ymax+0.4,font=2,p,pos=3,xpd=T,cex=AxisCex+0.5)
  ## plot Fst empPvalues
  points(fst.p[[p]][,"p"],fst.p[[p]][,"fst_empP_log"] ,pch=19,col="lightblue",cex=Varcex)
  migrIDX <- which( fst.p[[p]][,"p"] == 234825093 )
  points(fst.p[[p]][migrIDX,"p"],fst.p[[p]][migrIDX,"fst_empP_log"] ,pch=19,col="red",cex=TrgCex)
  abline(h=c(meanFstXP[p,'Fst_mean'],meanFstXP[p,'XP_mean']),col=c('cornflowerblue','grey'),lty=5,lwd=3)
  ## plot XPEHH empPval smoothed
  xdata <- xp_normP[[p]]
  roi=which(xdata[,"chr"] == 2 & xdata[,"pos"] >= 234806042 & xdata[,"pos"] <= 234948166)
  smoothingSpline = smooth.spline(xdata[roi,"pos"],-log10(xdata[roi,"normxp_pvalue"]),spar=0.01)
  lines(smoothingSpline,col="grey80",lwd=varlwd)
  points(xdata[roi,"pos"], -log10(xdata[roi,"normxp_pvalue"]),pch=23,bg="grey80",cex=Varcex)
  rsID <- which(xdata[,'pos'] == 234825093 )
  points(xdata[rsID,"pos"], -log10(xdata[rsID,"normxp_pvalue"]),pch=23,bg="red",cex=Varcex)
}
dev.off()

