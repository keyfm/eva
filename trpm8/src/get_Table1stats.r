##################
######### TABLE 1. Latitude, Temp, DAF, and Stats for allPOPs (excl. AME, ASW)
######### Stats for Migraine SNP only: 2_234825093       rs10166942 
##################
daf <- read.table('/home/felix_schulze/projects/trpm8_proj/daf/daf_roi_2_234810004_234875139_sglPops_target_allSNPs.bed',header=T)
daf <- daf[ which(daf[,'dbSNP']=='rs10166942') , ]
temp.lat <- read.table('~/projects/trpm8_proj/temp/temp_1000g_annual_under15_wLAT_w1kgAC_wCEU.tsv',header=T) # from trpm8_glm.R l.114
temp.lat <- unique(temp.lat[,1:4])
rownames(temp.lat) <- temp.lat[,1]
## fst from Muslihs data ~/from_muslih/TRPM8/fst_res/Code/TRPM8\ readme_extended40kbp.txt l.104
## NOTE: The Fst here that went into the table at line 60 are incorrecect (JS changed his calculations after Mulsihs departure). Correct Fsts generated at line 63+
load("/home/felix_schulze/from_muslih/TRPM8/fst_res/fst_list_all26_poppairs_candTRPM8_extended40kbp.RData")
idx <- which( out$positions[,3] == 'rs10166942')
fst <- out[ grep('YRI',names(out)) ]
myfst <- vector(length=19)
names(myfst) <- c("ESN","GWD","LWK","MSL","CEU","GBR","TSI","FIN","IBS","CDX","CHB","CHS","KHV","JPT","BEB","GIH","ITU","PJL","STU")
for (p in c("ESN","GWD","LWK","MSL","CEU","GBR","TSI","FIN","IBS","CDX","CHB","CHS","KHV","JPT","BEB","GIH","ITU","PJL","STU")){
  myfst[p] <- round(as.numeric(unlist(fst[ grep(p,names(fst)) ])[idx]),4)
}

## XPEHH from Mulsih: ~/from_muslih/TRPM8/fst_res/Code/XP_EHH.txt l.84++
load('/home/felix_schulze/from_muslih/TRPM8/fst_res/res_list_XP_extended_region.RData',header=T)
idx <- which( out$positions[,3] == 'rs10166942')
myxp <- vector(length=19)
names(myxp) <- c("ESN","GWD","LWK","MSL","CEU","GBR","TSI","FIN","IBS","CDX","CHB","CHS","KHV","JPT","BEB","GIH","ITU","PJL","STU")
for (p in c("ESN","GWD","LWK","MSL","CEU","GBR","TSI","FIN","IBS","CDX","CHB","CHS","KHV","JPT","BEB","GIH","ITU","PJL","STU")){
  if(length(grep(p,names(out)))!=0){
    myxp[p] <- round(as.numeric(unlist(out[ grep(p,names(out)) ])[idx]),3)
  } else {
    myxp[p] <- NA
  }
}

## iHS from Muslih: ~/from_muslih/TRPM8/fst_res/Code/iHS.txt l.9++
load('/home/felix_schulze/from_muslih/TRPM8/fst_res/ihs_res.RData')
ihs <- out
idx <- which( ihs$positions[,3] == 'rs10166942')
myihs <- vector(length=20)
names(myihs) <- c("YRI","ESN","GWD","LWK","MSL","CEU","GBR","TSI","FIN","IBS","CDX","CHB","CHS","KHV","JPT","BEB","GIH","ITU","PJL","STU")
for (p in c("YRI","ESN","GWD","LWK","MSL","CEU","GBR","TSI","FIN","IBS","CDX","CHB","CHS","KHV","JPT","BEB","GIH","ITU","PJL","STU")){
  myihs[p] <- round(as.numeric(unlist(ihs[ grep(p,names(ihs)) ])[idx]),3)
}
## write table
data <- matrix(ncol=8,nrow=20)
rownames(data) <- c("YRI","ESN","GWD","LWK","MSL","CEU","GBR","TSI","FIN","IBS","CDX","CHB","CHS","KHV","JPT","BEB","GIH","ITU","PJL","STU")
colnames(data) <- c("Pop","Latitude", "Temp","under15", "DAF", "Fst", "XP-EHH", "iHS")
for (p in rownames(data)){
  data[p,] <- cbind( p, temp.lat[p,'latitude'] , temp.lat[p,'annual'],temp.lat[p,'under15'] , round(daf[1,p],2) , myfst[p], myxp[p] ,myihs[p] )
}
data <- data[order(as.numeric(data[,'Latitude']),decreasing=T) , ]
write.table(data, file="/home/felix_schulze/projects/trpm8_proj/table/Lat_Temp_DAF_under15_fst_xp_ihs.txt",quote=F,sep="\t",row.names=F,col.names=T)

#########
### Fst vs YRI
#########
## files generated in trpm8_fst_ihs_xp.r l.230; here I extract the empPval and actual Fst vs YRI
load("/home/felix_schulze/projects/trpm8_proj/fst/allPop_vs_YRI_trpm8_refSeq_pm20kb_Fst_FstEmpP.RData")
## Fst target allel
rsd <- matrix(ncol=2,nrow=19)
rownames(rsd) <- c("ESN","GWD","LWK","MSL","CEU","GBR","TSI","FIN","IBS","CDX","CHB","CHS","KHV","JPT","BEB","GIH","ITU","PJL","STU")
colnames(rsd) <- c('fst.val','fst.p')
for (p in rownames(rsd)){
  d <- fst.val[[p]]
  idx <- which(d[,2]==234825093)
  rsd[p,'fst.val'] <- d[idx,'fst']
  d <- fst.p[[p]]
  idx <- which(d[,2]==234825093)
  rsd[p,'fst.p'] <- d[idx,'fst_empP']
}
 write.table(rsd,'~/projects/trpm8_proj/fst/raw_and_empP_fst_rs10166942_vs_YRI.tsv',quote=F,sep='\t',row.names=T,col.names=T)
## ratio of significant Fst values below 0.05 in target region (/home/felix_schulze/projects/trpm8_proj/bed/trpm8_target_roi.bed: 2	234810003	234875139)
load("/home/felix_schulze/projects/trpm8_proj/fst/allPop_vs_YRI_trpm8_refSeq_pm20kb_Fst_FstEmpP.RData")
rsd <- matrix(ncol=2,nrow=19)
rownames(rsd) <- c("ESN","GWD","LWK","MSL","CEU","GBR","TSI","FIN","IBS","CDX","CHB","CHS","KHV","JPT","BEB","GIH","ITU","PJL","STU")
colnames(rsd) <- c('fst.ratio005','fst.ratio001')
for (p in rownames(rsd)){
  d <- fst.p[[p]]
  idx <- which(d[,2] >= 234810003 & d[,2] <= 234875139)
  myfst <- d[idx,'fst_empP']
  rsd[p,'fst.ratio005'] <- sum(myfst < 0.05)/length(myfst)
  rsd[p,'fst.ratio001'] <- sum(myfst < 0.01)/length(myfst)
}
write.table(rsd,'~/projects/trpm8_proj/fst/ratio_below005_001_empP_fst_65trgROI_vs_YRI.tsv',quote=F,sep='\t',row.names=T,col.names=T)

