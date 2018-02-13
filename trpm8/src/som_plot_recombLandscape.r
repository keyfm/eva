############## 1000genomes Recombination rates based on illumina omni chips (~2M SNPs)
############## ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/
# extracted recombination rate for locus of interest and transformed >> /home/felix_schulze/projects/trpm8_proj/recrate/ABCtrg_3x65kb_RecRt1000gPops_noAME_cMpBP.bed

##plot mean all nonAME
d <- read.table('/home/felix_schulze/projects/trpm8_proj/recrate/ABCtrg_3x65kb_RecRt1000gPops_noAME_cMpBP.bed',header=T,sep="\t")
pdf("/home/felix_schulze/projects/trpm8_proj/pdf/rec_rate/1kg_omni_recRt_MEANnonAMEpops_ABCtrg_3x65kb_2rdmPOPpCont.pdf")
plot(d[,2],rowMeans(d[,c('YRI','LWK','GBR','TSI','CHB','GIH')])*1000000,typ='s',ylim=c(0,40),ylab="Recombination Rate (Cm/Mb)",xlab="Position on chromosome 2")
abline(v=c(234826043,234928166),lty=2,col="blue") # TRPM8 refSeq
abline(v=c(234810004,234875139),lty=2,col="red") # TRPM8 target
abline(h=5,lty=2,col="black") # RecRate cut off for peak consideration or background
dev.off()
