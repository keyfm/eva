##############################
###### Extract new mean temperatures
##############################
### Read in Grid File
cru <- read.table(gzfile('/home/felix_schulze/projects/trpm8_proj/temp/cru_ts3.23.2001.2010.tmp.dat.gz'),header=F)
## header is longitude from 179.75W to 179.75E (gives midpoint of 0.5degree step)
colnames(cru) <- paste( abs(seq(-179.75, 179.75,0.5)) , c(rep("W",360),rep("E",360)), sep="")
## rownames are similar: from 89.75S to 89.75N (360rows also in steps of 0.5
## the first 360 rows show the data for Jan 1901, next 360 rows the data for Feb 1901, next 360 rows for March 1901
one.month.set <- paste( abs(seq(-89.75,89.75,0.5)), c(rep("S",180),rep("N",180)), sep="")
one.year.set <- as.vector(sapply(1:12,function(x) paste(x,one.month.set,sep="_")))
one.decade.set <- as.vector(sapply(1:10,function(x) paste(x,one.year.set,sep="_")))
row.names(cru) <- one.decade.set

### Extract temperatures per month/ per year
pop.data <- read.table('/home/felix_schulze/projects/trpm8_proj/temp/temp_1000g_annual_under15_wLAT_wLON.tsv',stringsAsFactors=F,sep="\t",header=T)
tmp <- list()
for (i in 1:nrow(pop.data)){
  lat <- pop.data[i,'LatGrid']
  lon <- pop.data[i,'LonGrid']
  tmp[[ pop.data[i,'pop'] ]] <- matrix(nrow=10,ncol=12)
  for (m in 1:12){
    for (y in 1:10){
      tmp[[ pop.data[i,'pop'] ]][y,m] <- cru[ paste(y,m,lat,sep="_") , lon ]/10 # remember tmp values are *10 in CRU data
    }
  }
}
## worked. temps look good and no -999!

## get mean temp and ratio months below 15
meanT <- vector()
for (p in names(tmp)){ meanT[p] <- mean(tmp[[p]])}

