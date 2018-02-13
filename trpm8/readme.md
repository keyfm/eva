# README
### by Felix M Key

Contact the author at key [at] shh [dot] mpg [dot] de for any questions about this code
or correspondence about this work.

This README guides you through the process of analyses and figure/table creation presented in 
Key et al. Human local adaptation of the TRPM8 cold receptor along a latitudinal cline 
Please cite this paper if any meaningful piece of it was recycled for your analysis.

The scripts likely need to be modified to suit the file structure and source code paths on any other local computer.

The package is structured, by all that resembles data, can be found in the data folder; while all scripts and in-house dependencies are in the src folder. ABC is further sub-divided in its respective single analyses.
Please, note that some figures were made by hand for representative purposes of for improved representation. For those the underlying data is included (if any).
Enjoy!


## Latitude, Temperature and Signatures of Natural Selection

VCF-coded genotype info for all populations in our TRPM8 target region can be found in data/
Latitude and Termperature raw data was obtained as described in the manuscript
extract mean Temperature for geographic locations from CRU matrix

```
Rscript src/get_meanT_CRUdata.r
```
## Natural Selection Signatures and Genotype Population info (Table 1)
summary statisitics calculated as described in manuscript
all data for Table 1 can be generated with following R script

```
Rscript src/get_Table1stats.r
```


## PGLS analysis

NOTE: add plotting!
documented R code to run PGLS analysis

```
Rscript src/calcPGLS.r
```

## GLMM analysis
NOTE: add plotting!
1000Genomes based analysis 
documented R code to run  GLMM analysis

```
Rscript src/calcGLMM_1kg.r
```

## SGDP based analysis incl. 3D plot
documented R code to run  GLMM analysis
```
Rscript src/calcGLMM_sgdp.r
```

## ABC analysis
We reported three ABC analyses in the paper. 
\#1: Continues selection (internal run3v8)
\#2: Discontinues seleciton (internal run4v1)
\#3: Continues and Discontinues seleciton models combined (internal run5v1)
While \#1 and \#2 are independent analysis, in #3 uses the simulated data from #1 and #2 combined
Wrapping scripts include:
	- simulations and summary stats calculation
	- merging and filter of simulations
	- PLS transformation
	- ABC inference
	- plotting: posterior probabilities; prior and posterior distributions of parameters, cloud plots
	- power estimation
Requirements: msms needs to be set-up 

### 1: Continues selection (internal run3v8)
```
src/abc/continous_selection/fullABC_contSel.sh
Rscript src/abc/continous_selection/abcresults2table.r # get table with ABC estimates
```
### 2: Discontinues selection (internal run4v1)
```
src/abc/discontinous_selection/fullABC_discontSel.sh
```
### 3: Continous vs. Discontinues selection (internal run5v1)
```
src/abc/cont_vs_disc_selection/fullABC_contVSdiscontSel.sh
```
# SOM

## Recombination Rate Landscape: 1000g 2rdm Pop RecRate MEAN YRI,LWK;GBR,TSI;CHB,GIH
1000genomes Recombination rates based on illumina omni chips (~2M SNPs)
ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130507_omni_recombination_rates/
```
Rscript src/som_plot_recombLandscape.r
```

## Haploview data
based on vcf for target region: chr2:234806042-234948166 (TRPM8-refSeq +-20kb ) 
vcf transformed using public script: vcf_to_ped_convert.pl with minor-allele-frequency > 5%
*ped and *info files in data/

## Fst/XP-EHH plot 
XP-EHH and Fst calculated as described in manuscript, uploaded for plotting region to data/
```
Rscript src/som_plot_FstXPEHH.r
```

## SGDP data rs10166942 
extracted data in data/rs10166942_sgdp_daf_meta.txt
column names: c('ind','daf','pop','cont','ctry','src','sex','Latitude','Longitude','Coverage','HetRateAuto')
```
d=read.table('~/projects/trpm8_proj/vcf/sgdp/rs10166942_sgdp_daf_meta.txt',stringsAsFactors=F)
colnames(d)=c('ind','daf','pop','cont','ctry','src','sex','Latitude','Longitude','Coverage','HetRateAuto')
```

## African origin of rs10166942
based on pi calculation, and outlier haplotype characterization
inference (incl. table generation)
```
bash src/origin_rs10166942.sh
```
### plotting
```
Rscript src/plot_origin_rs10166942.r
```

