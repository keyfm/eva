#####################
### SOM Table 1 (A)
#####################
# Calculate pi for TRPM8 region (chr2_234810004_234875139)
for p in YRI LWK ESN GWD MSL ACB ASW CEU GBR TSI FIN IBS CDX CHB CHS KHV JPT BEB GIH ITU PJL STU CLM MXL PEL PUR;  do echo $p
 zcat /home/felix_schulze/projects/trpm8_proj/vcf/roi_2_234810004_234875139_sglPops_target_allSNPs/${p}_subset61_trpm8_target_roi_all.vcf.gz \
 |pairwiseDiff.py -p 2_234825093 -a 0 > /home/felix_schulze/projects/trpm8_proj/pairwiseDiff/migrHaps_allSNPs/fmk_${p}_pw.txt 
done

paste \
<(for p in ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI;  do echo $p; done) \
<(for p in ACB ASW BEB CDX CEU CHB CHS CLM ESN FIN GBR GIH GWD IBS ITU JPT KHV LWK MSL MXL PEL PJL PUR STU TSI YRI;  do \
   grep 'MeanPairwiseDiff' /home/felix_schulze/projects/trpm8_proj/pairwiseDiff/migrHaps_allSNPs/fmk_${p}_pw.txt  |cut -f2 ; done) \
| sort -k2,2n  |cat <(echo -e "pop\tTRPM8_trg_allSNPs") - \
> /home/felix_schulze/projects/trpm8_proj/pairwiseDiff/rs10166942_2_234825093_trg_allSNPs.tsv


#####################
### Get SNPs that are solely present on haplotypes that carry derived migraine allele with pairwise pi median>10
#####################
# pw_med10_cutoff_*.txt files are absed on the pairwiseDiff.py. Include haplotypes that have median pi > 10.
for p in YRI LWK ESN GWD MSL ACB ASW CEU GBR TSI FIN IBS CDX CHB CHS KHV JPT BEB GIH ITU PJL STU;  do \
 if [[ -e /home/felix_schulze/tmp/pw_med10_cutoff_${p}.txt ]]; then
  echo "process $p"
  zcat /home/felix_schulze/projects/trpm8_proj/vcf/roi_2_234810004_234875139_sglPops_target_allSNPs/${p}_subset61_trpm8_target_roi_all.vcf.gz > /tmp/fmk.bed     
  cat /tmp/fmk.bed |/home/felix_schulze/projects/trpm8_proj/scripts/extrctHapDAF_watterson_backup.py -p 2_234825093 -a 0 -k /home/felix_schulze/tmp/pw_med10_cutoff_${p}.txt |awk '{if($7>0){OFS="\t";print $1,$2-1,$2}}' > /tmp/fmk2.bed
  cat /tmp/fmk.bed |/home/felix_schulze/projects/trpm8_proj/scripts/extrctHapDAF_watterson_backup.py -p 2_234825093 -a 0 -e /home/felix_schulze/tmp/pw_med10_cutoff_${p}.txt |awk '{if($7>0){OFS="\t";print $1,$2-1,$2}}' > /tmp/fmk3.bed
  intersectBed -v -a /tmp/fmk2.bed -b /tmp/fmk3.bed > /home/felix_schulze/projects/trpm8_proj/pairwiseDiff/highPW_snps_only/highPWonly_${p}.bed
 fi
done

# Assess if SNPs linked to ancestral migraine allele
for p in YRI LWK ESN GWD MSL CEU GBR TSI FIN IBS CDX CHB CHS KHV JPT BEB GIH ITU PJL STU; do \
   zcat /home/felix_schulze/projects/trpm8_proj/vcf/roi_2_234810004_234875139_sglPops_target_allSNPs/${p}_subset61_trpm8_target_roi_all.vcf.gz \
   |/home/felix_schulze/projects/trpm8_proj/scripts/extrctHapDAF_watterson_backup.py -p 2_234825093 -a 1 |awk '{if($7>0){OFS="\t";print $1,$2-1,$2}}' \
   > /home/felix_schulze/projects/trpm8_proj/pairwiseDiff/highPW_snps_only/migrAncestralAllSNPs_${p}.bed
done

#####################
### Generate table that contains data generated above 
#####################
for p in YRI LWK ESN GWD MSL CEU GBR TSI FIN IBS CDX CHB CHS KHV JPT BEB GIH ITU PJL STU; do
 paste \
 <( echo $p) \
 <( cat /home/felix_schulze/projects/trpm8_proj/pairwiseDiff/highPW_snps_only/highPWonly_${p}.bed |wc) \
 <(intersectBed \
-a /home/felix_schulze/projects/trpm8_proj/pairwiseDiff/highPW_snps_only/highPWonly_${p}.bed \
-b /home/felix_schulze/projects/trpm8_proj/pairwiseDiff/highPW_snps_only/migrAncestralAllSNPs_${p}.bed |wc)
done > /home/felix_schulze/projects/trpm8_proj/pairwiseDiff/highPW_snps_only/highPWder_on_ancestral_perPop.txt


#####################
### SOM Table 1 (B)
#####################
# Calculate pi for TRPM8 region (chr2_234810004_234875139), but without haplotypes with median pi > 10 (and evidence of recombination)
paste \
<(for p in YRI LWK ESN GWD MSL ACB ASW CEU GBR TSI FIN IBS CDX CHB CHS KHV JPT BEB GIH ITU PJL STU CLM MXL PEL PUR;  do echo $p; done) \
<(for p in YRI LWK ESN GWD MSL ACB ASW CEU GBR TSI FIN IBS CDX CHB CHS KHV JPT BEB GIH ITU PJL STU CLM MXL PEL PUR;  do \
   if [[ -e /home/felix_schulze/tmp/pw_med10_cutoff_${p}.txt ]]; then
   zcat /home/felix_schulze/projects/trpm8_proj/vcf/roi_2_234810004_234875139_sglPops_target_allSNPs/${p}_subset61_trpm8_target_roi_all.vcf.gz \
    |pairwiseDiff.py -p 2_234825093 -a 0 -e /home/felix_schulze/tmp/pw_med10_cutoff_${p}.txt | tee  /home/felix_schulze/projects/trpm8_proj/pairwiseDiff/migrHaps_allSNPs/fmk_${p}_pw_noPw10med.txt \
    |grep 'MeanPairwiseDiff' |cut -f2 
   else
    zcat /home/felix_schulze/projects/trpm8_proj/vcf/roi_2_234810004_234875139_sglPops_target_allSNPs/${p}_subset61_trpm8_target_roi_all.vcf.gz \
    |pairwiseDiff.py -p 2_234825093 -a 0 | tee  /home/felix_schulze/projects/trpm8_proj/pairwiseDiff/migrHaps_allSNPs/fmk_${p}_pw_noPw10med.txt \
    |grep 'MeanPairwiseDiff' |cut -f2 
   fi
  done) \
| sort -k2,2n  |cat <(echo -e "pop\tTRPM8_trg_allSNPs_noPW10med") - \
> /home/felix_schulze/projects/trpm8_proj/watterson/rs10166942_2_234825093_trg_all_highPWmed10.tsv

