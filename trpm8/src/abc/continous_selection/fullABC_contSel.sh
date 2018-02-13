########################################
##### ABC continous selection
########################################
#### Simulations
# I calculate all sumstats on the flanking windows (wdw1+3) combined and on the target window2.
# code made for SGE platform
qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run3/ABCtrpm8_run3.sge sdn eur
qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run3/ABCtrpm8_run3.sge ssv eur
qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run3/ABCtrpm8_run3.sge ntr eur

qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run3/ABCtrpm8_run3.sge sdn asi
qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run3/ABCtrpm8_run3.sge ssv asi
qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run3/ABCtrpm8_run3.sge ntr asi


## Combine all Sumstat files for only XPEHH stat being present
rid=3;mod=3 
for p in eur asi; do echo $p
 for m in ntr sdn ssv; do echo $m
  # write header for sumstat results
  echo "seltime_bt,selstrP1_bt,selstrP2_bt,freqssP1,freqssP2,ihsP1,ihsP2,XPEHH,wdwNonTrg_p1_thetaPi,wdwNonTrg_p1_thetaW,wdwNonTrg_p1_thetaH,wdwNonTrg_p1_TD,wdwNonTrg_p1_FWH,wdwNonTrg_p2_thetaPi,wdwNonTrg_p2_thetaW,wdwNonTrg_p2_thetaH,wdwNonTrg_p2_TD,wdwNonTrg_p2_FWH,wdwNonTrg_freq_selsite_P1,wdwNonTrg_freq_selsite_P2,wdwNonTrg_FstAll,wdwNonTrg_Fst_selsite,wdw2_p1_thetaPi,wdw2_p1_thetaW,wdw2_p1_thetaH,wdw2_p1_TD,wdw2_p1_FWH,wdw2_p2_thetaPi,wdw2_p2_thetaW,wdw2_p2_thetaH,wdw2_p2_TD,wdw2_p2_FWH,wdw2_freq_selsite_P1,wdw2_freq_selsite_P2,wdw2_FstAll,wdw2_Fst_selsite" |tr ',' '\t' > /mnt/scratch/felix/trpm8/abc/run${rid}_out/all_sumstat_5ParamBT_mod${mod}_${p}_${m}.tsv
  # combine sge results into one file with header
  for i in {1..1000}; do echo $i
   if [[ -e /mnt/scratch/felix/trpm8/abc/run${rid}_out/sge_out/ss_${m}_${p}_${i}.txt ]]; then
    awk '{OFS="\t";if($8>-100){ print $0 }}' /mnt/scratch/felix/trpm8/abc/run${rid}_out/sge_out/ss_${m}_${p}_${i}.txt >> /mnt/scratch/felix/trpm8/abc/run${rid}_out/all_sumstat_5ParamBT_mod${mod}_${p}_${m}.tsv
   fi
  done
 done
done

######################################################################
############## Get observation results
######################################################################
### Run real obs on ABCtrg region with ABC python script from simulation run3

for p in CEU GBR TSI FIN IBS CDX CHB CHS KHV JPT BEB GIH ITU PJL STU; do echo $p
    rid=3
    echo "ihsP1,ihsP2,XPEHH,wdwNonTrg_p1_thetaPi,wdwNonTrg_p1_thetaW,wdwNonTrg_p1_thetaH,wdwNonTrg_p1_TD,wdwNonTrg_p1_FWH,wdwNonTrg_p2_thetaPi,wdwNonTrg_p2_thetaW,wdwNonTrg_p2_thetaH,wdwNonTrg_p2_TD,wdwNonTrg_p2_FWH,wdwNonTrg_freq_selsite_P1,wdwNonTrg_freq_selsite_P2,wdwNonTrg_FstAll,wdwNonTrg_Fst_selsite,wdw2_p1_thetaPi,wdw2_p1_thetaW,wdw2_p1_thetaH,wdw2_p1_TD,wdw2_p1_FWH,wdw2_p2_thetaPi,wdw2_p2_thetaW,wdw2_p2_thetaH,wdw2_p2_TD,wdw2_p2_FWH,wdw2_freq_selsite_P1,wdw2_freq_selsite_P2,wdw2_FstAll,wdw2_Fst_selsite" |tr ',' '\t' > /mnt/scratch/felix/trpm8/abc/obs/obs_abcRun${rid}scp_YRI_${p}.tsv 
    vcf2msms.py \
    -vcf_pop1 <(zcat /home/felix_schulze/projects/trpm8_proj/vcf/ABCtrg_3x65kb_2_234745004_234940139_sglPops_allSNPs/YRI_65ind_ABCtrg_3x65kb_all.vcf.gz) \
    -vcf_pop2 <(zcat /home/felix_schulze/projects/trpm8_proj/vcf/ABCtrg_3x65kb_2_234745004_234940139_sglPops_allSNPs/${p}_65ind_ABCtrg_3x65kb_all.vcf.gz) \
    | /mnt/scratch/felix/trpm8/abc/abc_scripts/run${rid}/abc_call_func.py -r > /tmp/fmk.out
     cat /tmp/fmk.out |cut -f6- |awk '{if($1<-100){$1=0};if($2<-100){$2=0};if($3<-100){$3=0};OFS="\t";print $0}' >> /mnt/scratch/felix/trpm8/abc/obs/obs_abcRun${rid}scp_YRI_${p}.tsv
done



################################
###### ABCinference
################################
# PLS transformation
# ABCtoolbox analysis
# Plotting
/mnt/scratch/felix/trpm8/abc/abc_scripts/run3/ABCwrap_findPLS_transf_ABCtbox_plots_run3v8.sh 3 8 all mod3


################################
###### Power analysis
################################
/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh eur 3 8 sdn 10000
/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh eur 3 8 ssv 10000
/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh eur 3 8 ntr 10000

/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh asi 3 8 sdn 10000
/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh asi 3 8 ssv 10000
/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh asi 3 8 ntr 10000

### Estimate True Positives, False Positives, False Negatives
Rscript abc_power2table_r3v8.r
