########################################
##### ABC discontinous selection
########################################
#### Simulations
# I calculate all sumstats on the flanking windows (wdw1+3) combined and on the target window2.
# code made for SGE platform
# selection stops 3000y ago
qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run4/ABCtrpm8_run4.sge sdn eur
qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run4/ABCtrpm8_run4.sge ssv eur
qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run4/ABCtrpm8_run4.sge ntr eur

qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run4/ABCtrpm8_run4.sge sdn asi
qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run4/ABCtrpm8_run4.sge ssv asi
qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/run4/ABCtrpm8_run4.sge ntr asi

### Combine all Sumstat files for only XPEHH stat being present
for p in eur asi; do echo $p
 rid=4;mod=3 ## CHANGE!!!!
 for m in ntr sdn ssv; do echo $m
  # write header for sumstat results
  echo "seltime_bt,selstrP1_bt,selstrP2_bt,freqssP1,freqssP2,ihsP1,ihsP2,XPEHH,wdwNonTrg_p1_thetaPi,wdwNonTrg_p1_thetaW,wdwNonTrg_p1_thetaH,wdwNonTrg_p1_TD,wdwNonTrg_p1_FWH,wdwNonTrg_p2_thetaPi,wdwNonTrg_p2_thetaW,wdwNonTrg_p2_thetaH,wdwNonTrg_p2_TD,wdwNonTrg_p2_FWH,wdwNonTrg_freq_selsite_P1,wdwNonTrg_freq_selsite_P2,wdwNonTrg_FstAll,wdwNonTrg_Fst_selsite,wdw2_p1_thetaPi,wdw2_p1_thetaW,wdw2_p1_thetaH,wdw2_p1_TD,wdw2_p1_FWH,wdw2_p2_thetaPi,wdw2_p2_thetaW,wdw2_p2_thetaH,wdw2_p2_TD,wdw2_p2_FWH,wdw2_freq_selsite_P1,wdw2_freq_selsite_P2,wdw2_FstAll,wdw2_Fst_selsite" |tr ',' '\t' > /mnt/scratch/felix/trpm8/abc/run${rid}_out/all_sumstat_5ParamBT_mod${mod}_${p}_${m}.tsv
  # combine sge results into one file with header
  for i in {1..100}; do echo $i
   if [[ -e /mnt/scratch/felix/trpm8/abc/run${rid}_out/sge_out/ss_${m}_${p}_${i}.txt ]]; then
    awk '{OFS="\t";if($8>-100){ print $0 }}' /mnt/scratch/felix/trpm8/abc/run${rid}_out/sge_out/ss_${m}_${p}_${i}.txt >> /mnt/scratch/felix/trpm8/abc/run${rid}_out/all_sumstat_5ParamBT_mod${mod}_${p}_${m}.tsv
   fi
  done
 done
done

rm /mnt/scratch/felix/trpm8/abc/run4_out/sge_out/*

######################################################################
############## Get observation results
######################################################################
### Same as for continous selection model
## >> just cp
for p in CEU GBR TSI FIN IBS CDX CHB CHS KHV JPT BEB GIH ITU PJL STU; do echo $p
 cp /mnt/scratch/felix/trpm8/abc/obs/obs_abcRun3scp_YRI_${p}.tsv /mnt/scratch/felix/trpm8/abc/obs/obs_abcRun4scp_YRI_${p}.tsv
done



################################
###### ABCinference
################################
# PLS transformation
# ABCtoolbox analysis
# Plotting
/mnt/scratch/felix/trpm8/abc/abc_scripts/run4/ABCwrap_findPLS_transf_ABCtbox_plots_run4v1.sh 4 1 all mod3 > /mnt/scratch/felix/trpm8/abc/run4_out/logs/ABCwrap_findPLS_transf_ABCtbox_plots_run4v1.log 2>&1


################################
###### Power analysis
################################
/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh eur 4 1 sdn 10000
/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh eur 4 1 ssv 10000
/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh eur 4 1 ntr 10000

/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh asi 4 1 sdn 10000
/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh asi 4 1 ssv 10000
/mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power.sh asi 4 1 ntr 10000

### Estimate True Positives, False Positives, False Negatives
Rscript abc_power2table_r4v1.r

