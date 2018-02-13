########################################
##### ABC continous vs discontinous selection
########################################
# Infer power to discriminate evolutionary models with either continous or discontinous selection
#### Simulations

## get sims as soft links from run3v8 (*nTrunc*) and run4v1 (*wTrunc*)
for i in eur asi; do
  for m in ssv sdn ntr; do
    ln -s -T /mnt/scratch/felix/trpm8/abc/run3_out/all_sumstat_5ParamBT_mod3_${i}_${m}.tsv.gz /mnt/scratch/felix/trpm8/abc/run5_out/all_sumstat_5ParamBT_mod3_${i}_${m}_nTrunc.tsv.gz
    ln -s -T /mnt/scratch/felix/trpm8/abc/run4_out/all_sumstat_5ParamBT_mod3_${i}_${m}.tsv.gz /mnt/scratch/felix/trpm8/abc/run5_out/all_sumstat_5ParamBT_mod3_${i}_${m}_wTrunc.tsv.gz
  done
done

######################################################################
############## Get observation results
######################################################################
### Same as for continous selection model
## >> just cp
for p in CEU GBR TSI FIN IBS CDX CHB CHS KHV JPT BEB GIH ITU PJL STU; do echo $p
 cp /mnt/scratch/felix/trpm8/abc/obs/obs_abcRun3scp_YRI_${p}.tsv /mnt/scratch/felix/trpm8/abc/obs/obs_abcRun5scp_YRI_${p}.tsv
done



################################
###### ABCinference
################################
# PLS transformation
# ABCtoolbox analysis
# Plotting
/mnt/scratch/felix/trpm8/abc/abc_scripts/run5/ABCwrap_findPLS_transf_ABCtbox_plots_run5v1.sh 5 1 all mod3 > /mnt/scratch/felix/trpm8/abc/run4_out/logs/ABCwrap_findPLS_transf_ABCtbox_plots_run5v1.log 2>&1


################################
###### Power analysis
################################
## SGE suited for faster computation
rid=5
v=1
cap=10000

for m in sdn ssv ntr; do
  for p in eur asi; do
    for t in wTrunc nTrunc; do
      qsub /mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power_wTnT.sge $p $rid $v $m $cap $t
done; done; done

### Estimate True Positives, False Positives, False Negatives
Rscript abc_power2table_r5v1.r

