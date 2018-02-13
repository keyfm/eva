#!/bin/bash
rid=${1}
v=${2}
pps=${3} # for plotting step implemented so far only (accounts if plots for eur,asi,all)
mod=${4}

## Things to be changed for each run:
# find_pls: 
#  - r script needs to be changed for desired parameters!
#  - boxcox transformation uses lambda -20->20 ... double check the Rplots output folder (versionX/) that peak is not lateral in plots!
# ABCtoolbox
#  - ABCest_template_wVctrl.input needs to be in run folder  /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/

## Things to know:
# The entire outputs overwrites any previous results (defined by run ID and version ID!)
# Make sure X11 connected...for the plots!
# The trasnformer script somehow always output ihsP1, rather than just the PLS components. I added an extra cut argument to get of those. Make the transformed files (obs,sims 5PLS) do carry the correct info. I double checked the transformation is not affected by that!
# After transformation I keep only PLS1-5. This is done correctly only when ihsP1 or ihsP2 is the first summary stat PLS transformed!

################################################
###### generate output folder / kill rpevious results
################################################
rm -rf /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}
mkdir -p /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}
mkdir -p /mnt/scratch/felix/trpm8/abc/run${rid}_out/transf/
mkdir -p /mnt/scratch/felix/trpm8/abc/obs/transf/run${rid}/
mkdir -p /home/felix_schulze/projects/trpm8_proj/pdf/abc/run${rid}/
################################################
###### find PLS
################################################
echo -e "Find PLS components."
for c in eur asi; do echo $c
 Rscript /mnt/scratch/felix/trpm8/abc/abc_scripts/run${rid}/find_pls_trpm8_run${rid}v${v}.r ${c} ${rid} ${v} ${mod}
 mv /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/Rplots.pdf /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/Rplots_lambda_${c}.pdf # otherwise lambda_Rplots for eur get overwritten by asi 
done
echo -e "Done."
# output: /mnt/scratch/felix/trpm8/abc/analysis/run",rid,"/version",v,"/find_pls_",c,"_allMod_v",v,".plsda

################################################
###### Transformation
################################################
# transformer prints all columns unchanged till the first statisti requested to be transformed. In case where ihsP1 is excluded and first stat is ihsP2 (affects cutting sims and obs to keep only PLS1-5
# get pos in header of ihsp1 and first stat of interest
# get first stat of interest in that run/version. Awful hack!
ss=$(head -n1 /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/find_pls_eur_allMod_v${v}.plsda |cut -f1) 
nums=$(zcat /mnt/scratch/felix/trpm8/abc/run${rid}_out/all_sumstat_5ParamBT_eur_sdn.tsv.gz |head -n1 |awk -v first=${ss} '{for(i=1;i<=NF;i++){if($i=="ihsP1"){ihsPOS=i};if($i==first){ssPOS=i};}print ihsPOS,ssPOS,ssPOS-ihsPOS+1,ssPOS-ihsPOS+5,ssPOS-ihsPOS+5+1,ssPOS-ihsPOS+5+5}')
ihsPOS=$(echo $nums |cut -d' ' -f1)
firstSS=$(echo $nums |cut -d' ' -f2)
diffPOS=$(echo $nums |cut -d' ' -f3)
obsEndPos=$(echo $nums |cut -d' ' -f4)
simStartPos=$(echo $nums |cut -d' ' -f5)
simEndPos=$(echo $nums |cut -d' ' -f6)

if [[ $ihsPOS == $firstSS ]]; then
    cutObs=1-5
    cutSims=1-10
else
    cutObs=${diffPOS}-${obsEndPos}
    cutSims=1-5,${simStartPos}-${simEndPos}
fi

### obs
echo -e "Transform obs."
plsda_file="/mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/find_pls_eur_allMod_v${v}.plsda"
for p in CEU GBR TSI FIN IBS CDX CHB CHS KHV JPT BEB GIH ITU PJL STU; do echo $p
 if [[ $p == 'CDX' ]]; then
  plsda_file="/mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/find_pls_asi_allMod_v${v}.plsda"
 fi
 /home/felix_schulze/src/ABCtoolbox/binaries/linux/transformer \
 ${plsda_file} \
 /mnt/scratch/felix/trpm8/abc/obs/obs_abcRun${rid}scp_YRI_${p}.tsv \
 /mnt/scratch/felix/trpm8/abc/obs/transf/run${rid}/obs_trans_abcRun${rid}v${v}scp_YRI_${p}.tsv \
 boxcox
 ## take only PLS 1-5 for inference (to actually reduce the stats!)
 sed -e 's/^[ \t]*//' /mnt/scratch/felix/trpm8/abc/obs/transf/run${rid}/obs_trans_abcRun${rid}v${v}scp_YRI_${p}.tsv |cut -f${cutObs} > /mnt/scratch/felix/trpm8/abc/obs/transf/run${rid}/obs_trans_abcRun${rid}v${v}scp_YRI_${p}_5PLS.tsv
 bgzip -f /mnt/scratch/felix/trpm8/abc/obs/transf/run${rid}/obs_trans_abcRun${rid}v${v}scp_YRI_${p}.tsv
done

### sims
echo -e "Transform sims."
for m in ntr sdn ssv; do echo $m
 for c in eur asi; do echo $c
  zcat /mnt/scratch/felix/trpm8/abc/run${rid}_out/all_sumstat_5ParamBT_${mod}_${c}_${m}.tsv.gz > /tmp/fmk_sims_${rid}_v${v}.tsv
  /home/felix_schulze/src/ABCtoolbox/binaries/linux/transformer \
  /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/find_pls_${c}_allMod_v${v}.plsda \
  /tmp/fmk_sims_${rid}_v${v}.tsv \
  /mnt/scratch/felix/trpm8/abc/run${rid}_out/transf/all_sumstat_5ParamBT_trans_${c}_${m}_run${rid}v${v}.tsv \
  boxcox
  rm /tmp/fmk_sims_${rid}_v${v}.tsv
 ## take only PLS 1-5 for inference (to actually reduce the stats!) and bgzip full file
 sed -e 's/^[ \t]*//' /mnt/scratch/felix/trpm8/abc/run${rid}_out/transf/all_sumstat_5ParamBT_trans_${c}_${m}_run${rid}v${v}.tsv |cut -f${cutSims} > /mnt/scratch/felix/trpm8/abc/run${rid}_out/transf/all_sumstat_5ParamBT_trans_${c}_${m}_run${rid}v${v}_5PLS.tsv
 bgzip -f /mnt/scratch/felix/trpm8/abc/run${rid}_out/transf/all_sumstat_5ParamBT_trans_${c}_${m}_run${rid}v${v}.tsv
 done
done


############################################################
############ 3. ABCtoolbox estimator
############################################################
echo -e "Run ABCtoolbox."
## get ABCtoolbox input template
cp /home/felix_schulze/projects/ifnan_project/abc/berkeley/reference/ceu_sdn3/ABCtoolbox /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/

for c in eur asi; do echo $c
 if [[ $c == 'eur' ]]; then
  pops="CEU GBR TSI FIN IBS"
 elif [[ $c == 'asi' ]]; then
  pops="CDX CHB CHS KHV JPT BEB GIH ITU PJL STU"
 fi
 for p in ${pops}; do echo $p
  rm -fr /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/YRI_${p}
  mkdir -p /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/YRI_${p}
  for m in ntr sdn ssv; do echo $m
    if [[ $m == 'sdn' ]]; then
     para_cols="1,2,3"
    elif [[ $m == 'ssv' ]]; then
     para_cols="1,3,5"
    elif [[ $m == 'ntr' ]]; then
     para_cols="1"
    fi
   cat /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/ABCest_template_wVctrl_in1M.input \
   |sed -e "s/-CONT-/${c}/g" -e "s/-MODEL-/${m}/g" -e "s/-POP1-/YRI/" -e "s/-POP2-/${p}/g" -e "s/-PARAM_COLS-/${para_cols}/" -e "s/-RID-/${rid}/g" -e "s/-VERSION-/${v}/g" \
   > /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/YRI_${p}/ABCest_YRI_${p}_${m}.input
   /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/ABCtoolbox /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/YRI_${p}/ABCest_YRI_${p}_${m}.input > /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/YRI_${p}/std_errout_${m}.txt 2>&1
  done
 done
done

### make table of marginal desnities
echo -e "Merge Tables."
rm -f /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/marginalDensities_ntr_sdn_ssv.txt
for p in CEU GBR TSI FIN IBS CDX CHB CHS KHV JPT BEB GIH ITU PJL STU; do echo $p
 cd /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/YRI_${p}/
 paste <(echo $p) <(cut -f2 out_ntr.marginalDensity.txt) <(cut -f2 out_sdn.marginalDensity.txt) <(cut -f2 out_ssv.marginalDensity.txt) \
 >> /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/version${v}/marginalDensities_ntr_sdn_ssv.txt
done

############################################################
############ 4. Plotting ModelChoice and PostEstimateDensities
############################################################
Rscript /mnt/scratch/felix/trpm8/abc/abc_scripts/ABCplot_PostSelMod_PostParam_PLS5clouds_wVctrl1.r ${rid} ${v} ${pps}
