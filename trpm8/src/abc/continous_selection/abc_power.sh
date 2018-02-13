#!/bin/bash

# Arguments
c=${1} # eur or asi
# p=${2} # pop name. Make sure that agrees with $c !
rid=${2} # rid...run ID!
v=${3} # version of analysis. 
model=${4} # Model to do the power analysis: SDN or SSV or NTR
cap=${5} # number of repetitive runs to do    

## control of working directory
mkdir -p /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/power

rm -fr /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/power/${c}_run${rid}v${v}_${model}_cap${cap}
mkdir -p /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/power/${c}_run${rid}v${v}_${model}_cap${cap}/tmp

## generate file with cap-number of random pseudo-obs sets
# tee command to keep input file so one can later check whoch parameters give rise to sims that are incorrectly assigned (ntr model!)
sort -R /mnt/scratch/felix/trpm8/abc/run${rid}_out/transf/all_sumstat_5ParamBT_trans_${c}_${model}_run${rid}v${v}_5PLS.tsv | head -n ${cap} |tee /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/power/${c}_run${rid}v${v}_${model}_cap${cap}/rdmInputSim_${model}_${cap}.tsv |cut -f6- \
> /mnt/scratch/felix/trpm8/abc/run${rid}_out/transf/all_sumstat_5ParamBT_trans_${c}_${model}_5PLS_run${rid}v${v}_PowerRdm${cap}.tsv

bgzip -f /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/power/${c}_run${rid}v${v}_${model}_cap${cap}/rdmInputSim_${model}_${cap}.tsv

## build static .input file for ntr,sdn,ssv. After each run I kill the output and calculate again. All I keep is MD values
for test_against_mod in ntr sdn ssv; do
    # get columns to estimate params from
    if [[ $test_against_mod == 'sdn' ]]; then
	para_cols="1,2,3"
    elif [[ $test_against_mod == 'ssv' ]]; then
	para_cols="1,3,5"
    elif [[ $test_against_mod == 'ntr' ]]; then
	para_cols="1"
    fi
    # make input file
    cat /mnt/scratch/felix/trpm8/abc/abc_scripts/ABCest_template_power.input \
	|sed -e "s/-CONT-/${c}/g" -e "s/-MODEL-/${model}/g" -e "s/-POP1-/YRI/g" -e "s/-POP2-/${p}/g" -e "s/-PARAM_COLS-/${para_cols}/" -e "s/-RID-/${rid}/g" -e "s/-VERSION-/${v}/g" -e "s/-NUMCAP-/${cap}/g" -e "s/-TESTAGAINSTMOD-/${test_against_mod}/g" \
	> /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/power/${c}_run${rid}v${v}_${model}_cap${cap}/ABCest_template_power_v${v}_vs_${test_against_mod}.input
done

## start loop with predefined cap
cd /mnt/scratch/felix/trpm8/abc/analysis/run${rid}/power/${c}_run${rid}v${v}_${model}_cap${cap}
ctr=1
while [[ $ctr -le ${cap} ]]; do 
    sed "${ctr}q;d" /mnt/scratch/felix/trpm8/abc/run${rid}_out/transf/all_sumstat_5ParamBT_trans_${c}_${model}_5PLS_run${rid}v${v}_PowerRdm${cap}.tsv \
	|cat <(echo -e "LinearCombination_0\tLinearCombination_1\tLinearCombination_2\tLinearCombination_3\tLinearCombination_4") - > tmp/sample.obs
    for test_against_mod in ntr sdn ssv; do
	/mnt/scratch/felix/trpm8/abc/analysis/run${rid}/ABCtoolbox ABCest_template_power_v${v}_vs_${test_against_mod}.input > tmp/st_out.txt 2>&1
    done
    # save only MD values
    paste <(cut -f2 tmp/out_ntr.marginalDensity.txt) <(cut -f2 tmp/out_sdn.marginalDensity.txt) <(cut -f2 tmp/out_ssv.marginalDensity.txt) \
	>> marginal_densities_ntr_sdn_ssv.tsv
    rm tmp/*
    ((ctr = $ctr + 1))
done

# remove random sim input file
rm /mnt/scratch/felix/trpm8/abc/run${rid}_out/transf/all_sumstat_5ParamBT_trans_${c}_${model}_5PLS_run${rid}v${v}_PowerRdm${cap}.tsv

# Vioplots
Rscript /mnt/scratch/felix/trpm8/abc/abc_scripts/abc_power_plot.r ${c} ${rid} ${v} ${model} ${cap}

exit 0

# //inputfile for the program ABCestimator
# estimationType standard
# //file with the simulations. Consists of first the parameters and then the stats. A header line is required!
# simName /mnt/scratch/felix/trpm8/abc/run-RID-_out/transf/all_sumstat_5ParamBT_trans_-CONT-_-TESTAGAINSTMOD-_run-RID-v-VERSION-_5PLS.tsv
# //file with obnserved statistics. A Header is required with names corresponding to the stats in the simfile!
# obsName /mnt/scratch/felix/trpm8/abc/analysis/run-RID-/power/-CONT-_run-RID-v-VERSION-_-MODEL-_cap-NUMCAP-/tmp/sample.obs
# //columns containg parameters for which estimates will be produced
# params -PARAM_COLS-
# //number of simulations to estimate the GLM on
# numRetained     1000
# maxReadSims     1000000
# //the width of the diracpeaks, affecting the smoothing..
# diracPeakWidth 0.001
# //number of points at which to estimate posterior density
# posteriorDensityPoints 100
# //should the statistics be standardized? values: 1 / 0 (default)
# standardizeStats 1
# //should the prior be written in a file? values: 0 (default) / 1
# writeRetained 1
# obsPValue       1
# outputPrefix /mnt/scratch/felix/trpm8/abc/analysis/run-RID-/power/-CONT-_run-RID-v-VERSION-_-MODEL-_cap-NUMCAP-/tmp/out_-TESTAGAINSTMOD-.
