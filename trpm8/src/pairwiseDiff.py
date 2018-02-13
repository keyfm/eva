#!/usr/bin/env python

'''
Created on Mar 13, 2016

@author: felix_key
'''
import argparse, sys
import numpy as np


''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    This script caclulates the pairwise differences for each haplotype in phased vcf.
    If required it is possible to condition the haplotypes used by defining a position and an allele to be present.
    It outputs the mean and the median in pairwise differences and pairwise differences for each haplotype pair.
    Input vcf can be read from STDIN and be written to STDOUT (pipeable).
                                ''',
                                epilog="Questions or comments? --> felix_key@eva.mpg.de")
parser.add_argument('VCF', help='Phased VCF input file (script build for 1000gPIII) (explicitly its AA coding)', type=argparse.FileType('rt'), nargs='?',default=sys.stdin)
parser.add_argument('-p', dest='chr_position',help='Provide the position to condition the haplotypes on. Format: chr_position. (1-based as in VCF)',default=None)
parser.add_argument('-a', dest='allele',type=int,help='Integer of 0 or 1 that refers to desired allele. 0 is Ref and 1 is Alt. (Script changes later into DAF)')
parser.add_argument('-e',dest='excl_hap',type=argparse.FileType('r'),help='List of Haplotype IDs to exclude. Separated by newline char. e.g. NX12345_a',default=None)
parser.add_argument('-o',dest='out', nargs='?', type=argparse.FileType('w'),default=sys.stdout)
args = parser.parse_args()


'''Functions'''

def extractHapIDs(line):
    line = line.split("\t")
    ids = line[9:]
    hapIDs = []
    for myid in ids:
        hapIDs.append(myid+"_a")
        hapIDs.append(myid+"_b")
    return hapIDs
    
def extractAA(infoList):
    for i in infoList:
        if i.startswith("AA="):
            i = i.split("=")[1]
            anc = i.split("|")[0]
    return anc

def extractGTs(line,dct_i,dct_a):
    line = line.split("\t")
    # extract AA
    aa = extractAA(line[7].split(";"))
    # get allele list for SNP
    gts = line[9:]
    allele_list = []
    for gt in gts:
        x = gt.split("|")
        allele_list.extend([int(x[0]),int(x[1])])   
    # seed data into dictionaries
    snp_id = line[0]+"_"+line[1]
    dct_i[snp_id] = line[0:5]+[aa]
    dct_a[snp_id] = allele_list
    return (dct_i,dct_a)

def findIDX(lst, a):
    result = []
    for i, x in enumerate(lst):
        if x == a:
            result.append(i)
    return result

def reduceDCT(dct_allele,lst_AlleleIDX):
    # get dict with allele for remaining haps
    dct_reducedAllele = {}
    for idx in lst_AlleleIDX:
        for key in dct_allele.keys():
            dct_reducedAllele.setdefault(key,[]).append(dct_allele[key][idx])
    return dct_reducedAllele

def transposeDCT(dct_reducedAll,dct_i,ls_hap):
    dct_c = {}
    for idx,val in enumerate(ls_hap):
        # loop over sorted dictionary (sorting based on info dct pos field)
        for key,unusedVal in sorted(dct_i.items(), key=lambda e: int(e[1][1])):
            # build new dct for new keys and add value (-> haplotype)
            dct_c.setdefault(val,[]).append(dct_reducedAll[key][idx])
    return dct_c
    
def pairwiseDiff(dct_hap,numAlleles):
    ls_id_pair = []
    ls_pws = []
    num_red_haps = len(dct_hap.keys())
    ls_hapIDs = dct_hap.keys()
    ctr = 1
    # start from all hapIDs
    for idA in ls_hapIDs:
        # start to-be-compared-ID from second ID (and each run start an ID later)
        for i in range(ctr,num_red_haps):
            idB = ls_hapIDs[i]
            # to avoid that last pair has same IDs
            if idA != idB:
                pw_ctr = 0
                # for each pair loop over all entries (i.e. haplotype alleles)
                for idx in range(numAlleles):
                    if dct_hap[idA][idx] != dct_hap[idB][idx]:
                        pw_ctr = pw_ctr + 1
                ls_id_pair.append(idA+"-"+idB)
                ls_pws.append(pw_ctr)    
        # push second loop one count forward to avoid redundancies
        ctr = ctr + 1
        
    return (ls_id_pair,ls_pws)
    
    

'''Handle arguments'''
vcf = args.VCF
out = args.out
pos = args.chr_position
allele = args.allele
exclude_haps = args.excl_hap    
    
'''variables'''
# create two dictionaries with same keys, but values are either info or alleles
dct_redAltFreq = {}
dct_info = {}
dct_allele = {}

'''main code'''
for line in vcf:
    line=line.strip('\n')
    if line.startswith('#CHR'):
        hapIDs = extractHapIDs(line)
    elif not line.startswith('#'):
        (dct_info,dct_allele) = extractGTs(line,dct_info,dct_allele)

# extract indices where desired SNP has desired allele
if pos == None:
    # when no pos defined keep all haplotypes (thus all IDX's)
    lst_AlleleIDX = range(len(dct_allele[dct_allele.keys()[0]]))
else:
    if pos in dct_allele.keys():
        lst_AlleleIDX = findIDX(dct_allele[pos],allele)        
    else:
        sys.exit("Error. Provided positional argument not in VCF.")
           

# remove desired Haps (if specified) from list with to-keep-items (lst_AlleleIDX)
if exclude_haps != None:
    ls_hap = []
    idx_excl = []
    for line in exclude_haps:
        line=line.strip('\n')
        ls_hap.append(line)
    # get idx of excl. haps from actual haps file
    for idx,val in enumerate(hapIDs):
        if val in ls_hap:
            idx_excl.append(idx)
    # remove idx of the haps desired from allele_idx list (if present)
    for idx in idx_excl:
        if idx in lst_AlleleIDX:
            lst_AlleleIDX.remove(idx)
    if len(lst_AlleleIDX) == 0:
        sys.exit("Error. After removal of desired haplotype(s) no haplotypes left.")

      
# loop through each allele list and extract the allele count to obtain the alt allele frequency
dct_cond = reduceDCT(dct_allele,lst_AlleleIDX)

# remove Hap IDs that did not carry conditional allele before transposing (if no condition all IDs kept!)
ls_condHapIDs = [ hapIDs[i] for i in lst_AlleleIDX ]   

# transpose dictionalry to have key as hapID and value as haplotype (so we can later calculate pairwise diff)
dct_condHaps = transposeDCT(dct_cond,dct_info, ls_condHapIDs)

# calculate paiwise diff for each unique hap pair and get mean and median
(ls_id_pair,ls_pws) = pairwiseDiff(dct_condHaps,len(dct_allele.keys()))

# print header with watterson and info
out.write("## MeanPairwiseDiff:\t"+ str(np.average(ls_pws))+"\n")
out.write("## MedianPairwiseDiff:\t"+ str(np.median(ls_pws))+"\n")
# write PW count for each pair
out.write("#HapPair\tPairwiseDiffCT")
for i in range(len(ls_id_pair)):
    out.write("\n"+ls_id_pair[i]+"\t"+str(ls_pws[i]))


out.close()
vcf.close()
sys.exit(1)
    
            
            
            