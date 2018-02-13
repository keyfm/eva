#!/usr/bin/env python

'''
Created on Feb 10, 2016

@author: felix_key
'''
import argparse, sys, random
from numpy.ma.core import abs

''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    This script extract a desired number of random haplotypes from a vcf based on the allele type of specific SNP and it calculates Watterson (Clark/Hartl pg.176 eq.4.23) (if >1 haplotypes).
    It provides Watterson estimate and the DAF of all SNPs that are variable in the remaining haplotypes (fixed ancestral and derived not included!; NOTE: When -k specified they are incl. so it prints DAF when e.g. only one haplotype repuested!).
    It provides a list of individuals that were used.
    Script can read vcf from STDIN and write output to STDOUT if not further defined. VCF has to be phased using pipe ('|') as GT delimiter. 
    For DAF BUT not Watterson calculation, only sites where used with AA == REF || AA == ALT. Number of SNPs used for Watterson and number of SNPs where DAF was inferred are indicated!
    Script also allows to exclude/keep specific haplotypes, which are encoded NAME-IN-VCF-HEADER and an added '_a' or '_b' depending if it is first or second allele respectively.
                                ''',
                                epilog="Questions or comments? --> felix_key@eva.mpg.de")
parser.add_argument('VCF', help='Phased VCF input file (script build for 1000gPIII) (explicitly its AA coding)', type=argparse.FileType('rt'), nargs='?',default=sys.stdin)
parser.add_argument('-n', dest='num',type=int,help='Define Number of Haplotypes. If not defined Max is used (encrypted within code by 0.',default=0)
parser.add_argument('-p', dest='chr_position',help='Provide the position to condition the haplotypes on. Format: chr_position. (1-based as in VCF)')
parser.add_argument('-a', dest='allele',type=int,help='Integer of 0 or 1 that refers to desired allele. 0 is Ref and 1 is Alt. (Script changes later into DAF)')
parser.add_argument('-e',dest='excl_hap',type=argparse.FileType('r'),help='List of Haplotype IDs to exclude. Separated by newline char. e.g. NX12345_a',default=None)
parser.add_argument('-k',dest='keep_hap',type=argparse.FileType('r'),help='List of Haplotype IDs to keep (exclusively!). Separated by newline char. e.g. NX12345_a',default=None)
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

def reducedAF(dct_allele,lst_AlleleIDX):
    # get dict with allele for remaining haps
    dct_reducedAllele = {}
    for idx in lst_AlleleIDX:
        for key in dct_allele.keys():
            dct_reducedAllele.setdefault(key,[]).append(dct_allele[key][idx])
    # get dict with Alt frequency
    dct_redAltFreq = {}
    n = len(lst_AlleleIDX)
    for key in dct_reducedAllele.keys():
        dct_redAltFreq[key] = round(float(sum(dct_reducedAllele[key]))/n,5)
    return dct_redAltFreq

def watterson(n,S):
    # following Clark/Hartl Principles of PopGen, pg.176 eq.:4.23
    a = 0
    for i in range(1,n): # NOTE: n-1 not needed becasue range is by default from x to y-1!!!
        a = a + 1/float(i)
    t = S/a
    return round(t,3)
        


'''Handle arguments'''
vcf = args.VCF
out = args.out
numHap = args.num
pos = args.chr_position
allele = args.allele
exclude_haps = args.excl_hap
keep_haps = args.keep_hap
        
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
if pos in dct_allele.keys():
    lst_AlleleIDX = findIDX(dct_allele[pos],allele)        
else:
    sys.exit("Error. Provided positional argument not in VCF.")
# control for max N if provided
if numHap != 0:
    # only shorten list for N if at least as many Haps are present then requested (of less Haps...just go with it!)
    if len(lst_AlleleIDX) >= numHap:
        lst_AlleleIDX = random.sample(lst_AlleleIDX,numHap)

# remove desired Haps (if specified) from list with to-keep-items
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

# keep desired Haps only (if specified) from list with to-keep-items
if keep_haps != None:
    ls_hap = []
    idx_keep = []
    for line in keep_haps:
        line=line.strip('\n')
        ls_hap.append(line)
    # get idx of keep haps from actual haps file
    for idx,val in enumerate(hapIDs):
        if val in ls_hap:
            idx_keep.append(idx)
    # make new list of idx for haps that carry conditional allele (if flag put) and are asked to be kept (exclusively...IDs not asked for will be removed)
    idx_keep_cond = []
    for idx in idx_keep:
        if idx in lst_AlleleIDX:
            idx_keep_cond.append(idx)
    lst_AlleleIDX = idx_keep_cond
    if len(lst_AlleleIDX) == 0:
        sys.exit("Error. None of the desired haplotype(s) to keep is present.")


# loop through each allele list and extract the allele count to obtain the alt allele frequency
dct_af = reducedAF(dct_allele,lst_AlleleIDX)
# change AF to DAF in dct_af and count number of positions (required for Watterson); use only sites with high confidence AA inference!
daf_ctr = 0
for key in dct_info.keys():
    '''# only remove nonvariable sites when keep flag not specified.(w/ -k we want to keep all sites and define presence absence outside of this script)
    if dct_af[key] == 0 or dct_af[key] == 1 and keep_haps == None:    
        del dct_info[key]
    else:'''
    if dct_info[key][4] == dct_info[key][5]:
        dct_af[key] = abs((dct_af[key]-1))
        daf_ctr = daf_ctr + 1
    elif dct_info[key][3] == dct_info[key][5]:
        daf_ctr = daf_ctr + 1
    else:
        del dct_info[key]        

# watterson can only be calculated with number of haps > 1 (when -k used maybe only one hap needed...than we do not need all that)
# following Clark/Hartl Principles of PopGen, pg.176 eq.:4.23
if len(lst_AlleleIDX) > 1:
    theta = watterson(len(lst_AlleleIDX),daf_ctr)
    # get haplotype IDs which were kept after pos filtering
    hap_names = [ hapIDs[i] for i in lst_AlleleIDX]    
    # print header with watterson and info
    out.write("## ThetaWatterson:\t"+ str(theta)+"\n")
    out.write("## N_Haplotypes:\t"+ str(len(lst_AlleleIDX))+"\n")
    out.write("## PolymorphicSites:\t"+ str(daf_ctr)+"\n")
    out.write("#X HaplotypeIDs:\t"+"\t".join(hap_names)+"\n")

# write DAF header and data for each pos
out.write("#chr\tpos\tdbSNP\tREF\tALT\tAA\tDAF")
for key, value in sorted(dct_info.items(), key=lambda e: int(e[1][1])):
    out.write("\n"+"\t".join(dct_info[key])+"\t"+str(dct_af[key]))            


out.close()
vcf.close()
if exclude_haps != None:
    exclude_haps.close()
if keep_haps != None:
    keep_haps.close()
sys.exit(1)

            
