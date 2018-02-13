#!/usr/bin/env python
'''
Created on Mar 24, 2016

@author: felix_key
'''
import argparse, sys
#from numpy.ma.core import abs



''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
                                 Script to modify a VCF file to msms file that includes all the info's needed to run abc_call_func.py for SumStat calculation!
                                 It requires that both input VCF's have the same number of individuals and the same variants encoded!
                                ''',
                                epilog="Questions or comments? --> felix_key@eva.mpg.de")
parser.add_argument('-vcf_pop1', help='Phased VCF input file for Pop1 (script build for 1000gPIII) (explicitly its AA coding)', type=argparse.FileType('rt'), nargs='?')
parser.add_argument('-vcf_pop2', help='Phased VCF input file for Pop2 (script build for 1000gPIII) (explicitly its AA coding)', type=argparse.FileType('rt'), nargs='?')
#parser.add_argument('-s', dest='sel_site',type=float,help='Integer that describes the location of the selected site. In the ouput that will be position 0.0')
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

def extract_derived_GTs(line,dct_i,dct_a):
    line = line.split("\t")
    # extract AA
    aa = extractAA(line[7].split(";"))
    # get allele list for SNP
    gts = line[9:]
    allele_list = []
         
    # AA == REF
    if line[3] == aa:
        for gt in gts:
            gt = gt.split(":")[0]
            x = gt.split("|")
            allele_list.extend([int(x[0]),int(x[1])])
    # AA == ALT        
    elif line[4] == aa:
        for gt in gts:
            gt = gt.split(":")[0]
            x = gt.split("|")
            for i in range(2):
                    if x[i] == "1":
                        x[i] = "0"
                    else:
                        x[i] = "1"
            allele_list.extend([int(x[0]),int(x[1])])
    # if none is met -> stop function and dont add new pos to any dictionary
    else:
        return (dct_i,dct_a)
                            
    # seed data into dictionaries
    snp_id = line[0]+"_"+line[1]
    dct_i[snp_id] = line[0:5]+[aa]
    dct_a[snp_id] = allele_list
    return (dct_i,dct_a)


'''Handle arguments'''
vcf_p1 = args.vcf_pop1
vcf_p2 = args.vcf_pop2
out = args.out
#sel_site = args.sel_site
        
'''variables'''
# create two dictionaries with same keys, but values are either info or alleles
dct_info_p1 = {}
dct_allele_p1 = {}
dct_info_p2 = {}
dct_allele_p2 = {}


'''main code'''
for line in vcf_p1:
    line=line.strip('\n')
    if line.startswith('#CHR'):
        hapIDs_p1 = extractHapIDs(line)
    elif not line.startswith('#'):
        (dct_info_p1,dct_allele_p1) = extract_derived_GTs(line,dct_info_p1,dct_allele_p1)


for line in vcf_p2:
    line=line.strip('\n')
    if line.startswith('#CHR'):
        hapIDs_p2 = extractHapIDs(line)
    elif not line.startswith('#'):
        (dct_info_p2,dct_allele_p2) = extract_derived_GTs(line,dct_info_p2,dct_allele_p2)

if len(hapIDs_p1) != len(hapIDs_p2):
    sys.exit("Error. Different number of individuals in vcf's.")

## Get sum stats for msms header

pos_id = []
for key, value in sorted(dct_info_p1.items(), key=lambda e: int(e[1][1])):
    pos_id.append( float(dct_info_p1[key][1]) )

length_seq = (pos_id[-1] - pos_id[0]) + 1
num_segsites = len(pos_id)

msms_string="ms -N X "+ str( len(hapIDs_p1) + len(hapIDs_p2) ) +" 1 -I 2 " + str(len(hapIDs_p1)) + " " + str(len(hapIDs_p2)) + " 0 -t <theta=4*N*mu*L> -r <rho=4*N*r*L> " + str(int(length_seq)) + " -I <number_of_populations> <sample_size_pop1> <sample_size_pop2> <migration_rate> -g <population> <growthrate> -n <population> <relative_size_of_pop1_relative_2_Ne> -n <population> <relative_size_of_pop2_relative_2_Ne> -ma x [m12] [m21] x -ej <divergencetime> 2 1 -N [Ne] -SFC -SI  <selection_time> <howmany populations> <frequency_in_pop1> <frequency_in_pop2> -Sc 0 <population> <selection_strength_in_populations> -Sc 0 <population> <selection_strength_in_population> -Sp <position_of_selected_site> -Smark\nrandom\n\n//\n"

## WRITE OUTPUT

out.write(msms_string)
out.write("segsites: "+str(num_segsites)+"\n")
out.write("positions: " + " ".join(map(str , pos_id)))


# write DAF header and data for each pos
# that should have really been better using an array in numpy. Clumpsy way to transform the dictioanries
for i in range(len(hapIDs_p1)):
    hap_list = []
    for key, value in sorted(dct_info_p1.items(), key=lambda e: int(e[1][1])):
        hap_list.append(dct_allele_p1[key][i])
    out.write("\n"+"".join(map(str,hap_list)))

for i in range(len(hapIDs_p2)):
    hap_list = []
    for key, value in sorted(dct_info_p2.items(), key=lambda e: int(e[1][1])):
        hap_list.append(dct_allele_p2[key][i])
    out.write("\n"+"".join(map(str,hap_list)))

out.close()
vcf_p1.close()
vcf_p2.close()

