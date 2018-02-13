#! /usr/bin/python

'''
Created on Mar 23, 2016

@author: felix_key
'''

import argparse,sys
import functions2read_msmsout

''' positional and optional argument parser'''

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description='''\
    This script calculates the summary statistics for the ABC based on msms or msms-like simulation files.
    It requires two files containing functions. (a) functions2read_msmsout and (b) functions2calculate_sumstat. 
    Current output:
    5 params (backtransformed), LD tests, wdw0+2 combined, wdw1
                                ''',
                                epilog="Questions or comments? --> felix_key@eva.mpg.de")
parser.add_argument('-msms',dest='msms', help='simulated data for two populations each 100 chromosomes in msms format', type=argparse.FileType('rt'), nargs='?',default=sys.stdin)
parser.add_argument('-ss_out',dest='ss_out', nargs='?', help='Output file with basic sumstats.', type=argparse.FileType('w'),default=sys.stdout)
#parser.add_argument('-fst_out',dest='fst_out', nargs='?', help='output file with fsts per SNP. Tobedone: only prints last window currently.', type=argparse.FileType('w'),default=None)
#parser.add_argument('-pi_out',dest='pi_out', nargs='?', help='output file with pi per SNP.Tobedone: only prints last window currently.', type=argparse.FileType('w'),default=None)
parser.add_argument('-w', dest='window_size',type=int,help='Window size for SumStats. Extended SumStats are 2*windowsize. Default is 0, calculates the entire region.',default=0)
#parser.add_argument('-c', dest='cutoff_selscan',type=str,help='EHH decay cutoff for selscan.',default="0.05")
#parser.add_argument('-s', dest='sel_site',type=float,help='Floating value that describes the location of the selected site. If not specified the middle SNP from all supplied positions is used.',default=-1)
parser.add_argument('-r', dest='recHS',action='store_false',help='Default True. Means recombination hotspots are simulated which then will be excluded based on pre-defined mask. Set flag if that is not the case, e.g. with real data.',default=True)
parser.add_argument('-mc', dest='modelChoice',type=str,help='Specify model (sdn, ssv, ntr) so the code extracts the parameters from the msms command line. If no model selected it prints "nan".',default=None)
args = parser.parse_args()


input_msms = args.msms # single values
filename1 = args.ss_out
#filename2 = args.fst_out # fsts per snp
#filename3 = args.pi_out # pi per snps
windowsize = args.window_size # windowsize (+- windowsize/2); extended is 2*windowsize
#cutoff = args.cutoff_selscan
#selsite_loc = args.sel_site # position ID of selected site. Default 0.5 because it likely is this position in the sims!
recHS_flag = args.recHS
model_choice = args.modelChoice

# from imported package read function:read_msms_output_simple
themat,selposition,freqs = functions2read_msmsout.read_msms_output_simple(filename1,input_msms,windowsize,recHS_flag,model_choice)

