'''
Created on Mar 23, 2016

@author: felix_key
'''

import numpy as np
import string
from subprocess import call
import tempfile



def get_all_stats(twopopmatrix,pop1size,posarray,sfs_constants,left,right,selposindex,Lsequence):
    # called: gss.get_all_stats(themat,popnochroms,positions,sfs_constants,left,right,selposition,Lsequence)
    pop1haplos = twopopmatrix[0:pop1size,]
    pop2haplos = twopopmatrix[pop1size:,]
    #posfloats = np.asarray(map(string.atof,posarray))  
    #left_ext = 2*left
    #right_ext = 2*right    
    # note that single_stats* is  list    
    
    single_stats_pop1,posint = gets_single_stats(pop1haplos,posarray,left,right,selposindex,Lsequence,sfs_constants)
    single_stats_pop2,posint = gets_single_stats(pop2haplos,posarray,left,right,selposindex,Lsequence,sfs_constants)
    
    pop1pop2_sglPopStats = single_stats_pop1 + single_stats_pop2 # this is thetaPi, thetaW,thetaH, TD, FWH per pop
    
    
    res_2PopStats = gets_twopop_stats(pop1haplos,pop2haplos,posarray,posint,selposindex) # freq_selsite P1, freq selsite P2, Fst All, Fst selsite
    
    #singlevalues = [str(selpos1freq[0]),str(selpos2freq[0]),str(FST),str(fstselsite[0]),str(fstaroundsite)]
    return pop1pop2_sglPopStats,res_2PopStats

    
def gets_LDtest_make_Hap_Map_files(twopopmatrix,pop1size,posarray, selposition,Lsequence,recHS_flag):
    pop1haplos = twopopmatrix[0:pop1size,]
    pop2haplos = twopopmatrix[pop1size:,]
    # get freqs in pop2 as floating points (for ratio's below) and get indices that fullfill criteria
    pop2freqs = np.asfarray( np.sum(pop2haplos[0:,],0) )
    pop1freqs = np.asfarray( np.sum(pop1haplos[0:,],0) )
    ihsMAF5pIDX_P2 = np.where( (pop2freqs/pop1size <= 0.05) | (pop2freqs/pop1size >= 0.95) ) # assume pop size for both pops is the same! >= <= required also by selscan
    ihsMAF5pIDX_P1 = np.where( (pop1freqs/pop1size <= 0.05) | (pop1freqs/pop1size >= 0.95) )
    xpDAF5pIDX = np.where(pop2freqs/pop1size <= 0.05)
     
    pos_ihs_p1 = np.delete( posarray , ihsMAF5pIDX_P1 )
    pos_ihs_p2 = np.delete( posarray , ihsMAF5pIDX_P2 )
    pos_xp = np.delete( posarray , xpDAF5pIDX )
    #selposition = 0.03023
    
    # calculate ihs for p1 IF target site still present
    # iHS pop1
    if selposition in pos_ihs_p1:
        # get hap file
        ihs_P1_input_hap = np.delete(pop1haplos,ihsMAF5pIDX_P1,1)
        tf_ihs_P1_input_hap = tempfile.NamedTemporaryFile()
        np.savetxt( tf_ihs_P1_input_hap.name , ihs_P1_input_hap ,delimiter=" ",fmt='%1i')      
        # get map file
        tf_ihs_P1_input_Map = tempfile.NamedTemporaryFile()
        if recHS_flag == True:
            phys_gen_pos = reset_phys_gen_pos_sims(pos_ihs_p1)
        else:
            phys_gen_pos = reset_phys_gen_pos_realData(pos_ihs_p1)
        idCtr = 1
        for i in range(len(phys_gen_pos)):
            tf_ihs_P1_input_Map.write( "2\t" + "id"+str(idCtr) + "\t"+ str(phys_gen_pos[i][0]) + "\t" + str(phys_gen_pos[i][1]) +"\n")
            idCtr = idCtr + 1
        tf_ihs_P1_input_Map.flush() # necessary to make the file bash accessable
        # get core pos for modified selscan
        coreSNP_idx = np.where(pos_ihs_p1==selposition)[0][0]
        # calculate ihs using selscan
        tf_ihs_P1_out = tempfile.NamedTemporaryFile()
        #print "/home/felix_schulze/projects/trpm8_proj/scripts/selscan_mod/selscan_sglP_2016_04_01_v3 --ihs --hap " + tf_ihs_P1_input_hap.name + " --map " + tf_ihs_P1_input_Map.name + " --loci " + str(coreSNP_idx) + " --out " + tf_ihs_P1_out.name
        call(["/home/felix_schulze/projects/trpm8_proj/scripts/selscan_mod/selscan_sglP_2016_04_01_v3","--ihs","--trunc-ok","--loci",str(coreSNP_idx),"--hap",tf_ihs_P1_input_hap.name,"--map",tf_ihs_P1_input_Map.name,"--out",tf_ihs_P1_out.name])
        #print "/Users/felix_schulze/Documents/Research/tools/selscan-master/bin/osx/selscan --ihs --hap " + tf_ihs_P1_input_hap.name + " --map " + tf_ihs_P1_input_Map.name + " --out " + tf_ihs_P1_out.name
        #call(["/Users/felix_schulze/Documents/Research/tools/selscan-master/bin/osx/selscan","--ihs","--trunc-ok","--hap",tf_ihs_P1_input_hap.name,"--map",tf_ihs_P1_input_Map.name,"--out",tf_ihs_P1_out.name])
        # extract ihs value and close (delete tmp files)
        selscan_out=open(tf_ihs_P1_out.name+".ihs.out")
        line = selscan_out.readline()
        if line != '':
            ihsP1 = line.split('\t')[5].strip()
        else:
            ihsP1 = -100.1
        call(["rm",tf_ihs_P1_out.name+".ihs.out",tf_ihs_P1_out.name+".ihs.log"])
        tf_ihs_P1_input_hap.close()
        tf_ihs_P1_input_Map.close()
        tf_ihs_P1_out.close()
    else:
        ihsP1 = -100.2
            
    # iHS pop2
    if selposition in pos_ihs_p2:
        # get hap file
        ihs_P2_input_hap = np.delete(pop2haplos,ihsMAF5pIDX_P2,1)
        tf_ihs_P2_input_hap = tempfile.NamedTemporaryFile()
        np.savetxt( tf_ihs_P2_input_hap.name , ihs_P2_input_hap ,delimiter=" ",fmt='%1i')  
        # get map file
        tf_ihs_P2_input_Map = tempfile.NamedTemporaryFile()
        if recHS_flag == True:
            phys_gen_pos = reset_phys_gen_pos_sims(pos_ihs_p2)
        else:
            phys_gen_pos = reset_phys_gen_pos_realData(pos_ihs_p2)
        idCtr = 1
        for i in range(len(phys_gen_pos)):
            tf_ihs_P2_input_Map.write( "2\t" + "id"+str(idCtr) + "\t"+ str(phys_gen_pos[i][0]) + "\t" + str(phys_gen_pos[i][1]) +"\n")
            idCtr = idCtr + 1
        tf_ihs_P2_input_Map.flush() # necessary to make the file bash accessable
        # get core pos for modified selscan
        coreSNP_idx = np.where(pos_ihs_p2==selposition)[0][0]
        # calculate ihs using selscan
        tf_ihs_P2_out = tempfile.NamedTemporaryFile()
        #print "/home/felix_schulze/projects/trpm8_proj/scripts/selscan_mod/selscan_sglP_2016_04_01_v3 --ihs --trunc-ok --hap " + tf_ihs_P2_input_hap.name + " --map " + tf_ihs_P2_input_Map.name + " --loci " + str(coreSNP_idx) + " --out " + tf_ihs_P2_out.name
        call(["/home/felix_schulze/projects/trpm8_proj/scripts/selscan_mod/selscan_sglP_2016_04_01_v3","--ihs","--trunc-ok","--loci",str(coreSNP_idx),"--hap",tf_ihs_P2_input_hap.name,"--map",tf_ihs_P2_input_Map.name,"--out",tf_ihs_P2_out.name])
        # extract ihs value and close (delete tmp files)
        selscan_out=open(tf_ihs_P2_out.name+".ihs.out")
        line = selscan_out.readline()
        if line != '':
            ihsP2 = line.split('\t')[5].strip()
        else:
            ihsP2 = -100.1
        call(["rm",tf_ihs_P2_out.name+".ihs.out",tf_ihs_P2_out.name+".ihs.log"])
        tf_ihs_P2_input_hap.close()
        tf_ihs_P2_input_Map.close()
        tf_ihs_P2_out.close()
    else:
        ihsP2 = -100.2
    
    # XPEHH pop1 vs pop2
    if selposition in pos_xp:
        # get hap file
        xp_P1_input_hap = np.delete(pop1haplos,xpDAF5pIDX,1)
        xp_P2_input_hap = np.delete(pop2haplos,xpDAF5pIDX,1)
        tf_xp_P1_input_hap = tempfile.NamedTemporaryFile()
        tf_xp_P2_input_hap = tempfile.NamedTemporaryFile()
        np.savetxt( tf_xp_P1_input_hap.name , xp_P1_input_hap ,delimiter=" ",fmt='%1i')
        np.savetxt( tf_xp_P2_input_hap.name , xp_P2_input_hap ,delimiter=" ",fmt='%1i') 
        # get map file
        tf_xp_input_Map = tempfile.NamedTemporaryFile()
        if recHS_flag == True:
            phys_gen_pos = reset_phys_gen_pos_sims(pos_xp)
        else:
            phys_gen_pos = reset_phys_gen_pos_realData(pos_xp)
        idCtr = 1
        for i in range(len(phys_gen_pos)):
            tf_xp_input_Map.write( "2\t" + "id"+str(idCtr) + "\t"+ str(phys_gen_pos[i][0]) + "\t" + str(phys_gen_pos[i][1]) +"\n")
            idCtr = idCtr + 1
        tf_xp_input_Map.flush() # necessary to make the file bash accessable
        # get core pos for modified selscan
        coreSNP_idx = np.where(pos_xp==selposition)[0][0]
        # calculate ihs using selscan
        tf_xp_out = tempfile.NamedTemporaryFile()
        #print "/home/felix_schulze/projects/trpm8_proj/scripts/selscan_mod/selscan_sglP_2016_04_01_v3 --xpehh --trunc-ok --ref " + tf_xp_P1_input_hap.name + " --hap " + tf_xp_P2_input_hap.name + " --map " + tf_xp_input_Map.name + " --loci " + str(coreSNP_idx) + " --out " + tf_xp_out.name
        call(["/home/felix_schulze/projects/trpm8_proj/scripts/selscan_mod/selscan_sglP_2016_04_01_v3","--xpehh","--trunc-ok","--loci",str(coreSNP_idx),"--ref",tf_xp_P1_input_hap.name,"--hap",tf_xp_P2_input_hap.name,"--map",tf_xp_input_Map.name,"--out",tf_xp_out.name])
        # extract xp value and close (delete tmp files)
        selscan_out=open(tf_xp_out.name+".xpehh.out")
        line = selscan_out.readline()
        line = selscan_out.readline()# output has header!
        if line != '':
            xp = line.split('\t')[7].strip()
        else:
            xp = -100.1
        call(["rm",tf_xp_out.name+".xpehh.out",tf_xp_out.name+".xpehh.log"])
        tf_xp_P1_input_hap.close()
        tf_xp_P2_input_hap.close()
        tf_xp_input_Map.close()
        tf_xp_out.close()
    else:
        xp = -100.2
    
    return [ihsP1,ihsP2,xp]
       
def reset_phys_gen_pos_sims(pos_arr):
    list_tuples = [] # a list with tuples, 0 index is genetic pos (in M) and 1 index is physical pos
    recrtBP = 7.441448e-07
    fullseq = 516923
    mask = [(0.2519041,0.3062575,28097) , (0.3215209,0.6244029,156567) , (0.6503236,0.6746158,12557) , (0.6960485,0.9370243,124566)] # each tuple is a window that corresponds to an extended window that has to be excluded (peak 1-4). See also: trpm8_rec_rate.R l.180+
    for i in range(len(pos_arr)):
        if pos_arr[i] < mask[0][1]:
            list_tuples.append( ( int(float(pos_arr[i])*fullseq)*recrtBP , int(float(pos_arr[i])*fullseq) ) )
        elif pos_arr[i] < mask[1][1]:
            list_tuples.append( ( int(float(pos_arr[i])*fullseq)*recrtBP , int(float(pos_arr[i])*fullseq)-mask[0][2] ) )
        elif pos_arr[i] < mask[2][1]:
            list_tuples.append( ( int(float(pos_arr[i])*fullseq)*recrtBP , int(float(pos_arr[i])*fullseq)-mask[0][2]-mask[1][2] ) )
        elif pos_arr[i] < mask[3][1]:
            list_tuples.append( ( int(float(pos_arr[i])*fullseq)*recrtBP , int(float(pos_arr[i])*fullseq)-mask[0][2]-mask[1][2]-mask[2][2] ) ) 
        else:
            list_tuples.append( ( int(float(pos_arr[i])*fullseq)*recrtBP , int(float(pos_arr[i])*fullseq)-mask[0][2]-mask[1][2]-mask[2][2]-mask[3][2] ) )             
    return list_tuples

def reset_phys_gen_pos_realData(pos_arr):
    list_tuples = [] # a list with tuples, 0 index is genetic pos (in M) and 1 index is physical pos
    start = pos_arr[0]
    recrtBP = 7.441448e-07
    mask = [(234875218,28097) , (234883108,156567) , (234896507,12557) , (234907586,124566)] # each tuple is midpoint of added distance2account for recHS, 2nd value length of insert based on trpm8_rec_rate.R l.181
    for i in range(len(pos_arr)):
        if pos_arr[i] < mask[0][0]:
            list_tuples.append( ( (pos_arr[i]-start)*recrtBP , int(pos_arr[i]) ) )
        elif pos_arr[i] < mask[1][0]:
            list_tuples.append( ( (pos_arr[i]-start+mask[0][1])*recrtBP , int(pos_arr[i]) ) )
        elif pos_arr[i] < mask[2][0]:
            list_tuples.append( ( (pos_arr[i]-start+mask[0][1]+mask[1][1])*recrtBP , int(pos_arr[i]) ) )
        elif pos_arr[i] < mask[3][0]:
            list_tuples.append( ( (pos_arr[i]-start+mask[0][1]+mask[1][1]+mask[2][1])*recrtBP , int(pos_arr[i]) ) )
        else:
            list_tuples.append( ( (pos_arr[i]-start+mask[0][1]+mask[1][1]+mask[2][1]+mask[3][1])*recrtBP , int(pos_arr[i]) ) )     
    return list_tuples
  
def gets_twopop_stats(pop1haplos,pop2haplos,posarray,posint,selpos):
    pop2freqs = np.sum(pop2haplos[0:,],0)
    pop1freqs = np.sum(pop1haplos[0:,],0)
    no_chroms = pop1haplos.shape[0]
    no_segsites = pop1haplos.shape[1]
        
    FST,wholenum,wholeden = gets_fst_pi_perSNP(pop1freqs,pop2freqs,no_chroms,no_segsites)
    # calculate Fst and freq for selected site (only the case in wdw2!)
    if selpos != -100:
        selpos1freq = pop1freqs[selpos]/no_chroms
        selpos2freq = pop2freqs[selpos]/no_chroms
        if wholeden[selpos] != 0.0:
            fstselsite = [wholenum[selpos]/wholeden[selpos]][0]
        else:
            fstselsite=0.0
    else:
        selpos1freq = 'NA'
        selpos2freq = 'NA'
        fstselsite = 'NA'
    
    #fstaroundsite = np.sum(wholenum[posint])/np.sum(wholeden[posint])
    #fstpersnp,modifiednum,modifiedden = gets_fsts_per_snp(posarray,wholenum,wholeden)
    return [str(selpos1freq),str(selpos2freq),str(FST),str(fstselsite)]

def gets_single_stats(pophaplos,posarray,left,right,selpos,Lsequence,sfs_constants):
    if right == 0.0:
        flaggo ="all"
        freq_stats,posint = gets_stats_inside_interval(flaggo,pophaplos,posarray,sfs_constants)
    else:
        flaggo ="subset"
        freq_stats,posint = gets_stats_inside_interval(flaggo,pophaplos,posarray,sfs_constants,left,right,selpos,Lsequence)

    return freq_stats,posint

def gets_stats_inside_interval(flaggo,haplotmat,posarray,sfs_constants,left='NA',right='NA',selpos='NA',Lsequence='NA'):
    if flaggo == "subset":
        pos = string.atof(posarray[selpos])
        allpos = np.asarray(map(string.atof,posarray))
        # to define window in data its necessary to differentiated between observed data (positions > 1) or simulated data (positions range(0,1))
        if allpos[0] > 1:
            intervalL = left
            intervalR = right
        else:
            intervalL = left/float(Lsequence)
            intervalR = right/float(Lsequence)
        # get array with left/right roi boundary length-posarray-Times 
        xleft = np.asarray(len(posarray)*[pos-intervalL])
        xright = np.asarray(len(posarray)*[pos+intervalR])
        # boolean arrays for both tails (NOTE: > and < BUT NOT >= and <=!)
        toleft = np.greater(allpos,xleft)
        toright = np.less(allpos,xright)
        # obtains an array of same length but TRUE are only the overalpping TRUE's 
        posint = np.logical_and(toleft,toright)
        # get array of pos IDs that are within roi
        posinside = posarray[posint]
        
        windowmat = haplotmat[0:,posint]
        singlepopstats = get_sfs(windowmat,sfs_constants)
    else:
        windowmat = np.copy(haplotmat)

        # singlepopstats = [thetaPi,thetaW,thetaH,TD,FWH]
        singlepopstats = get_sfs(windowmat,sfs_constants)
        posint =[]
    # singlepopstats is a list
    return singlepopstats,posint

def get_sfs(windowmat,sfs_constants):
    # sfs_constants =[a1,a2,e1,e2] as in Rick Durretts book
    sumder = np.sum(windowmat[0:,],0)
    sizes = windowmat.shape
    chroms = np.arange(sizes[0])
    sfs = np.asarray([len(np.where(sumder==x)[0]) for x in chroms])

    thetaPi,thetaW,thetaH,TD,FWH = get_all_singlepop_stats(sfs,sfs_constants)
    allstats = [thetaPi,thetaW,thetaH,TD,FWH]
    return allstats 

def get_all_singlepop_stats(sfs,sfs_constants):
    n = len(sfs)
    vector = np.arange(1,n)
    reversedvec = vector[::-1]    
    # calculate number of pairwise combinations for n chromosomes
    nchoose2 = n*(n-1.0)/2.0 
    
    a1 = sfs_constants[0]
    # calc ss. Note: sfs first value is count of zeros...exclude this value by sfs[1:]
    thetaW = float(np.sum(sfs[1:]))/a1 # checked;
    thetaPi = (1.0/nchoose2)*np.sum(vector*reversedvec*sfs[1:]) # checked; thetaPi:all pairwise diffs per numPossPairsChr
    thetaH = (1.0/nchoose2)*np.sum(vector*vector*sfs[1:]) # checked
    #print "SFS: ",sfs
    NoSegsites = np.sum(sfs[1:])
    TD = tajimasD(thetaPi,thetaW,NoSegsites,sfs_constants)
    FayWuH = Fay_and_Wu_H(thetaH,thetaPi)
    return thetaPi,thetaW,thetaH,TD,FayWuH

def tajimasD(thetaPi,thetaW,Noss,sfs_constants):
    # sfs_constants =[a1,a2,e1,e2] as in Rick Durretts book
    denominator = sfs_constants[2]*Noss + sfs_constants[3]*Noss*(Noss-1)
    # denominator is the expected variance based on Tajima '89
    #print "thetaPi: ",thetaPi
    #print "thetaW: ",thetaW
    #print "Number of Segsites: ",Noss
    #print "constants ",sfs_constants
    #print "denominator: ",denominator
    if denominator != 0.0:
        TD = (thetaPi-thetaW)/np.sqrt(denominator)
    else:
        TD = -100.0
    #print "Tajima's D: ",TD
    return TD
    
def Fay_and_Wu_H(thetaH,thetaPi):
    FWH = thetaH-thetaPi
    return FWH
    

def get_2D_sfs(pop1freqs,pop2freqs,selpos,n1):
    twodsfs = np.zeros(81*81)
    #maxdiff = np.max(
    for i in range(len(pop1freqs)):
        freq1 = pop1freqs[i]
        freq2 = pop2freqs[i]
        if i<10:
            print "freq1 ",freq1
            print "freq2 ",freq2

        theindex = freq1*(n1 +1) + freq2
        twodsfs[int(theindex)] +=1.0

    return twodsfs

def gets_fsts_per_snp(posarray,wholenum,wholeden):
    denzeros = np.where(wholeden==0.0)[0]
    #denzeros are all the entries that have ZERO in the denominator
    for i in denzeros:
        #print "do i go into this loop ? HOPE NOT!"
        wholenum[i]=0.0; wholeden[i]=1.0
    # NOTE: for some individual snps, the denominator can be zero.
    # in that case, I am letting the fst=0
    fstpersnp = wholenum/wholeden
    #subsetpos = posarray[denzeros]
    return fstpersnp,wholeden,wholenum

def gets_fst_pi_perSNP(pop1freqs,pop2freqs,no_chroms,no_segsites):
    no_ind = no_chroms/2.0 
    sizearray = no_ind*np.ones(no_segsites)
    FST,wholenum,wholeden = FST_reynolds_num_den(sizearray,sizearray,pop1freqs,pop2freqs)
    #pi_pop1,pi_pop2 = pi_per_snp(pop1freqs,pop2freqs,no_chroms)

    return FST,wholenum,wholeden 

def FST_reynolds_num_den(noT,noH,T1,H1):
    # Reynolds/Weir/Cockerham '83: weigthed measure of PopDiff
    # weighting of AF difference and difference in number of chr between pops
    
    # Numerator
    # NOTE: noT and noH are the number of individuals, not the number of chromosomes
    T1 = T1/(2.0*noT); H1 = H1/(2.0*noH)
    sharednum = noT*(2.0*T1-2.0*T1**2) + noH*(2.0*H1-2.0*H1**2)
    num1 = (H1-T1)**2
    fracnum = (noT + noH)*(sharednum)
    fracden = 4.0*noT*noH*(noT + noH - 1.0)
    frac = fracnum/fracden
    
    wholenum = num1 - frac

    # Denominator
    denfracnum = (4.0*noT*noH -noT-noH)*sharednum
    denfrac = denfracnum/fracden
    
    wholeden = num1 + denfrac
    
    # set negative values to 0
    negvals = np.where(wholenum<0)
    wholenum[negvals]=0.0
    wholeden[negvals]=0.0

    if np.sum(wholeden) != 0.0:       
        FST = np.sum(wholenum)/np.sum(wholeden)
    else:
        FST=0.0
    
    return FST,wholenum,wholeden

def pi_per_snp(pop1freqs,pop2freqs,no_chroms):
    constant = 2.0/(no_chroms*(no_chroms-1.0))
    pi_pop1 = pi_per_pop(pop1freqs,constant,no_chroms)
    pi_pop2 = pi_per_pop(pop2freqs,constant,no_chroms)
    return pi_pop1,pi_pop2

def pi_per_pop(freqs,constant,no_chroms):
    narray = no_chroms*np.ones(len(freqs))
    nminusi = narray-freqs
    pi = constant*freqs*nminusi
    return pi

def avg_distance_between_positions(positions):
    pos = np.asarray(positions)
    dist=[]
    for i in range(1,len(pos)):
        if i==1:
            dist.append(pos[i]-pos[0])
        else:
            dist.append(pos[i]-pos[i-1])
    
    return dist

            
