'''
Created on Mar 23, 2016

@author: felix_key
'''

import numpy,sys,string
import random
import functions2calculate_sumstat as gss



def read_msms_output_simple(filestat1,input_msms,windowsize,recHS_flag,model_choice):
    '''Read in msms Header'''
    firstlines=[]
    for i in range(6):
        line = input_msms.readline()
        spline = line.split()
        if i==0:
            Lsequence = string.atoi(spline[ spline.index('-r')+2 ]) # extracts length of sim seq from -r option in msms command @ position -r+2; turns string into integer; NOTE overwrite Lsequence in case of recombination 
            isselection = spline.count("-SI")
            nochroms = string.atoi(spline[3]) #nochroms = string.atoi(spline[1])
            ms_cmd = numpy.array(spline) 
            if model_choice == 'sdn':
                para = ms_cmd[ [53,72,78] ] # seltime, selstrP1, selstrP2
                para = para.tolist()
                # backtransform time and selstr
                para[0] = float(para[0])*4*25*10000
                para[1] = float(para[1])/40000
                para[2] = float(para[2])/40000
                if para[0] >= 51000:
                    freq_p1 = float(1)/(2*14474)
                    freq_p2 = 0
                else:
                    freq_p1 = 0
                    freq_p2 = float(1)/(2*1861)
                para.append(freq_p1)
                para.append(freq_p2)
                my_params = para
            elif model_choice == 'ssv':
                para = ms_cmd[ [53,55,66] ] # seltime, initialFreqP1P2, selstrP2
                para = para.tolist()
                # backtransform time and selstr
                para[0] = float(para[0])*4*25*10000
                para[2] = float(para[2])/40000
                my_params = [para[0],0,para[2],para[1],para[1]]
            elif model_choice == 'ntr':
                para = ms_cmd[ [53,55,56] ] # seltime, initialFreqP1, initialFreqP2
                para = para.tolist()
                # backtransform time and selstr
                para[0] = float(para[0])*4*25*10000
                my_params = [para[0],0,0,para[1],para[2]]
            else:
                para = numpy.array(['nan','nan','nan','nan','nan'])
                my_params = para.tolist()
            outline = string.join(spline,sep=" ") + "\n"
            firstlines.append(outline)

        elif i>0 and i<=3:
            justincase=line.strip() + "\n"
            firstlines.append(justincase)

        elif i == 4:
            nsegs = string.atoi(line.split(" ")[1].strip("\n")) #integer with segregating sites; NOTE not changed below if recHS flag set!
            outsegsites = line.strip() + "\n"
            firstlines.append(outsegsites)

        elif i == 5:
            positions = spline[1:] # NOTE not changed below if recHS flag set!
            outpositions = line.strip() + "\n"
            firstlines.append(outpositions)

            posarray=numpy.asarray(map(string.atof,positions))
            
            #print "Selposition ",selposition
            #print "howmany ",numpy.where(posarray==0.5)

    
    '''Read in (simulated) Genetic Data'''
    # create array of 0 with NumChr==Lines and NumSegsites==Columns (as in msms out!) and fill it with simulated data
    themat = numpy.zeros([nochroms,nsegs])
    for i in range(nochroms):
        line = input_msms.readline()
        spline = line.split()
        for k in range(nsegs):
            # loop through each pos on each chr and fill array accordingly.
            # Note: derived alleles can be encoded 1 and [2:9,a-z], especially for sites that require a certain freq at a specific time. The number indicates their origin and can be re-coded as 1
            # Further info: https://github.com/delt0r/msms/issues/34
            if spline[0][k] == "0":
                themat[i,k] = string.atoi(spline[0][k])
            else:
                themat[i,k] = 1
                 
    # THEMAT contains all the haplotypes in pop1 and pop2..
    # close input file
    input_msms.close()
    
    '''mask extended simulated regions that mirror recombination hotspots'''
    if recHS_flag == True:
        mask_man = [(0.2519041,0.3062575) , (0.3215209,0.6244029) , (0.6503236,0.6746158) , (0.6960485,0.9370243)] # each tuple is a window that corresponds to an extended window that has to be excluded (peak 1-4). See also: trpm8_rec_rate.R l.180+
        posid2exclude = []
        for wdw in mask_man:
            idx = numpy.where( (posarray >= wdw[0]) & (posarray <= wdw[1]) )[0]
            posid2exclude.append(idx)
        posid2exclude = numpy.hstack(posid2exclude) # flatten to 1D array
        # exclude pos from posarray and from haplotype data
        posarray = numpy.delete( posarray , posid2exclude )
        themat = numpy.delete(themat, posid2exclude ,1)    
        Lsequence = 195136 # actual Lsequence
        # msms miraculously has sometimes positional duplicates in there. I do not know why (nothing in raised issues). Thus I break those analyses and continue to next iteration
        
    ''' msms miraculously has sometimes positional duplicates in there...I keep the first position and remove all others'''
    if len(numpy.unique(posarray)) != len(posarray):
        idx_sort = numpy.argsort(posarray)
        sorted_pos = posarray[idx_sort]
        vals, idx_start = numpy.unique(sorted_pos, return_index=True)
        counts_pos = numpy.split(idx_sort, idx_start[1:])
        duplicate_pos = [i[1:] for i in counts_pos if i.shape[0]>1] # NOTE: list comprehensions in python are much faster than for loops!
        duplicate_pos = numpy.hstack(duplicate_pos)
        posarray = numpy.delete( posarray , duplicate_pos )
        themat = numpy.delete(themat, duplicate_pos ,1) 
        
    '''get sum stat calculation rolling'''
    # calculates SNP freq (each col in array, if '1' for rows)
    freqs = numpy.sum(themat,0)
    # get pop info based on 2 pops in sim data with same num chr
    pop1size = int(nochroms/2.0)
    popnochroms=pop1size
        
    # constants required for TajD normalisation
    sfs_constants = get_constants(popnochroms)
    left = windowsize/2
    right= windowsize/2
    
    '''call sum stat function'''    
    # data to divide input into windows and calculate sum stats on those
    if recHS_flag == True:
        split_wdw = [ (0,0.1257441,65000,0) , (0.1257441,0.2517512,65136,0.18875) , (0.2517512,1,65000,0) ] #[0,1] start end roi; [2] Lsequence, [3] pos selsite (only 2nd window carries selsite, which corresponds to the middle position of the 2nd wdw!)
    else:
        split_wdw = [(234745003.0,234810003.0,65000,0) , (234810003.0,234875139.0,65136,234825093.0) , (234875139.0,234940139.0,65000,0)] # Note: Selsite is migraine SNP
    
    # first LD based stats which are calculated on full haplotype data
    selposition = split_wdw[1][3]
    #allstats = []
    allstats = gss.gets_LDtest_make_Hap_Map_files(themat, popnochroms, posarray, selposition,Lsequence,recHS_flag)
    allstats = my_params + allstats # add model parameter first!
    # get sumstats on flanking windows combined   
    wdw = split_wdw[1]
    idx = numpy.where( (posarray >= wdw[0]) & (posarray < wdw[1]) )[0] # gets all positions that fall in wdw1 (target window) and remove them subsequently
    wdw_posarray = numpy.delete( posarray , idx )
    wdw_themat = numpy.delete(themat, idx ,1)    
    Lsequence = wdw[2]    
    selposition = -100
    pop1pop2_sglPopStats,res_2PopStats = gss.get_all_stats(wdw_themat,popnochroms,wdw_posarray,sfs_constants,left,right,selposition,Lsequence)
    allstats = allstats + pop1pop2_sglPopStats + res_2PopStats  
    # get sumstats only on 2nd window (target window)
    wdw = split_wdw[1]
    idx = numpy.where( (posarray < wdw[0]) | (posarray >= wdw[1]) )[0] # gets all positions that fall NOT in wdw1
    wdw_posarray = numpy.delete( posarray , idx )
    wdw_themat = numpy.delete(themat, idx ,1)    
    Lsequence = wdw[2]
    if wdw[3] != 0:
        selposition = numpy.where(wdw_posarray==wdw[3])[0][0]
    else:
        selposition = -100
    pop1pop2_sglPopStats,res_2PopStats = gss.get_all_stats(wdw_themat,popnochroms,wdw_posarray,sfs_constants,left,right,selposition,Lsequence)
    allstats = allstats + pop1pop2_sglPopStats + res_2PopStats    
    
    outstring = string.join(map(str,allstats),sep="\t")
    #outfsts = string.join(map(str,fstpersnp),sep=" ")
    #outpi =  string.join(map(str,pipersnp2),sep=" ")
    
    # WRITE output files
    #fp1 = open(filestat1,"w")
    filestat1.write(outstring + "\n")
    filestat1.close()

    '''if fstfile != None:
        #fp2 = open(fstfile,"w")
        fstfile.write(outfsts + "\n")
        fstfile.close()

    if pifile != None:
    #fp3 = open(pifile,"w")
        pifile.write(outpi + "\n")
        pifile.close()'''
    
    return themat,selposition,freqs


def get_constants(popnochroms):
    # this defines constants as in R. Durrets' book  and Tajima '89
    # e1 and e2 required to calculate Variance of TajD (pi-watterson)
    vector = numpy.arange(1,popnochroms)
    n=popnochroms
    a1 = numpy.sum(1.0/vector) # corresponds to a in eq4.22 Clark/Hartl pg.176; Tajima '89 eq3
    a2 = numpy.sum(1.0/(vector**2)) # x**y means x^y ; Tajima '89 eq4
    b1 = float((popnochroms + 1))/(3.0*(popnochroms-1))
    b2 = (2.0*((popnochroms**2) + popnochroms + 3.0))/(9.0*popnochroms*(popnochroms-1.0))
    c1 = b1 - (1.0/a1)
    c2 = b2 - ((n+2.0)/(a1*n)) + (a2/(a1**2))
    e1 = c1/a1
    e2 = c2/((a1**2) + a2) # all checked
    #print "a1,a2,e1,e2 :"    
    #print a1,a2,e1,e2    
    
    return [a1,a2,e1,e2]

