## function for msms command buid
"%+%" = function(x,y){paste(x,y,sep="")}

complete_gsub = function(command,paramlist)
{
		for (p in names(paramlist))
		{
			command = gsub('[' %+% p %+% ']',paramlist[[p]],command,fixed=T)
		}
	command
}

## Command line options
args <- commandArgs(TRUE)
pop <- args[1]
model <- args[2]

## Read in Demography
#msmstemplate = "msms -ms <samplesize> <number of replicates> -t <theta=4*N*mu*L> -r <rho=4*N*r*L> <length_of_sequence=L> -I <number_of_populations> <sample_size_pop1> <sample_size_pop2> <migration_rate> -g <population> <growthrate> -n <population> <relative_size_of_pop1_relative_2_Ne> -n <population> <relative_size_of_pop2_relative_2_Ne> -ma x [m12] [m21] x -ej <divergencetime> 2 1 -N [Ne] -SFC -SI  <selection_time> <howmany populations> <frequency_in_pop1> <frequency_in_pop2> -Sc 0 <population> <selection_strength_in_populations> -Sc 0 <population> <selection_strength_in_population> -Sp <position_of_selected_site> -Smark"
# Sp flag sets the position of the selected site
# SFC conditions on the selected site not being lost
# Smark flag is to include the selected site in the ouput

## theta: -t <theta=4*N*mu*L> : 4 * 10000 * 0.0000000235 * 516923 = 485.9076
## mu based on gutentkunst SOM pg. 8: μ=0.0113·25/(2·6×106)=2.35×10−8 per generation.
## L: 516923 (including artificial recombination regions)
## recombination: -r <rho=4*N*r*L> 
## r = 7.441448e-07 cM/bp (average w/o recombination hotspot regions included) [trpm8_rec_rate.R l.90]
## -> 7.441448e-09 M/bp
## -> 4*10000*7.441448e-09*516923 = 153.8662

### SDN MODEL
if (model == 'sdn'){
  ## Select demography for SDN model
  if (pop == 'eur') {
    ## Africa vs Europe
    alpha = 151.712
    curr_Psize <- 33813
    msmstemplate = " -N 10000 -ms 244 1 -I 2 122 122 0 -t 485.9076 -r 153.8662 516923 -g 2 151.7119 -n 2 3.3813 -n 1 1.4474 -m 1 2 1 -m 2 1 1 -en 0.023 2 0.1861 -ema 0.023 2 x 6 6 x -ej 0.051 2 1 -en 0.148 1 0.731 -SFC -SI [seltime] 2 [initialfrequencyP1] [initialfrequencyP2] -Sc 0 1 [selstrengthP1] [halfselstrengthP1] 0 -Sc 0 2 [selstrengthP2] [halfselstrengthP2] 0 -Sp 0.18875 -Smark"
  } else if (pop == 'asi') {
    ## Africa vs Asia
    alpha = 191.541
    curr_Psize <- 45369.72
    msmstemplate = " -N 10000 -ms 244 1 -I 2 122 122 0 -t 485.9076 -r 153.8662 516923 -g 2 191.541 -n 2 4.536972 -n 1 1.4474 -m 1 2 0.312 -m 2 1 0.312 -ema 0.023 2 x 6 6 x -en 0.023 2 0.1861 -ej 0.051 2 1 -en 0.148 1 0.731 -SFC -SI [seltime] 2 [initialfrequencyP1] [initialfrequencyP2] -Sc 0 1 [selstrengthP1] [halfselstrengthP1] 0 -Sc 0 2 [selstrengthP2] [halfselstrengthP2] 0 -Sp 0.18875 -Smark"
  } else {
    stop("None of the required populations selected. (eur or asi)")
  }
  ## Build msms command line
  msmscommand = msmstemplate
  paramlist = list()
  paramlist[['seltime']] = runif(1,min=.03,max=.06)# between 30kya and 60kya
  time_of_selection = paramlist[['seltime']]
  ## calculation of SDN frequency 1/2N in pop2 given the demography	
  if (time_of_selection < 0.023){
    ## this never happens actually (>0.03!)...but maybe it changes in a repetition...so I leave it here
    N_time_sel = curr_Psize*exp(-alpha*time_of_selection)
    freq_time_selP2 = 1.0/(2.0*N_time_sel)
    freq_time_selP1 = 0
  } else if (time_of_selection >= 0.023 & time_of_selection < 0.051){
    N_time_sel = 1861.0
    freq_time_selP2= 1.0/(2.0*N_time_sel)
    freq_time_selP1 = 0
  } else{
    ## if seltime older than split, pop2 inexistent and subsequently frequency in pop2 = 0
    N_time_sel = 14474.0
    freq_time_selP2 = 0
    freq_time_selP1 = 1.0/(2.0*N_time_sel) #in SDN model it has freq 1/2N only if var appeared before outofAfrica
  }
  paramlist[['initialfrequencyP1']] = freq_time_selP1 #1/2N ...if before out of Africa
  paramlist[['initialfrequencyP2']] = freq_time_selP2 # 1 over the number of chromosome at the time of selection 
  paramlist[['selstrengthP1']] = runif(1,min=0,max=300)
  paramlist[['halfselstrengthP1']] = paramlist[['selstrengthP1']]*0.5
  paramlist[['selstrengthP2']] = runif(1,min=100,max=1000)
  paramlist[['halfselstrengthP2']] = paramlist[['selstrengthP2']]*0.5
  msmscommand = complete_gsub(msmscommand,paramlist)
  cat(c(msmscommand,"\n"))

### SSV MODEL
} else if (model == 'ssv') {
  ## Select demography for SDN model
  if (pop == 'eur') {
    ## Africa vs Europe
    msmstemplate = " -N 10000 -ms 244 1 -I 2 122 122 0 -t 485.9076 -r 153.8662 516923 -g 2 151.7119 -n 2 3.3813 -n 1 1.4474 -m 1 2 1 -m 2 1 1 -en 0.023 2 0.1861 -ema 0.023 2 x 6 6 x -ej 0.051 2 1 -en 0.148 1 0.731 -SFC -SI [seltime] 2 [initialfrequencyP1] [initialfrequencyP2] -Sc 0 2 [selstrengthP2] [halfselstrengthP2] 0 -Sp 0.18875 -Smark"
  } else if (pop == 'asi') {
    ## Africa vs Asia
    msmstemplate = " -N 10000 -ms 244 1 -I 2 122 122 0 -t 485.9076 -r 153.8662 516923 -g 2 191.541 -n 2 4.536972 -n 1 1.4474 -m 1 2 0.312 -m 2 1 0.312 -ema 0.023 2 x 6 6 x -en 0.023 2 0.1861 -ej 0.051 2 1 -en 0.148 1 0.731 -SFC -SI [seltime] 2 [initialfrequencyP1] [initialfrequencyP2] -Sc 0 2 [selstrengthP2] [halfselstrengthP2] 0 -Sp 0.18875 -Smark"
  } else {
      stop("None of the required populations selected. (eur or asi)")
  }
  ## Build msms command line
  msmscommand = msmstemplate
  paramlist = list()
  paramlist[['seltime']] = runif(1,min=.0209,max=.0509)# between 21kya and 51kya
  paramlist[['initialfrequencyP1']] = runif(1,min=1.0/15000.0,max=.2)
  paramlist[['initialfrequencyP2']] = paramlist[['initialfrequencyP1']]
  paramlist[['selstrengthP2']] = runif(1,min=1, max=1000)
  paramlist[['halfselstrengthP2']] = paramlist[['selstrengthP2']]*0.5
  msmscommand = complete_gsub(msmscommand,paramlist)
  cat(c(msmscommand,"\n")) # quotes in printed command cause a crash of msms and the vector index cat() solved both

### NTR MODEL
} else if (model == 'ntr') {
  ## Select demography for SDN model
  if (pop == 'eur') {
    ## Africa vs Europe
    msmstemplate = " -N 10000 -ms 244 1 -I 2 122 122 0 -t 485.9076 -r 153.8662 516923 -g 2 151.7119 -n 2 3.3813 -n 1 1.4474 -m 1 2 1 -m 2 1 1 -en 0.023 2 0.1861 -ema 0.023 2 x 6 6 x -ej 0.051 2 1 -en 0.148 1 0.731 -SFC -SI [seltime] 2 [initialfrequencyP1] [initialfrequencyP2] -Sc 0 -1 0 0 0 -Sp 0.18875 -Smark"
  } else if (pop == 'asi') {
    ## Africa vs Asia
    msmstemplate = " -N 10000 -ms 244 1 -I 2 122 122 0 -t 485.9076 -r 153.8662 516923 -g 2 191.541 -n 2 4.536972 -n 1 1.4474 -m 1 2 0.312 -m 2 1 0.312 -ema 0.023 2 x 6 6 x -en 0.023 2 0.1861 -ej 0.051 2 1 -en 0.148 1 0.731 -SFC -SI [seltime] 2 [initialfrequencyP1] [initialfrequencyP2] -Sc 0 -1 0 0 0 -Sp 0.18875 -Smark"
  } else {
    stop("None of the required populations selected. (eur or asi)")
  }
  ## Build msms command line
  msmscommand = msmstemplate
  paramlist = list()
  paramlist[['seltime']] = runif(1,min=.03,max=.06)
  time_of_selection = paramlist[['seltime']]
  ## calculation of de novo frequency 1/2N in pop2 given the demography	
  if (time_of_selection >= 0.023 & time_of_selection < 0.051){
    N_time_sel = 1861.0
    freq_time_selP2= 1.0/(2.0*N_time_sel)
    freq_time_selP1 = 0
  } else { ## if seltime older than split, pop2 inexistent and subsequently frequency in pop2 = 0
    N_time_sel = 14474.0
    freq_time_selP2 = 0
    freq_time_selP1 = 1.0/(2.0*N_time_sel) #in SDN model it has freq 1/2N only if var appeared before outofAfrica
  }
  paramlist[['initialfrequencyP1']] = freq_time_selP1 #1/2N ...if before out of Africa
  paramlist[['initialfrequencyP2']] = freq_time_selP2
  msmscommand = complete_gsub(msmscommand,paramlist)
  cat(c(msmscommand,"\n"))
} else {
  stop("None of the required models selected. (sdn, ssv, ntr)")
}

