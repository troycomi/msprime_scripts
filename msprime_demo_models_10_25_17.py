from __future__ import print_function
import msprime
import math
import numpy as np
import pybedtools
import random
import scipy
import subprocess
import glob
import gzip
from scipy import spatial
from optparse import OptionParser
#
#usage = "usage: python <script> -p nonAfr -o Tenn -s 1 -i 1 -n 0.02 -d 0.01 -t 350"
#parser = OptionParser(usage=usage)
#parser.add_option("-p", "--population", action = "store", type = "string", dest = "pop", default="nonAfr", help="call introg_haplotypes in EUR, EAS, or all nonAfr; default=nonAfr")
#parser.add_option("-o", "--outdir", action = "store", type = "string", dest = "outdir", help="output directory name, model name")
#parser.add_option("-s", "--seed", action = "store", type = "int", dest = "seed", default=1, help="Set random seed for replicate chromosomes; default=1")
#parser.add_option("-i", "--introgress_pulses", action = "store", type = "int", dest = "pulses", default=1, help="Set number of introgression pulses; default=1")
#parser.add_option("-n", "--neand1_admixture_proportion", action = "store", type = "float", dest = "n1_admix_prop", default=0.02, help="Set N1 admixture proportion; default=0.02")
#parser.add_option("-d", "--deni_admixture_propotion", action = "store", type = "float", dest = "n2_admix_prop", default=0.02, help="Set N2 admixture proportion; default=0.01")
#parser.add_option("-t", "--time_N1_N2_split", action = "store", type = "int", dest = "t_n1_n2", default=350, help="Set N1 N2 split time in kya; default=350")
#parser.add_option("-h", "--haplotype_or_stats", action="store", type = "string", dest = "haplo", default="F4Dstat", help="Pick haplo or F4Dstat output; default=F4Dstat")
##parser.add_option("-t", "--taskid", action = "store", type = "string", dest = "taskid")

#length = 1e6	# simulate 1Mb chromosomes for faster F4 and D calculation
recombination_rate = 1e-8  # 2e-8 original from Rajiv, not sure where it is from.
mutation_rate = 1e-8 # see Shendure and Akey, 2015


def Tenn_demography(S_N1, S_N2, S_AF, S_EU, S_AS, pulses, seed, n1_admix_prop, n2_admix_prop, outdir, t_n1_n2, haplo, length):
    N_A = 7310 # effective population size used for scaling
    N_N1 = 1000 # Altai/Vindija lineage effective population size
    N_N2 = 1000 # Eastern Neandertal effective population size
    N_B = 1861 # Population size after EUR + ASN join
    N_B_BN = 100 # Population size of EUR_ASN bottleneck (following initial admixture event)
    N_AF = 14474
    N_EU0 = 1032
    N_AS0 = 554
    N_CH = 1000
    N_DE = 1000
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_MH_CH = 7e6 / generation_time # Chimp_Human lineage split
    T_MH_N = 700e3 / generation_time # Neanderthal / modern human split
    T_DE_N = 500e3 / generation_time # Denisovan / Neandertal split time
    T_N1_N2 = ((t_n1_n2) * 1e3) / generation_time # Altai/Vindija and Eastern Neandertal popualation split, default set to 350kya
    T_AF = 148e3 / generation_time # African population expansion
    T_B = 100e3 / generation_time # out of Africa
    T_PULSE1 = 55e3 / generation_time # Neandertal introgression pulse 1
    #T_B_BN = 45e3 / generation_time # A 3rd Bottleneck following shortly after initial admixture event
    T_EU_AS = 23e3 / generation_time # European / East Asian split
    T_PULSE2 = 18e3 / generation_time # Neandertal introgression pulse 2
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 15e-5
    m_AF_EU = 2.5e-5
    m_AF_AS = 0.78e-5
    m_EU_AS = 3.11e-5
    m_N1_N2 = 0
    m_N1_AF = 0
    m_N1_EU = 0
    m_N1_AS = 0
    m_N2_AF = 0
    m_N2_EU = 0
    m_N2_AS = 0
    m_N1_CH = 0
    m_N2_CH = 0
    m_AF_CH = 0
    m_EU_CH = 0
    m_AS_CH = 0
    m_DE_N1 = 0
    m_DE_N2 = 0
    m_DE_AF = 0
    m_DE_EU = 0
    m_DE_AS = 0
    m_DE_CH = 0
    # Set Percent introgression level
    m_PULSE1 = n1_admix_prop # an instantaneous movement of population from EU_AS to Neand1 (backward in time)
    m_PULSE2 = n2_admix_prop
    #print(m_PULSE1)
    
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
        initial_size = N_N1), # population 0
        msprime.PopulationConfiguration(
        initial_size = N_N2), # population 1
        msprime.PopulationConfiguration(
        initial_size = N_AF), # population 2
        msprime.PopulationConfiguration(
        initial_size = N_EU, growth_rate = r_EU), # population 3
        msprime.PopulationConfiguration(
        initial_size = N_AS, growth_rate = r_AS), # population 4
	msprime.PopulationConfiguration(
	initial_size = N_CH), #popultaion 5
	msprime.PopulationConfiguration(
	initial_size = N_DE)  #population 6
    ]
    # specify desired sample
    N1_sample = [msprime.Sample(population = 0, time = 100e3 / generation_time)] * S_N1
    N2_sample = [msprime.Sample(population = 1, time = 100e3 / generation_time)] * S_N2
    AF_sample = [msprime.Sample(population = 2, time = 0)] * S_AF
    EU_sample = [msprime.Sample(population = 3, time = 0)] * S_EU
    AS_sample = [msprime.Sample(population = 4, time = 0)] * S_AS
    CH_sample = [msprime.Sample(population = 5, time = 0)] * 1
    DE_sample = [msprime.Sample(population = 6, time = 100e3/ generation_time)] * 1
    samples = N1_sample + N2_sample + AF_sample + EU_sample + AS_sample + CH_sample + DE_sample
    # set up the initial migration matrix
    migration_matrix = [
        [      0, m_N1_N2, m_N1_AF, m_N1_EU, m_N1_AS, m_N1_CH, m_DE_N1],
        [m_N1_N2,       0, m_N2_AF, m_N2_EU, m_N2_AS, m_N2_CH, m_DE_N2],
        [m_N1_AF, m_N2_AF,       0, m_AF_EU, m_AF_AS, m_AF_CH, m_DE_AF],
        [m_N1_EU, m_N2_EU, m_AF_EU,       0, m_EU_AS, m_EU_CH, m_DE_EU],
        [m_N1_AS, m_N2_AS, m_AF_AS, m_EU_AS,       0, m_AS_CH, m_DE_AS],
	[m_N1_CH, m_N2_CH, m_AF_CH, m_EU_CH, m_AS_CH,       0, m_DE_CH],
	[m_DE_N1, m_DE_N2, m_DE_AF, m_DE_EU, m_DE_AS, m_DE_CH,       0]
    ]
    one_pulse = [
        msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (2, 3)), # migration between "B" and Africa begins
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)),
        msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
	#msprime.PopulationParametersChange(time = T_B_BN, initial_size = N_B_BN, growth_rate = 0, population_id = 3),  # population bottleneck of B begins
	#msprime.PopulationParametersChange(time = T_B_BN + 20, initial_size = N_B, growth_rate = 0, population_id = 3), # population bottleneck ends shortly before initial admixture
	msprime.MassMigration(time = T_PULSE1, source = 3, destination = 0, proportion = m_PULSE1 ), # Neand1 to EUR_EAS pulse of introgression
	msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
        msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
        msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
	msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
	msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
        msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
	msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_Ch
	msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]
    two_pulse = [
        msprime.MassMigration(time = T_PULSE2, source = 4, destination = 1, proportion = m_PULSE2), # Neand2 to EAS pulse of introgression
	msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (2, 3)), # migration between "B" and Africa begins
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)),
        msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
	#msprime.PopulationParametersChange(time = T_B_BN, initial_size = N_B_BN, growth_rate = 0, population_id = 3), # population bottleneck of B begins
	#msprime.PopulationParametersChange(time = T_B_BN + 20, initial_size = N_B, growth_rate = 0, population_id = 3), # population bottleneck of B ends shortly before initial admixture
	msprime.MassMigration(time = T_PULSE1, source = 3, destination = 0, proportion = m_PULSE1), # Neand1 to EUR_EAS pulse of introgression
	msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
        msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
        msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
	msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
	msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
        msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
	msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_CH
	msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]
    
   #    dp = msprime.DemographyDebugger(
#    		    Ne = N_A,
#    		    population_configurations = [
#					        msprime.PopulationConfiguration(
#					        initial_size = N_N1, sample_size=1), # population 0
#					        msprime.PopulationConfiguration(
#					        initial_size = N_N2, sample_size=1), # population 1
#					        msprime.PopulationConfiguration(
#					        initial_size = N_AF, sample_size=1), # population 2
#					        msprime.PopulationConfiguration(
#					        initial_size = N_EU, growth_rate = r_EU, sample_size=1), # population 3
#					        msprime.PopulationConfiguration(
#					        initial_size = N_AS, growth_rate = r_AS, sample_size=1), # population 4
#						msprime.PopulationConfiguration(
#						initial_size = N_CH, sample_size=1),
#						msprime.PopulationConfiguration(
#						initial_size = N_DE, sample_size=1)
#					    ],
#    		    migration_matrix = migration_matrix,
#    		    demographic_events = demographic_events
#    		    )
    #dp.print_history()

    if (pulses == 1):
        demographic_events = one_pulse
    elif (pulses == 2):
        demographic_events = two_pulse


    if (haplo == "F4Dstat"):
####### HAPLOTYPE SIMULATION FOR F4 and DSTAT CALCULATIONS ############
	    replicates = msprime.simulate(
			    Ne = N_A,
			    length = length,
			    recombination_rate = recombination_rate,
			    mutation_rate = mutation_rate,
			    samples = samples,
			    population_configurations = population_configurations,
			    migration_matrix = migration_matrix,
			    demographic_events = demographic_events,
			    num_replicates = 20,
			    random_seed = seed
			    )
	    geno_outfile = open(outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')  ##  sim.eigenstratgeno.n1_0.01_n2_0.05_200
	    snp_outfile = open(outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	    ind_outfile = open(outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	
	    rs_num = 0
	
	
	    for t, tree_sequence in enumerate(replicates): # enumerate creates a list of tuples, where each chrom is a tuple with an index and the tree e.g. [(1,Tree1), (2,Tree2)...]
	        ###  WRITE THE .IND FILE  ###
	        if t<1:                                                        ## Write .ind file based on output of Tree1
	                for i in range(tree_sequence.get_sample_size()):        ## Get the sample size from the tree
	                        if tree_sequence.get_population(i)==0:          ## For each individual in Tree1, get the corresponding population
	                                pop='Neand1'
	                        elif tree_sequence.get_population(i)==1:
	                                pop='Neand2'
	                        elif tree_sequence.get_population(i)==2:
	                                pop='AFR'
	                        elif tree_sequence.get_population(i)==3:
	                                pop='EUR'
	                        elif tree_sequence.get_population(i)==4:
	                                pop='ASN'
	                        elif tree_sequence.get_population(i)==5:
	                                pop='Chimp'
	                        elif tree_sequence.get_population(i)==6:
	                                pop='Deni'
	                        ind_entry = 'Sample_'+str(i)+'\t'+'U'+'\t'+pop  ## Write, as a string, the sample entry using the population read from Tree1
	                        ind_outfile.write(ind_entry+'\n')
		###  WRITE THE .EIGENSTRATGENO AND .SNP FILES  ###
		chr_num = t+1
	        for variant in tree_sequence.variants(as_bytes=True):
		    rs_num+=1
		    geno_outfile.write(variant.genotypes+'\n') # write genotypes to .eigenstratgeno file
		    line = str('rs'+str(rs_num)+'\t'+str(chr_num)+'\t'+str(variant.position/length)+'\t'+str(int(variant.position))+'\t'+'A'+'\t'+'T'+'\n')
		    snp_outfile.write(line) # write snp_allele info to .snp file
	
	
	    parF4_outfile = open('parfile.F4stat.'+outdir+'.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	    parF4_outfile.write('genotypename: '+outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parF4_outfile.write('snpname: '+outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parF4_outfile.write('indivname: '+outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parF4_outfile.write('popfilename: sim.popfile_F4stat'+'\n')
	  
	    parD_outfile = open('parfile.Dstat.'+outdir+'.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	    parD_outfile.write('genotypename: '+outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parD_outfile.write('snpname: '+outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parD_outfile.write('indivname: '+outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parD_outfile.write('popfilename: sim.popfile_Dstat'+'\n')
	
	
	    geno_outfile.close()
	    snp_outfile.close()
	    ind_outfile.close()
	    parF4_outfile.close()
	    parD_outfile.close()


    elif(haplo == "haplo"):
####### HAPLOTYPE SIMULATION FOR DESERT DISTRIBUTION ############
	    return msprime.simulate(
			    Ne = N_A,
			    length = length,	# Simulate 10Mb chromosomes for looking at desert distribution
			    recombination_rate = recombination_rate,
			    mutation_rate = mutation_rate,
			    samples = samples,
			    population_configurations = population_configurations,
			    migration_matrix = migration_matrix,
			    demographic_events = demographic_events,
			    num_replicates = 1,
			    random_seed = seed
			    )
	
################################
################################

def Sriram_demography(S_N1, S_N2, S_AF, S_EU, S_AS, pulses, seed, n1_admix_prop, n2_admix_prop, outdir, t_n1_n2, haplo, length):
    N_A = 10000 # effective population size used for scaling
    N_N1 = 10000 # Altai/Vindija lineage effective population size
    N_N2 = 10000 # Eastern Neandertal effective population size
    N_N_BN = 100 # Neandertal population bottleneck size
    N_B = 10000 # Population size after EUR + ASN join
    N_B_BN = 100 # Population size of EUR_ASN bottleneck (following initial admixture event)
    N_AF = 10000
    N_EU0 = 1000
    N_AS0 = 510
    N_CH = 10000
    N_DE = 10000
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_MH_CH = 7e6 / generation_time # Chimp_Human lineage split (7 MYA)
    T_MH_N = 700e3 / generation_time # Neanderthal / modern human split (700 KYA)
    T_DE_N = 500e3 / generation_time # Denisovan / Neandertal split time
    T_N1_N2 = ((t_n1_n2) * 1e3) / generation_time # Altai/Vindija and Eastern Neandertal popualation split, default set to 350kya
    #T_N1_N2 = 175e3 / generation_time # Altai/Vindija and Eastern Neandertal popualation split	### DIFFERENT FROM TENNESSEN
    #T_AF = 220e3 / generation_time # African population expansion
    T_N_BN = 150e3 / generation_time # Neandertal bottleneck			### DIFFERENT FROM TENNESSEN
    T_B = 75e3 / generation_time # out of Africa				### DIFFERENT FROM TENNESSEN
    T_PULSE1 = 50e3 / generation_time # Neandertal introgression pulse 1	### DIFFERENT FROM TENNESSEN
    T_B_BN = 45e3 / generation_time # A 3rd Bottleneck following shortly after initial admixture event	### SET AT 1800 GEN
    T_EU_AS = 23e3 / generation_time # European / East Asian split
    T_PULSE2 = 18e3 / generation_time # Neandertal introgression pulse 2
    ## We need to work out the starting (diploid) population sizes based on
    ## the growth rates provided for these two populations
    #r_EU = 0.004
    #r_AS = 0.0055
    #N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    #N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    N_EU = 10000
    N_AS = 10000
    
    ## Migration rates during the various epochs.
    #m_AF_B = 25e-5
    #m_AF_EU = 3e-5
    #m_AF_AS = 1.9e-5
    #m_EU_AS = 9.6e-5
    ####
    m_AF_B = 0
    m_AF_EU = 0
    m_AF_AS = 0
    m_EU_AS = 0
    ####
    m_N1_N2 = 0
    m_N1_AF = 0
    m_N1_EU = 0
    m_N1_AS = 0
    m_N2_AF = 0
    m_N2_EU = 0
    m_N2_AS = 0
    m_N1_CH = 0
    m_N2_CH = 0
    m_AF_CH = 0
    m_EU_CH = 0
    m_AS_CH = 0
    m_DE_N1 = 0
    m_DE_N2 = 0
    m_DE_AF = 0
    m_DE_EU = 0
    m_DE_AS = 0
    m_DE_CH = 0
    # Set Percent introgression level
    m_PULSE1 = n1_admix_prop # an instantaneous movement of population from EU_AS to Neand1 (backward in time)
    m_PULSE2 = n2_admix_prop
    #print(m_PULSE1)
    
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
        initial_size = N_N1), # population 0
        msprime.PopulationConfiguration(
        initial_size = N_N2), # population 1
        msprime.PopulationConfiguration(
        initial_size = N_AF), # population 2
        msprime.PopulationConfiguration(
        initial_size = N_EU), # population 3
        msprime.PopulationConfiguration(
        initial_size = N_AS), # population 4
	msprime.PopulationConfiguration(
	initial_size = N_CH), #popultaion 5
	msprime.PopulationConfiguration(
	initial_size = N_DE)  #population 6
    ]
    # specify desired sample
    N1_sample = [msprime.Sample(population = 0, time = 100e3 / generation_time)] * S_N1
    N2_sample = [msprime.Sample(population = 1, time = 100e3 / generation_time)] * S_N2
    AF_sample = [msprime.Sample(population = 2, time = 0)] * S_AF
    EU_sample = [msprime.Sample(population = 3, time = 0)] * S_EU
    AS_sample = [msprime.Sample(population = 4, time = 0)] * S_AS
    CH_sample = [msprime.Sample(population = 5, time = 0)] * 1
    DE_sample = [msprime.Sample(population = 6, time = 100e3/ generation_time)] * 1
    samples = N1_sample + N2_sample + AF_sample + EU_sample + AS_sample + CH_sample + DE_sample
    # set up the initial migration matrix
    migration_matrix = [
        [      0, m_N1_N2, m_N1_AF, m_N1_EU, m_N1_AS, m_N1_CH, m_DE_N1],
        [m_N1_N2,       0, m_N2_AF, m_N2_EU, m_N2_AS, m_N2_CH, m_DE_N2],
        [m_N1_AF, m_N2_AF,       0, m_AF_EU, m_AF_AS, m_AF_CH, m_DE_AF],
        [m_N1_EU, m_N2_EU, m_AF_EU,       0, m_EU_AS, m_EU_CH, m_DE_EU],
        [m_N1_AS, m_N2_AS, m_AF_AS, m_EU_AS,       0, m_AS_CH, m_DE_AS],
	[m_N1_CH, m_N2_CH, m_AF_CH, m_EU_CH, m_AS_CH,       0, m_DE_CH],
	[m_DE_N1, m_DE_N2, m_DE_AF, m_DE_EU, m_DE_AS, m_DE_CH,       0]
    ]
    one_pulse = [
        msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
	msprime.PopulationParametersChange(time = T_B_BN, initial_size = N_B_BN, growth_rate = 0, population_id = 3),  # population bottleneck of B begins
	msprime.PopulationParametersChange(time = T_B_BN + 20, initial_size = N_B, growth_rate = 0, population_id = 3), # population bottleneck ends shortly before initial admixture
	msprime.MassMigration(time = T_PULSE1, source = 3, destination = 0, proportion = m_PULSE1 ), # Neand1 to EUR_EAS pulse of introgression
	msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
        msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        msprime.PopulationParametersChange(time = T_N_BN, initial_size = N_N_BN, growth_rate = 0, population_id = 0), # Neand1 has short bottleneck 
	msprime.PopulationParametersChange(time = T_N_BN + 20, initial_size = N_N1, growth_rate = 0, population_id = 0), # Neand1 bottleneck ends
	msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
	msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
	msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
        msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
	msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_Ch
	msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]
    two_pulse = [
        msprime.MassMigration(time = T_PULSE2, source = 4, destination = 1, proportion = m_PULSE2), # Neand2 to EAS pulse of introgression
	msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
       	msprime.PopulationParametersChange(time = T_B_BN, initial_size = N_B_BN, growth_rate = 0, population_id = 3), # population bottleneck of B begins
	msprime.PopulationParametersChange(time = T_B_BN + 20, initial_size = N_B, growth_rate = 0, population_id = 3), # population bottleneck of B ends shortly before initial admixture
	msprime.MassMigration(time = T_PULSE1, source = 3, destination = 0, proportion = m_PULSE1), # Neand1 to EUR_EAS pulse of introgression
	msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
        msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        msprime.PopulationParametersChange(time = T_N_BN, initial_size = N_N_BN, growth_rate = 0, population_id = 0), # Neand1 has short bottleneck 
	msprime.PopulationParametersChange(time = T_N_BN + 20, initial_size = N_A, growth_rate = 0, population_id = 0), # Neand1 bottleneck ends
        msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_A, population_id = 0), # set parameters of ancestral Neandertal population
	msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
	msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_A, population_id = 0),
        msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
	msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_CH
	msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]


    if (pulses == 1):
        demographic_events = one_pulse
    elif (pulses == 2):
        demographic_events = two_pulse


    if (haplo == "F4Dstat"):
####### HAPLOTYPE SIMULATION FOR F4 and DSTAT CALCULATIONS ############
	    replicates = msprime.simulate(
			    Ne = N_A,
			    length = length,
			    recombination_rate = recombination_rate,
			    mutation_rate = mutation_rate,
			    samples = samples,
			    population_configurations = population_configurations,
			    migration_matrix = migration_matrix,
			    demographic_events = demographic_events,
			    num_replicates = 20,
			    random_seed = seed
			    )
	    geno_outfile = open(outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')  ##  sim.eigenstratgeno.n1_0.01_n2_0.05_200
	    snp_outfile = open(outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	    ind_outfile = open(outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	
	    rs_num = 0
	
	
	    for t, tree_sequence in enumerate(replicates): # enumerate creates a list of tuples, where each chrom is a tuple with an index and the tree e.g. [(1,Tree1), (2,Tree2)...]
	        ###  WRITE THE .IND FILE  ###
	        if t<1:                                                        ## Write .ind file based on output of Tree1
	                for i in range(tree_sequence.get_sample_size()):        ## Get the sample size from the tree
	                        if tree_sequence.get_population(i)==0:          ## For each individual in Tree1, get the corresponding population
	                                pop='Neand1'
	                        elif tree_sequence.get_population(i)==1:
	                                pop='Neand2'
	                        elif tree_sequence.get_population(i)==2:
	                                pop='AFR'
	                        elif tree_sequence.get_population(i)==3:
	                                pop='EUR'
	                        elif tree_sequence.get_population(i)==4:
	                                pop='ASN'
	                        elif tree_sequence.get_population(i)==5:
	                                pop='Chimp'
	                        elif tree_sequence.get_population(i)==6:
	                                pop='Deni'
	                        ind_entry = 'Sample_'+str(i)+'\t'+'U'+'\t'+pop  ## Write, as a string, the sample entry using the population read from Tree1
	                        ind_outfile.write(ind_entry+'\n')
		###  WRITE THE .EIGENSTRATGENO AND .SNP FILES  ###
		chr_num = t+1
	        for variant in tree_sequence.variants(as_bytes=True):
		    rs_num+=1
		    geno_outfile.write(variant.genotypes+'\n') # write genotypes to .eigenstratgeno file
		    line = str('rs'+str(rs_num)+'\t'+str(chr_num)+'\t'+str(variant.position/length)+'\t'+str(int(variant.position))+'\t'+'A'+'\t'+'T'+'\n')
		    snp_outfile.write(line) # write snp_allele info to .snp file
	
	
	    parF4_outfile = open('parfile.F4stat.'+outdir+'.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	    parF4_outfile.write('genotypename: '+outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parF4_outfile.write('snpname: '+outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parF4_outfile.write('indivname: '+outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parF4_outfile.write('popfilename: sim.popfile_F4stat'+'\n')
	  
	    parD_outfile = open('parfile.Dstat.'+outdir+'.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	    parD_outfile.write('genotypename: '+outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parD_outfile.write('snpname: '+outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parD_outfile.write('indivname: '+outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parD_outfile.write('popfilename: sim.popfile_Dstat'+'\n')
	
	
	    geno_outfile.close()
	    snp_outfile.close()
	    ind_outfile.close()
	    parF4_outfile.close()
	    parD_outfile.close()


    elif(haplo == "haplo"):
####### HAPLOTYPE SIMULATION FOR DESERT DISTRIBUTION ############
	    return msprime.simulate(
			    Ne = N_A,
			    length = length,	# Simulate 10Mb chromosomes for looking at desert distribution
			    recombination_rate = recombination_rate,
			    mutation_rate = mutation_rate,
			    samples = samples,
			    population_configurations = population_configurations,
			    migration_matrix = migration_matrix,
			    demographic_events = demographic_events,
			    num_replicates = 1,
			    random_seed = seed
			    )
	
################################
################################

def SplitPop_demography(S_N1, S_N2, S_AF, S_EU, S_AS, pulses, seed, n1_admix_prop, n2_admix_prop, outdir, t_n1_n2, haplo, length):
    N_A = 7310 # effective population size used for scaling
    N_N1 = 1000 # Altai/Vindija lineage effective population size
    N_N2 = 1000 # Eastern Neandertal effective population size
    N_B = 1861 # Population size after EUR + ASN join
    #N_B_BN = 100 # Population size of EUR_ASN bottleneck (following initial admixture event)
    N_AF = 14474
    N_EU0 = 1032
    N_AS0 = 554
    N_CH = 1000
    N_DE = 1000
    N_SP = 100
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_MH_CH = 7e6 / generation_time # Chimp_Human lineage split
    T_MH_N = 700e3 / generation_time # Neanderthal / modern human split
    T_DE_N = 500e3 / generation_time # Denisovan / Neandertal split time
    T_N1_N2 = ((t_n1_n2) * 1e3) / generation_time # Altai/Vindija and Eastern Neandertal popualation split, default set to 350kya
    T_AF = 148e3 / generation_time # African population expansion
    T_B = 100e3 / generation_time # out of Africa
    T_PULSE1 = 55e3 / generation_time # Neandertal introgression pulse 1
    T_SP = 54e3 /generation_time # split population off of B
    #T_B_BN = 45e3 / generation_time # A 3rd Bottleneck following shortly after initial admixture event
    T_EU_AS = 23e3 / generation_time # European / East Asian split
    T_PULSE2 = 18e3 / generation_time # Neandertal introgression pulse 2
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 15e-5
    m_AF_EU = 2.5e-5
    m_AF_AS = 0.78e-5
    m_EU_AS = 3.11e-5
    m_N1_N2 = 0
    m_N1_AF = 0
    m_N1_EU = 0
    m_N1_AS = 0
    m_N2_AF = 0
    m_N2_EU = 0
    m_N2_AS = 0
    m_N1_CH = 0
    m_N2_CH = 0
    m_AF_CH = 0
    m_EU_CH = 0
    m_AS_CH = 0
    m_DE_N1 = 0
    m_DE_N2 = 0
    m_DE_AF = 0
    m_DE_EU = 0
    m_DE_AS = 0
    m_DE_CH = 0
    m_SP_N1 = 0
    m_SP_N2 = 0
    m_SP_AF = 0
    m_SP_EU = 0
    m_SP_AS = 0
    m_SP_CH = 0
    m_SP_DE = 0

    # Set Percent introgression level
    m_PULSE1 = n1_admix_prop # an instantaneous movement of population from EU_AS to Neand1 (backward in time)
    m_PULSE2 = n2_admix_prop
    #print(m_PULSE1)
    
    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
        initial_size = N_N1), # population 0
        msprime.PopulationConfiguration(
        initial_size = N_N2), # population 1
        msprime.PopulationConfiguration(
        initial_size = N_AF), # population 2
        msprime.PopulationConfiguration(
        initial_size = N_EU, growth_rate = r_EU), # population 3
        msprime.PopulationConfiguration(
        initial_size = N_AS, growth_rate = r_AS), # population 4
	msprime.PopulationConfiguration(
	initial_size = N_CH), # popultaion 5
	msprime.PopulationConfiguration(
	initial_size = N_DE),  # population 6
	msprime.PopulationConfiguration(
	initial_size = N_SP)	# population 7
    ]
    # specify desired sample
    N1_sample = [msprime.Sample(population = 0, time = 100e3 / generation_time)] * S_N1
    N2_sample = [msprime.Sample(population = 1, time = 100e3 / generation_time)] * S_N2
    AF_sample = [msprime.Sample(population = 2, time = 0)] * S_AF
    EU_sample = [msprime.Sample(population = 3, time = 0)] * S_EU
    AS_sample = [msprime.Sample(population = 4, time = 0)] * S_AS
    CH_sample = [msprime.Sample(population = 5, time = 0)] * 1
    DE_sample = [msprime.Sample(population = 6, time = 100e3/ generation_time)] * 1
    samples = N1_sample + N2_sample + AF_sample + EU_sample + AS_sample + CH_sample + DE_sample
    # set up the initial migration matrix
    migration_matrix = [
        [      0, m_N1_N2, m_N1_AF, m_N1_EU, m_N1_AS, m_N1_CH, m_DE_N1, m_SP_N1],
        [m_N1_N2,       0, m_N2_AF, m_N2_EU, m_N2_AS, m_N2_CH, m_DE_N2, m_SP_N2],
        [m_N1_AF, m_N2_AF,       0, m_AF_EU, m_AF_AS, m_AF_CH, m_DE_AF, m_SP_AF],
        [m_N1_EU, m_N2_EU, m_AF_EU,       0, m_EU_AS, m_EU_CH, m_DE_EU, m_SP_EU],
        [m_N1_AS, m_N2_AS, m_AF_AS, m_EU_AS,       0, m_AS_CH, m_DE_AS, m_SP_AS],
	[m_N1_CH, m_N2_CH, m_AF_CH, m_EU_CH, m_AS_CH,       0, m_DE_CH, m_SP_CH],
	[m_DE_N1, m_DE_N2, m_DE_AF, m_DE_EU, m_DE_AS, m_DE_CH,       0, m_SP_DE],
	[m_SP_N1, m_SP_N2, m_SP_AF, m_SP_EU, m_SP_AS, m_SP_CH, m_SP_DE,       0]
    ]
    one_pulse = [
        msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (2, 3)), # migration between "B" and Africa begins
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)),
        msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
	msprime.MassMigration(time = T_SP, source = 3, destination = 7, proportion = 0.1),  # Population splits off of EUR_EAS, 10% migration rate to SplitPop
	msprime.PopulationParametersChange(time = T_SP, initial_size = N_SP, growth_rate = 0, population_id = 7), # SplitPop is set at size Ne=100
	msprime.MassMigration(time = T_PULSE1, source = 7, destination = 0, proportion = m_PULSE1 ), # Neand1 to SplitPop pulse of introgression
	msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into African at T_B
        msprime.MassMigration(time = T_B, source = 7, destination = 2, proportion = 1.0), # SplitPop merges into African at T_B
	msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
        msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
	msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
	msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
        msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
	msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_Ch
	msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]

    two_pulse = [
        msprime.MassMigration(time = T_PULSE2, source = 4, destination = 1, proportion = m_PULSE2), # Neand2 to EAS pulse of introgression
	msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (2, 3)), # migration between "B" and Africa begins
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)),
        msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
	msprime.MassMigration(time = T_SP, source = 3, destination = 7, proportion = 0.1),  # Population splits off of EUR_EAS, 10% migration rate to SplitPop
	msprime.PopulationParametersChange(time = T_SP, initial_size = N_SP, growth_rate = 0, population_id = 7), # SplitPop is set at size Ne=100
	msprime.MassMigration(time = T_PULSE1, source = 7, destination = 0, proportion = m_PULSE1 ), # Neand1 to SplitPop pulse of introgression
	msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into African at T_B
        msprime.MassMigration(time = T_B, source = 7, destination = 2, proportion = 1.0), # SplitPop merges into African at T_B
        msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
        msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
	msprime.MassMigration(time = T_DE_N, source = 6, destination = 0, proportion = 1.0), # DE merges with N1
	msprime.PopulationParametersChange(time = T_DE_N, initial_size = N_N1, population_id = 0),
        msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2), # set parameters of ancetral hominin population
	msprime.MassMigration(time = T_MH_CH, source = 5, destination = 2, proportion = 1.0), # Chimp lineage merges into ancestral hominin population at time T_MH_CH
	msprime.PopulationParametersChange(time = T_MH_CH, initial_size = N_A, population_id = 2) # set parameters of ancestral hominin population
    ]
    
   #    dp = msprime.DemographyDebugger(
#    		    Ne = N_A,
#    		    population_configurations = [
#					        msprime.PopulationConfiguration(
#					        initial_size = N_N1, sample_size=1), # population 0
#					        msprime.PopulationConfiguration(
#					        initial_size = N_N2, sample_size=1), # population 1
#					        msprime.PopulationConfiguration(
#					        initial_size = N_AF, sample_size=1), # population 2
#					        msprime.PopulationConfiguration(
#					        initial_size = N_EU, growth_rate = r_EU, sample_size=1), # population 3
#					        msprime.PopulationConfiguration(
#					        initial_size = N_AS, growth_rate = r_AS, sample_size=1), # population 4
#						msprime.PopulationConfiguration(
#						initial_size = N_CH, sample_size=1),
#						msprime.PopulationConfiguration(
#						initial_size = N_DE, sample_size=1)
#					    ],
#    		    migration_matrix = migration_matrix,
#    		    demographic_events = demographic_events
#    		    )
    #dp.print_history()

    if (pulses == 1):
        demographic_events = one_pulse
    elif (pulses == 2):
        demographic_events = two_pulse

    if (haplo == "F4Dstat"):
####### HAPLOTYPE SIMULATION FOR F4 and DSTAT CALCULATIONS ############
	    replicates = msprime.simulate(
			    Ne = N_A,
			    length = length,
			    recombination_rate = recombination_rate,
			    mutation_rate = mutation_rate,
			    samples = samples,
			    population_configurations = population_configurations,
			    migration_matrix = migration_matrix,
			    demographic_events = demographic_events,
			    num_replicates = 20,
			    random_seed = seed
			    )
	    geno_outfile = open(outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')  ##  sim.eigenstratgeno.n1_0.01_n2_0.05_200
	    snp_outfile = open(outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	    ind_outfile = open(outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	
	    rs_num = 0
	
	
	    for t, tree_sequence in enumerate(replicates): # enumerate creates a list of tuples, where each chrom is a tuple with an index and the tree e.g. [(1,Tree1), (2,Tree2)...]
	        ###  WRITE THE .IND FILE  ###
	        if t<1:                                                        ## Write .ind file based on output of Tree1
	                for i in range(tree_sequence.get_sample_size()):        ## Get the sample size from the tree
	                        if tree_sequence.get_population(i)==0:          ## For each individual in Tree1, get the corresponding population
	                                pop='Neand1'
	                        elif tree_sequence.get_population(i)==1:
	                                pop='Neand2'
	                        elif tree_sequence.get_population(i)==2:
	                                pop='AFR'
	                        elif tree_sequence.get_population(i)==3:
	                                pop='EUR'
	                        elif tree_sequence.get_population(i)==4:
	                                pop='ASN'
	                        elif tree_sequence.get_population(i)==5:
	                                pop='Chimp'
	                        elif tree_sequence.get_population(i)==6:
	                                pop='Deni'
	                        ind_entry = 'Sample_'+str(i)+'\t'+'U'+'\t'+pop  ## Write, as a string, the sample entry using the population read from Tree1
	                        ind_outfile.write(ind_entry+'\n')
		###  WRITE THE .EIGENSTRATGENO AND .SNP FILES  ###
		chr_num = t+1
	        for variant in tree_sequence.variants(as_bytes=True):
		    rs_num+=1
		    geno_outfile.write(variant.genotypes+'\n') # write genotypes to .eigenstratgeno file
		    line = str('rs'+str(rs_num)+'\t'+str(chr_num)+'\t'+str(variant.position/length)+'\t'+str(int(variant.position))+'\t'+'A'+'\t'+'T'+'\n')
		    snp_outfile.write(line) # write snp_allele info to .snp file
	
	
	    parF4_outfile = open('parfile.F4stat.'+outdir+'.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	    parF4_outfile.write('genotypename: '+outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parF4_outfile.write('snpname: '+outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parF4_outfile.write('indivname: '+outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parF4_outfile.write('popfilename: sim.popfile_F4stat'+'\n')
	  
	    parD_outfile = open('parfile.Dstat.'+outdir+'.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed), 'w')
	    parD_outfile.write('genotypename: '+outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parD_outfile.write('snpname: '+outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parD_outfile.write('indivname: '+outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_t_'+str(t_n1_n2)+'_'+str(seed)+'\n')
	    parD_outfile.write('popfilename: sim.popfile_Dstat'+'\n')
	
	
	    geno_outfile.close()
	    snp_outfile.close()
	    ind_outfile.close()
	    parF4_outfile.close()
	    parD_outfile.close()


    elif(haplo == "haplo"):
####### HAPLOTYPE SIMULATION FOR DESERT DISTRIBUTION ############
	    return msprime.simulate(
			    Ne = N_A,
			    length = length,	# Simulate 10Mb chromosomes for looking at desert distribution
			    recombination_rate = recombination_rate,
			    mutation_rate = mutation_rate,
			    samples = samples,
			    population_configurations = population_configurations,
			    migration_matrix = migration_matrix,
			    demographic_events = demographic_events,
			    num_replicates = 1,
			    random_seed = seed
			    )

	    
	    
#######################################
#######################################
