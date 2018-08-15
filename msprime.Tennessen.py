from __future__ import print_function
import msprime
import math
import numpy as np
import pybedtools
import random
import scipy
import subprocess
import glob
from scipy import spatial
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-p", "--population", action = "store", type = "string", dest = "pop")
parser.add_option("-o", "--outdir", action = "store", type = "string", dest = "outdir")
parser.add_option("-s", "--seed", action = "store", type = "int", dest = "seed")
parser.add_option("-i", "--introgress_pulses", action = "store", type = "int", dest = "pulses")
parser.add_option("-n", "--neand1_admixture_proportion", action = "store", type = "float", dest = "n1_admix_prop")
parser.add_option("-d", "--deni_admixture_propotion", action = "store", type = "float", dest = "n2_admix_prop")
#parser.add_option("-t", "--taskid", action = "store", type = "string", dest = "taskid")
(options, args) = parser.parse_args()

def simulate_demography(S_N1, S_N2, S_AF, S_EU, S_AS, pulses, seed, n1_admix_prop, n2_admix_prop, outdir):
    N_A = 7310 # effective population size used for scaling
    N_N1 = 1000 # Altai/Vindija lineage effective population size
    N_N2 = 1000 # Eastern Neandertal effective population size
    N_B = 2100 # Population size after EUR + ASN join
    N_B_BN = 100 # Population size of EUR_ASN bottleneck (following initial admixture event)
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    N_CH = 1000
    N_DE = 1000
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_MH_CH = 7e6 / generation_time # Chimp_Human lineage split
    T_MH_N = 700e3 / generation_time # Neanderthal / modern human split
    T_DE_N = 500e3 / generation_time # Denisovan / Neandertal split time
    T_N1_N2 = 350e3 / generation_time # Altai/Vindija and Eastern Neandertal popualation split
    T_AF = 220e3 / generation_time # African population expansion
    T_B = 140e3 / generation_time # out of Africa
    T_PULSE1 = 55e3 / generation_time # Neandertal introgression pulse 1
    T_B_BN = 45e3 / generation_time # A 3rd Bottleneck following shortly after initial admixture event
    T_EU_AS = 21.2e3 / generation_time # European / East Asian split
    T_PULSE2 = 18e3 / generation_time # Neandertal introgression pulse 2
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU = 0.004
    r_AS = 0.0055
    N_EU = N_EU0 / math.exp(-r_EU * T_EU_AS)
    N_AS = N_AS0 / math.exp(-r_AS * T_EU_AS)
    # Migration rates during the various epochs.
    m_AF_B = 25e-5
    m_AF_EU = 3e-5
    m_AF_AS = 1.9e-5
    m_EU_AS = 9.6e-5
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
    print(m_PULSE1)
    
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
	msprime.PopulationParametersChange(time = T_B_BN, initial_size = N_B_BN, growth_rate = 0, population_id = 3),  # population bottleneck of B begins
	msprime.PopulationParametersChange(time = T_B_BN + 20, initial_size = N_B, growth_rate = 0, population_id = 3), # population bottleneck ends shortly before initial admixture
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
       	msprime.PopulationParametersChange(time = T_B_BN, initial_size = N_B_BN, growth_rate = 0, population_id = 3), # population bottleneck of B begins
	msprime.PopulationParametersChange(time = T_B_BN + 20, initial_size = N_B, growth_rate = 0, population_id = 3), # population bottleneck of B ends shortly before initial admixture
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


####### HAPLOTYPE SIMULATION FOR F4 and DSTAT CALCULATIONS ############

    length = 1e6	# simulate 1Mb chromosomes for faster F4 and D calculation
    recombination_rate = 2e-8
    mutation_rate = 1e-8 # see Shendure and Akey, 2015
#
#
#    replicates = msprime.simulate(
#		    Ne = N_A,
#		    length = length,
#		    recombination_rate = recombination_rate,
#		    mutation_rate = mutation_rate,
#		    samples = samples,
#		    population_configurations = population_configurations,
#		    migration_matrix = migration_matrix,
#		    demographic_events = demographic_events,
#		    num_replicates = 20,
#		    random_seed = seed
#		    )
#    geno_outfile = open(outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_'+str(seed), 'w')  ##  sim.eigenstratgeno.n1_0.01_n2_0.05_200
#    snp_outfile = open(outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_'+str(seed), 'w')
#    ind_outfile = open(outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_'+str(seed), 'w')
#
#    rs_num = 0
#
#
#    for t, tree_sequence in enumerate(replicates): # enumerate creates a list of tuples, where each chrom is a tuple with an index and the tree e.g. [(1,Tree1), (2,Tree2)...]
#        ###  WRITE THE .IND FILE  ###
#        if t<1:                                                        ## Write .ind file based on output of Tree1
#                for i in range(tree_sequence.get_sample_size()):        ## Get the sample size from the tree
#                        if tree_sequence.get_population(i)==0:          ## For each individual in Tree1, get the corresponding population
#                                pop='Neand1'
#                        elif tree_sequence.get_population(i)==1:
#                                pop='Neand2'
#                        elif tree_sequence.get_population(i)==2:
#                                pop='AFR'
#                        elif tree_sequence.get_population(i)==3:
#                                pop='EUR'
#                        elif tree_sequence.get_population(i)==4:
#                                pop='ASN'
#                        elif tree_sequence.get_population(i)==5:
#                                pop='Chimp'
#                        elif tree_sequence.get_population(i)==6:
#                                pop='Deni'
#                        ind_entry = 'Sample_'+str(i)+'\t'+'U'+'\t'+pop  ## Write, as a string, the sample entry using the population read from Tree1
#                        ind_outfile.write(ind_entry+'\n')
#	###  WRITE THE .EIGENSTRATGENO AND .SNP FILES  ###
#	chr_num = t+1
#        for variant in tree_sequence.variants():
#	    rs_num+=1
#	    geno_outfile.write(variant.genotypes+'\n') # write genotypes to .eigenstratgeno file
#	    line = str('rs'+str(rs_num)+'\t'+str(chr_num)+'\t'+str(variant.position/length)+'\t'+str(int(variant.position))+'\t'+'A'+'\t'+'T'+'\n')
#	    snp_outfile.write(line) # write snp_allele info to .snp file
#
#
#    parF4_outfile = open('parfile.F4stat.'+outdir+'.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_'+str(seed), 'w')
#    parF4_outfile.write('genotypename: '+outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_'+str(seed)+'\n')
#    parF4_outfile.write('snpname: '+outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_'+str(seed)+'\n')
#    parF4_outfile.write('indivname: '+outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_'+str(seed)+'\n')
#    parF4_outfile.write('popfilename: sim.popfile_F4stat'+'\n')
#  
#    parD_outfile = open('parfile.Dstat.'+outdir+'.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_'+str(seed), 'w')
#    parD_outfile.write('genotypename: '+outdir+'.eigenstratgeno.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_'+str(seed)+'\n')
#    parD_outfile.write('snpname: '+outdir+'.snp.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_'+str(seed)+'\n')
#    parD_outfile.write('indivname: '+outdir+'.ind.n1_'+str(n1_admix_prop)+'_n2_'+str(n2_admix_prop)+'_'+str(seed)+'\n')
#    parD_outfile.write('popfilename: sim.popfile_Dstat'+'\n')
#
#
#    geno_outfile.close()
#    snp_outfile.close()
#    ind_outfile.close()
#    parF4_outfile.close()
#    parD_outfile.close()

####### HAPLOTYPE SIMULATION FOR DESERT DISTRIBUTION ############

    return msprime.simulate(
		    Ne = N_A,
		    length = 1e7,	# Simulate 10Mb chromosomes for looking at desert distribution
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

# run the simulation
#def simulate_demography(S_N1, S_N2, S_AF, S_EU, S_AS, pulses, seed, n1_admix_prop, n2_admix_prop, outdir)
simulation = simulate_demography(1, 1, 1, 1006, 1008, options.pulses, options.seed, options.n1_admix_prop, options.n2_admix_prop, options.outdir)

#for t, tree in enumerate(simulation):
#	size = tree.get_sample_size()
#	print(size)


##################################

# define the sample indices
N_samples  = range(0, 2)
AF_samples = range(2, 3)
nonAfr_samples = range(3, 2017)
EU_samples = range(3, 1009)
AS_samples = range(1009, 2017)
Chimp_samples = range(2017, 2018)
Deni_samples = range(2018, 2019)

class introgression_object(object):
    def __init__(self, mutations, haplotypes, intervals):
        self.mutations = mutations
        self.haplotypes = haplotypes
        self.intervals = intervals

# this function requires that the archaic samples come first in the simulation
def get_introgressed_intervals(tree_sequence, human_indices, neand_indices, outdir, taskid):
    # store mutation positions in a list
    mutations = []
    for mutation in tree_sequence.mutations():
        mutations.append(mutation.position)
    if (len(mutations) == 0):
	    #quit()
	    return
    
   #test = print(mutations)
   #return test
##  PASSED  ##    

    # store the haplotypes in a list
    haplotypes = []
    for haplotype in tree_sequence.haplotypes():
        haplotypes.append(haplotype)
    if (len(haplotypes) == 0):
	    #quit()
	    return

    #test = print(haplotypes)
    #return test
##  PASSED  ##

    # store introgressed intervals in a dictionary of BedTools lists
    introgressed_haplotypes = {}
    # loop through trees and samples
    for tree in tree_sequence.trees():
        interval = tree.get_interval()
        start = interval[0]
        end = interval[1]
        introgressed_haplotypes[interval] = []
        # for each tree, get the samples that are introgressed
        for leaf in human_indices:
            if tree.get_population(tree.get_mrca(0, leaf)) in neand_indices:
                # if the common ancestor with Neanderthal leaf is in the Neandertal population (i.e. is introgressed)
                introgressed_haplotypes[interval].append(leaf)
   #test = print(introgressed_haplotypes)
   #return test
##  PASSED  ##


    # for each interval, record introgressed samples in an array
    introgression_list = []
    for interval in introgressed_haplotypes:
        n_leaves = len(introgressed_haplotypes[interval])
        if (n_leaves > 0):
            c = np.array([int(1), int(math.ceil(interval[0])), int(math.ceil(interval[1]))])
            bed_init = np.tile(c, (n_leaves, 1))
            bed_array = np.append(bed_init, np.asarray(introgressed_haplotypes[interval])[..., None], axis = 1)
            introgression_list.append(bed_array)
    #test = print(introgression_list)
    #return test
##  PASSED  ##

    if (len(introgression_list) == 0):
	    #quit()
	    return
    interval_array = np.concatenate(introgression_list)
    
    #test = print(interval_array)
    #return test
##  FAILED ; remove 'quit()' replace with 'return None' ; PASSED  ##

   # use bedtools to operate on the introgressed intervals
    prefix = outdir + str(taskid) + "/"
    COMMAND = "mkdir " + prefix
    return_code = subprocess.call(COMMAND, shell = True)
    np.savetxt(prefix + outdir + str(taskid) + ".bed", interval_array, delimiter = "\t", newline = "\n", fmt = "%s")
    COMMAND = '''awk '{print >> "''' + prefix + outdir + str(taskid) + '''."$4".bed"}' ''' + prefix + outdir + str(taskid) + ".bed; rm " + prefix + outdir + str(taskid) + ".bed"
    return_code = subprocess.call(COMMAND, shell = True)
    COMMAND = '''for i in ''' + prefix + outdir + str(taskid) + '''*.bed; do sort -k2,2n $i > $i.sorted; rm $i; done'''
    return_code = subprocess.call(COMMAND, shell = True)
    COMMAND = '''for i in ''' + prefix + outdir + str(taskid) + '''*.sorted; do bedtools merge -i $i -c 4 -o distinct > $i.merged; rm $i; done'''
    return_code = subprocess.call(COMMAND, shell = True)
    # concatenate the individual sample files
    files = prefix + outdir + str(taskid) + "*.merged"
    read_files = glob.glob(files)
    fid = prefix + outdir + str(taskid) + ".bed"
    with open(fid, "wb") as outfile:
        for f in read_files:
            with open(f, "rb") as infile:
                outfile.write(infile.read())
    # then load them back into python
    interval_bed = pybedtools.BedTool(fid)
    return introgression_object(mutations, haplotypes, interval_bed)

if (options.pop == "EAS"):
    human_samples = AS_samples
elif (options.pop == "EUR"):
    human_samples = EU_samples
elif (options.pop == "nonAfr"):
    human_samples = nonAfr_samples


###  FOR EACH SIMULATED CHROMOSOME, PRINT ALL THE INTROGRESSED HAPLOTYPES BELONGING TO THE SPECIFIED NON-AFR POPULATION IN BED FORMAT  ###
for t, tree in enumerate(simulation):
	print(str(tree.get_sample_size())+'_'+str(t))
	#int_results = get_introgressed_intervals(tree, human_samples, N_samples, options.outdir, iteration)
	#def get_introgressed_intervals(tree_sequence, human_indices, neand_indices, outdir, taskid):
	get_introgressed_intervals(tree, human_samples, N_samples, options.outdir+'_'+options.pop+'_'+str(options.seed), '_n1_'+str(options.n1_admix_prop)+'_n2_'+str(options.n2_admix_prop))
	print('DONE')

print('fin')

#	
#	# choose a random introgressed haplotype and get all overlapping haplotypes
#	sample_index = random.sample(range(len(int_results.intervals)), 1)[0]
#	random_interval = pybedtools.BedTool(str(int_results.intervals[sample_index]), from_string = True)
#	olaps = random_interval.intersect(int_results.intervals, f = 1, wb = True)
#	start = random_interval[0].start
#	end = random_interval[0].end
#	length = end - start
#	
#	# find the start and end coordinates of the mutations
#	m_indices = [index for index, item in enumerate(int_results.mutations) if item >= start and item <= end]
#	if (len(m_indices) == 0):
#	    quit()
#	
#	m_start = m_indices[0]
#	m_end = m_indices[-1]
#	
#	# get all overlapping haplotypes
#	int_haps = []
#	for sample in olaps:
#	    s_index = int(sample.fields[7])
#	    int_hap = map(int, int_results.haplotypes[s_index][m_start:m_end + 1])
#	    int_haps.append(int_hap)
#	
#	int_haps = np.matrix(int_haps)
#	n_haps = int_haps.shape[0]
#	if (n_haps < 2):
#	    quit()
#	
#	manhattan = scipy.spatial.distance.pdist(int_haps, metric = "cityblock")
#	pi_mean = np.mean(manhattan)
#	pi_min = np.min(manhattan)
#	pi_max = np.max(manhattan)
#	
#	# get divergence of each haplotype to the Neandertal haplotype
#	neand = map(int, int_results.haplotypes[N_samples[0]][m_start:m_end + 1])
#	ndists = []
#	for row in int_haps[range(n_haps)]:
#	    ndists.append(scipy.spatial.distance.cityblock(row, neand))
#	
#	ndist_mean = np.mean(ndists)
#	ndist_min = np.min(ndists)
#	ndist_max = np.max(ndists)
#	
#	print(length, pi_mean, pi_min, pi_max, ndist_mean, ndist_min, ndist_max, n_haps, options.pop, task, sep = "\t")
