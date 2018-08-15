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
import msprime_demo_models

usage = "usage: python <script> -p nonAfr -o Tenn -s 1 -i 1 -n 0.02 -d 0.02 -t 350 -c F4Dstat"
parser = OptionParser(usage=usage)
parser.add_option("-p", "--population", action = "store", type = "string", dest = "pop", default="nonAfr", help="call introg_haplotypes in EUR, EAS, or all nonAfr; default=nonAfr")
parser.add_option("-o", "--outdir", action = "store", type = "string", dest = "outdir", default="Tenn", help="output directory name, best to include model name; default=Tenn")
parser.add_option("-s", "--seed", action = "store", type = "int", dest = "seed", default=1, help="Set random seed for replicate chromosomes; default=1")
parser.add_option("-i", "--introgress_pulses", action = "store", type = "int", dest = "pulses", default=1, help="Set number of introgression pulses; default=1")
parser.add_option("-n", "--neand1_admixture_proportion", action = "store", type = "float", dest = "n1_admix_prop", default=0.02, help="Set N1 admixture proportion; default=0.02")
parser.add_option("-d", "--deni_admixture_propotion", action = "store", type = "float", dest = "n2_admix_prop", default=0.02, help="Set N2 admixture proportion; default=0.01")
parser.add_option("-t", "--time_N1_N2_split", action = "store", type = "int", dest = "t_n1_n2", default=350, help="Set N1 N2 split time in kya; default=350")
parser.add_option("-c", "--calls_or_stats", action="store", type = "string", dest = "haplo", default="F4Dstat", help="Pick haplotype calls or F4Dstat output; default=F4Dstat")
(options, args) = parser.parse_args()



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
    
    # store the haplotypes in a list
    haplotypes = []
    for haplotype in tree_sequence.haplotypes():
        haplotypes.append(haplotype)
    if (len(haplotypes) == 0):
	    #quit()
	    return

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

    # for each interval, record introgressed samples in an array
    introgression_list = []
    for interval in introgressed_haplotypes:
        n_leaves = len(introgressed_haplotypes[interval])
        if (n_leaves > 0):
            c = np.array([int(1), int(math.ceil(interval[0])), int(math.ceil(interval[1]))])
            bed_init = np.tile(c, (n_leaves, 1))
            bed_array = np.append(bed_init, np.asarray(introgressed_haplotypes[interval])[..., None], axis = 1)
            introgression_list.append(bed_array)

    if (len(introgression_list) == 0):
	    #quit()
	    return
    interval_array = np.concatenate(introgression_list)
    
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




####### SIMULATION RUNS ############

if (options.outdir == "Tenn"):
	#Tenn_demography(S_N1, S_N2, S_AF, S_EU, S_AS, pulses, seed, n1_admix_prop, n2_admix_prop, outdir, t_n1_n2, haplo):
	simulation = msprime_demo_models.Tenn_demography(1, 1, 1, 1006, 1008, options.pulses, options.seed, options.n1_admix_prop, options.n2_admix_prop, options.outdir, options.t_n1_n2, options.haplo)

elif (options.outdir == "Sriram"):
	simulation = msprime_demo_models.Sriram_demography(1, 1, 1, 1006, 1008, options.pulses, options.seed, options.n1_admix_prop, options.n2_admix_prop, options.outdir, options.t_n1_n2, options.haplo)

elif (options.outdir == "SplitPop"):
	simulation = msprime_demo_models.SplitPop_demography(1, 1, 1, 1006, 1008, options.pulses, options.seed, options.n1_admix_prop, options.n2_admix_prop, options.outdir, options.t_n1_n2, options.haplo)


###### GET HAPLOTYPES FROM SIMULATED TREES #######

if (options.haplo == "haplo"):
	# define the sample indices
	N_samples  = range(0, 2)
	AF_samples = range(2, 3)
	nonAfr_samples = range(3, 2017)
	EU_samples = range(3, 1009)
	AS_samples = range(1009, 2017)
	Chimp_samples = range(2017, 2018)
	Deni_samples = range(2018, 2019)
	if (options.pop == "EAS"):
		human_samples = AS_samples
	elif (options.pop == "EUR"):
		human_samples = EU_samples
	elif (options.pop == "nonAfr"):
		human_samples = nonAfr_samples

##  FOR EACH SIMULATED CHROMOSOME, PRINT ALL THE INTROGRESSED HAPLOTYPES BELONGING TO THE SPECIFIED NON-AFR POPULATION IN BED FORMAT  ###
	for t, tree in enumerate(simulation):
		print(str(tree.get_sample_size())+'_'+str(t))
		#def get_introgressed_intervals(tree_sequence, human_indices, neand_indices, outdir, taskid):
		get_introgressed_intervals(tree, human_samples, N_samples, options.outdir+'_'+options.pop+'_'+str(options.seed), '_n1_'+str(options.n1_admix_prop)+'_n2_'+str(options.n2_admix_prop))
		print('DONE')

print('fin')
