from __future__ import print_function
import sys
import msprime
import math
import numpy as np
import pybedtools
import random
import scipy
import subprocess
import glob
import collections
import gzip
from scipy import spatial
from optparse import OptionParser
import msprime_demo_models

usage = "usage: python <script> -p nonAfr -o Tenn -s 1 -i 1 -n 0.02 -d 0.02 -t 350 -c F4Dstat"
parser = OptionParser(usage=usage)
parser.add_option("-p", "--population", action = "store", type = "string", dest = "pop", default="nonAfr", help="call introg_haplotypes in EUR, EAS, or all nonAfr; default=nonAfr")
parser.add_option("-o", "--outdir", action = "store", type = "string", dest = "outdir", default="Tenn", help="output directory name, best to include model name; default=Tenn")
parser.add_option("-s", "--seed", action = "store", type = "int", dest = "seed", default=1, help="Set random seed for replicate chromosomes; default=1")
parser.add_option("-i", "--introgress_pulses", action = "store", type = "int", dest = "pulses", default=2, help="Set number of introgression pulses; default=2")
parser.add_option("-n", "--neand1_admixture_proportion", action = "store", type = "float", dest = "n1_admix_prop", default=0.02, help="Set N1 admixture proportion; default=0.02")
parser.add_option("-d", "--deni_admixture_propotion", action = "store", type = "float", dest = "n2_admix_prop", default=0.0, help="Set N2 admixture proportion; default=0.0")
parser.add_option("-t", "--time_N1_N2_split", action = "store", type = "int", dest = "t_n1_n2", default=350, help="Set N1 N2 split time in kya; default=350")
parser.add_option("-c", "--calls_or_stats", action="store", type = "string", dest = "haplo", default="F4Dstat", help="Pick haplotype calls ('haplo') or 'F4Dstat' output; default=F4Dstat")
parser.add_option("-l", "--length_chrom", action="store", type= "float", dest = "length", default=1e6, help="Define length of simulated chromosome; default=1e6")
(options, args) = parser.parse_args()



class introgression_object(object):
    def __init__(self, mutations, haplotypes, intervals):
        self.mutations = mutations
        self.haplotypes = haplotypes
        self.intervals = intervals

########## No Longer Use This Function ############
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
########################################################################

def introgressed_samples_fn(ts, neanderthal_mrca, neanderthal_samples, segments):
    # Define the samples that carry introgressed segments?
    trees = ts.trees()
    tree = next(trees)
    for left, right in segments:
##        print("\tINTERVAL:", left, right)
        # Skip ahead to the tree that intersects with this segment.
	# NOTE: Regading tree data; The "Records" group consists of four pieces of information: 
	# 1) the left and 2) right coordinates of the coalescing interval, 3) the list of child nodes (modern human samples) and 4) the parent node.
	# Each record returned can be accessed via the attributes 'left', 'right', 'node', 'children', 'time' and 'population'
	# A record represents the assignment of a pair of children 'c' to a parent 'u'. This assignment happens at 't' generations in the past within the population with ID 'd'
        while tree.get_interval()[0] < left:
	    ## Pass through the potential trees until you get to the one where the start_position matches the start of the defined segment
	    tree = next(trees)
	    #print(tree.get_interval()[0], left)
        assert tree.get_interval()[0] == left
        start = None
        last_human_leaves = None
##	print(start, last_human_leaves)
        while tree is not None and tree.get_interval()[1] <= right:
	    #print( set(tree.leaves(neanderthal_mrca)) )
	    #print( neanderthal_samples )
	    human_leaves = set(tree.leaves(neanderthal_mrca)) - set(neanderthal_samples)
	    #print( human_leaves )
	    ## Create a set() of all the human leaves represented in this tree that share an mrca in the Neandertal lineage
	    # (i.e. mrca node is in Neand population)
	    ## neanderthal_mrca is supplied from the node_map dictionary as the parent node for the tree present in the Neandertal population 
	    ## In python set() is an unordered collection of unique values; 
	    # x = [1, 1, 2, 2, 3, 3] --> set(x) --> set([1, 2, 3])
	    ## tree.leaves(n) returns an iterator over all the leaves in this tree underneath the specified node (in this case the neanderthal_mrca)
##          print("\t\tTREE:", tree.get_interval(), human_leaves)
            if start is None:
                last_human_leaves = human_leaves
                start = tree.get_interval()[0]
            elif human_leaves != last_human_leaves:
                end = tree.get_interval()[0]
                yield start, end, last_human_leaves
                start = end
                last_human_leaves = human_leaves
            tree = next(trees, None)
        yield start, right, last_human_leaves


####### SIMULATION RUNS ############

print("Model: ", options.outdir, "Seed: ", options.seed, "Neand_Pulse1: ", options.n1_admix_prop, "Neand_Pulse2: ", options.n2_admix_prop, "Length: ", options.length, file=sys.stderr, sep ='\t')

if (options.outdir == "Tenn"):
	#Tenn_demography(S_N1, S_N2, S_AF, S_EU, S_AS, pulses, seed, n1_admix_prop, n2_admix_prop, outdir, t_n1_n2, haplo, length):
	simulation = msprime_demo_models.Tenn_demography(1, 1, 1, 1006, 1008, options.pulses, options.seed, options.n1_admix_prop, options.n2_admix_prop, options.outdir, options.t_n1_n2, options.haplo, options.length)

elif (options.outdir == "Sriram"):
	simulation = msprime_demo_models.Sriram_demography(1, 1, 1, 1006, 1008, options.pulses, options.seed, options.n1_admix_prop, options.n2_admix_prop, options.outdir, options.t_n1_n2, options.haplo, options.length)

elif (options.outdir == "SplitPop"):
	simulation = msprime_demo_models.SplitPop_demography(1, 1, 1, 1006, 1008, options.pulses, options.seed, options.n1_admix_prop, options.n2_admix_prop, options.outdir, options.t_n1_n2, options.haplo, options.length)


###### GET HAPLOTYPES FROM SIMULATED TREES #######

if (options.haplo == "haplo"):
	# Define the sample indices
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
	# Create a .bed file to write to for the simulation
	haplo_outfile = gzip.open( options.outdir+'_'+options.pop+'_'+str(options.seed)+'_n1_'+str(options.n1_admix_prop)+'_n2_'+str(options.n2_admix_prop)+'.bed.gz' , 'wb' )
##  FOR EACH SIMULATED CHROMOSOME, PRINT ALL THE INTROGRESSED HAPLOTYPES BELONGING TO THE SPECIFIED NON-AFR POPULATION IN BED FORMAT  ###
	for t, ts in enumerate(simulation):
	    ## t is the tree ID, ts is a given tree from the simulation
	    node_map = collections.defaultdict(list)
	    ## 'collections' = pythons high-performance container types, 'defaultdict()' is one such container type
	    # defaultdict() = dict subclass that calls a factory function to supply missing values
	    # Using 'list' as the default_factory, it is easy to group a sequence of key-value pairs into a dictionary of lists:
	    ## Defines "node_map" as an empty dictionary of lists ; [('n1', [0, 10]), ('n2', [11, 14]), ('n3', [15, 20])]
	    for record in ts.records():
		## for a single record from all the records of a given tree ts
		# We are interested in coalescence events that occured in the Neandertal population
	        if record.population == 0:
		    ## If the population of the tree is Neand (0)
	            # Since our Neanderthal sample is only two, we can easily exclude events that just concern Neanderthals
		    if set(record.children) != N_samples:
			## If the child/leaf nodes for the record are not Neandertal samples
			# add the node ID and the position of the segments to node_map
			# print("record:", record)
	                node_map[record.node].append((record.left, record.right))
	    for neanderthal_mrca, segments in node_map.items():
		## >>> node_map.items()
		# [('n1', [0, 10]), ('n2', [11, 14]), ('n3', [15, 20])]
		# where the node is the neanderthal_mrca (e.g. 'n1'), and the segments is the tree interval (e.g. [0,10])
#		print(neanderthal_mrca, "->", segments)
	        iterator = introgressed_samples_fn(ts, neanderthal_mrca, N_samples, segments)
		## Run the introgressed_samples function, 
		# using the given tree, the defined Neandertal_mrca/node, the defined tree interval/'segments', and the defined Neand_samples
	        for left, right, samples in iterator:
#			print("\t", left, right, samples)
			for s in samples:
				if s in human_samples:
					if int(math.ceil(left)) < int(math.ceil(right)):
						#print(s, int(math.ceil(left)), int(math.ceil(right)), s)
						haplo_entry=str(s)+'\t'+str(int(math.ceil(left)))+'\t'+str(int(math.ceil(right)))+'\t'+str(s)
						haplo_outfile.write(haplo_entry+'\n')
						# Only print nonAfr samples
						## NOTE: After printing to a bed file, still need to sort and merge the bedfile.
		
	haplo_outfile.close()
	
print('fin',file=sys.stderr)
