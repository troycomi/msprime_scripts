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
parser.add_option("-t", "--taskid", action = "store", type = "string", dest = "taskid")
parser.add_option("-i", "--introgress_pulses", action = "store", type = "int", dest = "pulses")
(options, args) = parser.parse_args()

def simulate_demography(S_N1, S_N2, S_AF, S_EU, S_AS, pulses):
    N_A = 7310 # effective population size used for scaling
    N_N1 = 1000 # Altai/Vindija lineage effective population size
    N_N2 = 1000 # Eastern Neandertal effective population size
    N_B = 2100
    N_AF = 12300
    N_EU0 = 1000
    N_AS0 = 510
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_MH_N = 700e3 / generation_time # Neanderthal / modern human split
    T_N1_N2 = 350e3 / generation_time # Altai/Vindija and Eastern Neandertal popualation split
    T_AF = 220e3 / generation_time # African population expansion
    T_B = 140e3 / generation_time # out of Africa
    T_PULSE1 = 55e3 / generation_time # Neandertal introgression pulse 1
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
        initial_size = N_AS, growth_rate = r_AS) # population 4
    ]
    # specify desired sample
    N1_sample = [msprime.Sample(population = 0, time = 100e3 / generation_time)] * S_N1
    N2_sample = [msprime.Sample(population = 1, time = 100e3 / generation_time)] * S_N2
    AF_sample = [msprime.Sample(population = 2, time = 0)] * S_AF
    EU_sample = [msprime.Sample(population = 3, time = 0)] * S_EU
    AS_sample = [msprime.Sample(population = 4, time = 0)] * S_AS
    samples = N1_sample + N2_sample + AF_sample + EU_sample + AS_sample
    # set up the initial migration matrix
    migration_matrix = [
        [      0, m_N1_N2, m_N1_AF, m_N1_EU, m_N1_AS],
        [m_N1_N2,       0, m_N2_AF, m_N2_EU, m_N2_AS],
        [m_N1_AF, m_N2_AF,       0, m_AF_EU, m_AF_AS],
        [m_N1_EU, m_N2_EU, m_AF_EU,       0, m_EU_AS],
        [m_N1_AS, m_N2_AS, m_AF_AS, m_EU_AS,       0],
    ]
    one_pulse = [
        msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (2, 3)), # migration between "B" and Africa begins
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)),
        msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
        msprime.MigrationRateChange(time = T_PULSE1, rate = 1e-3, matrix_index = (3, 0)), # introgression
        msprime.MigrationRateChange(time = T_PULSE1 + 100, rate = 0, matrix_index = (3, 0)), # introgression ends
        msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
        msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
        msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
        msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2) # set parameters of ancetral hominin population
    ]
    two_pulse = [
        msprime.MigrationRateChange(time = T_PULSE2, rate = 2e-04, matrix_index = (4, 1)), # introgression
        msprime.MigrationRateChange(time = T_PULSE2 + 100, rate = 0, matrix_index = (4, 1)), # introgression ends
        msprime.MassMigration(time = T_EU_AS, source = 4, destination = 3, proportion = 1.0), # AS merges into EU, now termed "B"
        msprime.MigrationRateChange(time = T_EU_AS, rate = 0), # set all migration rates to zero
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (2, 3)), # migration between "B" and Africa begins
        msprime.MigrationRateChange(time = T_EU_AS, rate = m_AF_B, matrix_index = (3, 2)),
        msprime.PopulationParametersChange(time = T_EU_AS, initial_size = N_B, growth_rate = 0, population_id = 3), # set parameters of population "B"
        msprime.MigrationRateChange(time = T_PULSE1, rate = 1e-3, matrix_index = (3, 0)), # introgression
        msprime.MigrationRateChange(time = T_PULSE1 + 100, rate = 0, matrix_index = (3, 0)), # introgression ends
        msprime.MassMigration(time = T_B, source = 3, destination = 2, proportion = 1.0), # Population B merges into Africa at T_B
        msprime.MigrationRateChange(time = T_B, rate = 0), # set all migration rates to zero
        msprime.PopulationParametersChange(time = T_AF, initial_size = N_A, population_id = 2), # set parameters of ancestral modern human population
        msprime.MassMigration(time = T_N1_N2, source = 1, destination = 0, proportion = 1.0), # N_2 merges with N_1 at T_N1_N2
        msprime.PopulationParametersChange(time = T_N1_N2, initial_size = N_N1, population_id = 0), # set parameters of ancestral Neandertal population
        msprime.MassMigration(time = T_MH_N, source = 0, destination = 2, proportion = 1.0), # Neandertals merge into modern human lineage at time T_MH_N
        msprime.PopulationParametersChange(time = T_MH_N, initial_size = N_A, population_id = 2) # set parameters of ancetral hominin population
    ]
    if (pulses == 1):
        demographic_events = one_pulse
    elif (pulses == 2):
        demographic_events = two_pulse
    length = 1e6
    recombination_rate = 2e-8
    mutation_rate = 1e-8 # see Shendure and Akey, 2015
    return msprime.simulate(Ne = N_A,
    length = length,
    recombination_rate = recombination_rate,
    mutation_rate = mutation_rate,
    samples = samples,
    population_configurations = population_configurations,
    migration_matrix = migration_matrix,
    demographic_events = demographic_events)

# run the simulation
simulation = simulate_demography(1, 1, 10, 1006, 1008, options.pulses)
# define the sample indices
N_samples  = range(0, 2)
AF_samples = range(2, 12)
EU_samples = range(12, 1018)
AS_samples = range(1018, 2026)

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
        quit()
    # store the haplotypes in a list
    haplotypes = []
    for haplotype in tree_sequence.haplotypes():
        haplotypes.append(haplotype)
    if (len(haplotypes) == 0):
        quit()
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
        quit()
    interval_array = np.concatenate(introgression_list)
    # use bedtools to operate on the introgressed intervals
    prefix = outdir + str(taskid) + "/"
    COMMAND = "mkdir " + prefix
    return_code = subprocess.call(COMMAND, shell = True)
    np.savetxt(prefix + str(taskid) + ".bed", interval_array, delimiter = "\t", newline = "\n", fmt = "%s")
    COMMAND = '''awk '{print >> "''' + prefix + str(taskid) + '''."$4".bed"}' ''' + prefix + str(taskid) + ".bed; rm " + prefix + str(taskid) + ".bed"
    return_code = subprocess.call(COMMAND, shell = True)
    COMMAND = '''for i in ''' + prefix + str(taskid) + '''*.bed; do sort -k2,2n $i > $i.sorted; rm $i; done'''
    return_code = subprocess.call(COMMAND, shell = True)
    COMMAND = '''for i in ''' + prefix + str(taskid) + '''*.sorted; do bedtools merge -i $i -c 4 -o distinct > $i.merged; rm $i; done'''
    return_code = subprocess.call(COMMAND, shell = True)
    # concatenate the individual sample files
    files = prefix + str(taskid) + "*.merged"
    read_files = glob.glob(files)
    fid = prefix + str(taskid) + ".bed"
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

int_results = get_introgressed_intervals(simulation, human_samples, N_samples, options.outdir, options.taskid)

# choose a random introgressed haplotype and get all overlapping haplotypes
sample_index = random.sample(range(len(int_results.intervals)), 1)[0]
random_interval = pybedtools.BedTool(str(int_results.intervals[sample_index]), from_string = True)
olaps = random_interval.intersect(int_results.intervals, f = 1, wb = True)
start = random_interval[0].start
end = random_interval[0].end
length = end - start

# find the start and end coordinates of the mutations
m_indices = [index for index, item in enumerate(int_results.mutations) if item >= start and item <= end]
if (len(m_indices) == 0):
    quit()

m_start = m_indices[0]
m_end = m_indices[-1]

# get all overlapping haplotypes
int_haps = []
for sample in olaps:
    s_index = int(sample.fields[7])
    int_hap = map(int, int_results.haplotypes[s_index][m_start:m_end + 1])
    int_haps.append(int_hap)

int_haps = np.matrix(int_haps)
n_haps = int_haps.shape[0]
if (n_haps < 2):
    quit()

manhattan = scipy.spatial.distance.pdist(int_haps, metric = "cityblock")
pi_mean = np.mean(manhattan)
pi_min = np.min(manhattan)
pi_max = np.max(manhattan)

# get divergence of each haplotype to the Neandertal haplotype
neand = map(int, int_results.haplotypes[N_samples[0]][m_start:m_end + 1])
ndists = []
for row in int_haps[range(n_haps)]:
    ndists.append(scipy.spatial.distance.cityblock(row, neand))

ndist_mean = np.mean(ndists)
ndist_min = np.min(ndists)
ndist_max = np.max(ndists)

print(length, pi_mean, pi_min, pi_max, ndist_mean, ndist_min, ndist_max, n_haps, options.pop, sep = "\t")
