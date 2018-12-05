import sys
import math
import pybedtools
import collections
import gzip
import os.path
from AdmixtureOptionParser import admixture_option_parser
import Demography_Models


class introgression_object(object):
    def __init__(self, mutations, haplotypes, intervals):
        self.mutations = mutations
        self.haplotypes = haplotypes
        self.intervals = intervals


def introgressed_samples_fn(ts, neanderthal_mrca,
                            neanderthal_samples, segments):
    # Define the samples that carry introgressed segments?
    trees = ts.trees()
    tree = next(trees)
    for left, right in segments:
        # Skip ahead to the tree that intersects with this segment.
        # NOTE: Regading tree data:
        #  The "Records" group consists of four pieces of information:
        #   1) the left and 2) right coordinates of the coalescing interval,
        #   3) the list of child nodes (modern human samples) and \
        #   4) the parent node.
        #  Each record returned can be accessed via the attributes:
        #   'left', 'right', 'node', 'children', 'time' and 'population'
        #  A record represents the assignment of a pair of children 'c'
        #   to a parent 'u'. This assignment happens at 't' generations
        #   in the past within the population with ID 'd'
        while tree.get_interval()[0] < left:
            # Pass through the potential trees until you
            # get to the one where the start_position matches
            # the start of the defined segment
            tree = next(trees)

        assert tree.get_interval()[0] == left
        start = None
        last_human_leaves = None

        while tree is not None and tree.get_interval()[1] <= right:
            human_leaves = set(tree.leaves(neanderthal_mrca))\
                           - set(neanderthal_samples)
            # Create a set() of all the human leaves represented in this tree
            #  that share an mrca in the Neandertal lineage
            #  (i.e. mrca node is in Neand population)
            # neanderthal_mrca is supplied from the node_map dictionary
            #  as the parent node for the tree present
            #  in the Neandertal population
            # tree.leaves(n) returns an iterator over all the leaves
            #  in this tree underneath the specified node (in this
            #  case the neanderthal_mrca)
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


# SIMULATION RUNS ############
if __name__ == "__main__":
    parser = admixture_option_parser()
    options = parser.parse_args()

    print("Model: ", options.outdir, "Seed: ", options.seed, "Neand_Pulse1: ",
          options.n1_admix_prop, "Neand_Pulse2: ", options.n2_admix_prop,
          "Length: ", options.length, file=sys.stderr, sep='\t')

    S_N1 = 2
    S_N2 = 2

    if (options.outdir == "Tenn"):
        simulation = Demography_Models.Tenn_demography(
                S_N1, S_N2, options).simulate()

    elif (options.outdir == "Sriram"):
        simulation = Demography_Models.Sriram_demography(
                S_N1, S_N2, options).simulate()

    elif (options.outdir == "SplitPop"):
        simulation = Demography_Models.SplitPop_demography(
                S_N1, S_N2, options).simulate()


    # Define the sample indices
    EU_count = options.EU_sample_size
    AS_count = options.AS_sample_size
    AF_count = options.AF_sample_size
    N_samples = range(0, S_N1+S_N2)        # range(0,4) --> 0,1,2,3
    AF_samples = range(S_N1+S_N2, S_N1+S_N2+AF_count)   # range(4,6) --> 4,5
    non_trgt = len(N_samples) + len(AF_samples)      # 6
    nonAfr_samples = range(non_trgt, (EU_count + AS_count + non_trgt))
    EU_samples = range(non_trgt, (EU_count + non_trgt))
    AS_samples = range((EU_count + non_trgt), (EU_count+AS_count+non_trgt))
    Chimp_samples = range((EU_count+AS_count+non_trgt),
                          (EU_count+AS_count+non_trgt)+2)
    Deni_samples = range((EU_count+AS_count+non_trgt) + 2,
                         (EU_count+AS_count+non_trgt) + 2 + 2)

    if (options.pop == "EAS"):
        human_samples = AS_samples
    elif (options.pop == "EUR"):
        human_samples = EU_samples
    elif (options.pop == "nonAfr"):
        human_samples = nonAfr_samples


    # GET HAPLOTYPES FROM SIMULATED TREES #######

    if (options.haplo == "haplo"):
        # Create a .bed file to write to for the simulation
        haplo_outfile = gzip.open(options.outdir + '_' + options.pop + '_' +
                                  str(options.seed) + '_n1_' +
                                  str(options.n1_admix_prop) + '_n2_' +
                                  str(options.n2_admix_prop) +
                                  '.bed.merged.gz', 'wb')
        haplo_entry_list = []
        # FOR EACH SIMULATED CHROMOSOME,
        #  PRINT ALL THE INTROGRESSED HAPLOTYPES BELONGING TO
        #   THE SPECIFIED NON-AFR POPULATION IN BED FORMAT
        # t is the tree ID, ts is a given tree from the simulation
        for t, ts in enumerate(simulation):

            # 'collections' = pythons high-performance container types,
            # 'defaultdict()' is one such container type
            # defaultdict() = dict subclass that calls a factory
            # function to supply missing values
            # Using 'list' as the default_factory, it is easy to
            # group a sequence of key-value pairs into a dictionary of lists:
            # Defines "node_map" as an empty dictionary of lists :
            # [('n1', [0, 10]), ('n2', [11, 14]), ('n3', [15, 20])]
            node_map = collections.defaultdict(list)

            # for a single record from all the records of a given tree ts
            # We are interested in coalescence events that occured
            #  in the Neandertal population
            for record in ts.records():

                # If the population of the tree is Neand (0 or 1)
                # Since our Neanderthal sample is only two,
                #  we can easily exclude events that just concern Neanderthals
                if record.population <= 1:
                    node_map[record.node].append((record.left, record.right))

            for neanderthal_mrca, segments in node_map.items():
                # >>> node_map.items()
                # [('n1', [0, 10]), ('n2', [11, 14]), ('n3', [15, 20])]
                # where the node is the neanderthal_mrca (e.g. 'n1'),
                #  and the segments is the tree interval (e.g. [0,10])

                # Run the introgressed_samples function,
                # using the given tree, the defined Neandertal_mrca/node,
                # the defined tree interval/'segments', and the Neand_samples
                iterator = introgressed_samples_fn(ts, neanderthal_mrca,
                                                   N_samples, segments)

                for left, right, samples in iterator:
                    for s in samples:
                        if s in human_samples and \
                                math.ceil(left) < math.ceil(right):
                            haplo_entry = str(s)\
                                    + '\t' + str(int(math.ceil(left)))\
                                    + '\t' + str(int(math.ceil(right)))\
                                    + '\t' + str(s)
                            haplo_entry_list.append(haplo_entry + '\n')
        # NOTE: After printing to a bed file, still need to sort and merge
        # Join together the list of introgressed haplotypes into a string,
        #  and convert this to a BED file using pybedtools
        # Then perform the sort() and merge() functions in python
        haplo_entry_string = ''.join(haplo_entry_list)
        pybedtools.set_tempdir('/scratch/')
        BEDFILE = pybedtools.BedTool(haplo_entry_string, from_string=True)
        BEDFILE_SORTED_MERGED = pybedtools.BedTool.sort(BEDFILE).merge()

        # Read the BEDfile line by line, add in the chr#,
        # reorder so that the columns are: chr, strt, end, ind
        # write to haplo_outfile in gzip
        for bed_line in BEDFILE_SORTED_MERGED:
            new_bed_line = str(bed_line).strip() + '\t' + str(options.seed)
            new_bed_line = new_bed_line.split('\t')
            new_bed_line[0], new_bed_line[3] = new_bed_line[3], new_bed_line[0]
            haplo_outfile.write(str.encode('\t'.join(new_bed_line)+'\n'))

        pybedtools.cleanup(verbose=True)
        haplo_outfile.close()

#    GET EIGENSTRATGENO FILES AND SNP FILES    #########

    elif (options.haplo == "F4Dstat"):
            outdir = options.outdir
            n1_admix_prop = options.n1_admix_prop
            n2_admix_prop = options.n2_admix_prop
            t_n1_n2 = options.t_n1_n2
            length = options.length
            seed = options.seed

            #  sim.eigenstratgeno.n1_0.01_n2_0.05_200
            geno_outfile = gzip.open(outdir + '.eigenstratgeno.n1_' +
                                     str(n1_admix_prop) + '_n2_' +
                                     str(n2_admix_prop) + '_t_' +
                                     str(t_n1_n2) + '_' + str(seed) +
                                     '.gz', 'wb')
            snp_outfile = gzip.open(outdir + '.snp.n1_' +
                                    str(n1_admix_prop) + '_n2_' +
                                    str(n2_admix_prop) + '_t_' +
                                    str(t_n1_n2) + '_' + str(seed) +
                                    '.gz', 'wb')
            ind_outfile = gzip.open(outdir + '.ind.n1_' + str(n1_admix_prop)
                                    + '_n2_' + str(n2_admix_prop) +
                                    '_t_' + str(t_n1_n2) + '_' +
                                    str(seed) + '.gz', 'wb')
            rs_num = 0

            # enumerate creates a list of tuples, where each chrom
            # is a tuple with an index and the tree
            # e.g. [(1,Tree1), (2,Tree2)...]
            for t, tree_sequence in enumerate(simulation): 
                #  WRITE THE .IND FILE  ###
                # Write .ind file based on output of Tree1
                if t<1:
                    # Get the sample size from the tree
                    for i in range(tree_sequence.get_sample_size()):        
                        # For each individual in Tree1,
                        # get the corresponding population
                        if tree_sequence.get_population(i) == 0:
                            pop = 'Neand1'
                        elif tree_sequence.get_population(i) == 1:
                            pop = 'Neand2'
                        elif tree_sequence.get_population(i) == 2:
                            pop = 'AFR'
                        elif tree_sequence.get_population(i) == 3:
                            pop = 'EUR'
                        elif tree_sequence.get_population(i) == 4:
                            pop = 'ASN'
                        elif tree_sequence.get_population(i) == 5:
                            pop = 'Chimp'
                        elif tree_sequence.get_population(i) == 6:
                            pop = 'Deni'
                        elif tree_sequence.get_population(i) == 7:
                            pop = 'Split'

                        # Write the sample entry using the population from Tree1
                        ind_entry = 'Sample_' + str(i) + '\t' + 'U' + '\t' + pop
                        ind_outfile.write(ind_entry+'\n')
            #  WRITE THE .EIGENSTRATGENO AND .SNP FILES  ###
            chr_num = t+1
            for variant in tree_sequence.variants(as_bytes=True):
                rs_num += 1
                # write genotypes to .eigenstratgeno file
                geno_outfile.write(variant.genotypes+'\n')
                line = str('rs' + str(rs_num) + '\t' + str(chr_num) +
                           '\t' + str(variant.position/length) + '\t' +
                           str(int(variant.position)) + '\t' + 'A' +
                           '\t' + 'T' + '\n')
                # write snp_allele info to .snp file
                snp_outfile.write(line)

            parF4_outfile = gzip.open('parfile.F4stat.' + outdir + '.n1_'
                                      + str(n1_admix_prop) + '_n2_' +
                                      str(n2_admix_prop) + '_t_' +
                                      str(t_n1_n2) + '_' + str(seed) +
                                      '.gz', 'wb')
            parF4_outfile.write('genotypename: ' + outdir +
                                '.eigenstratgeno.n1_' + str(n1_admix_prop) +
                                '_n2_' + str(n2_admix_prop) + '_t_' +
                                str(t_n1_n2) + '_' + str(seed) + '\n')
            parF4_outfile.write('snpname: ' + outdir + '.snp.n1_' +
                                str(n1_admix_prop) + '_n2_' +
                                str(n2_admix_prop) + '_t_' +
                                str(t_n1_n2) + '_' +
                                str(seed) + '\n')
            parF4_outfile.write('indivname: ' + outdir + '.ind.n1_' +
                                str(n1_admix_prop) + '_n2_' +
                                str(n2_admix_prop) + '_t_' +
                                str(t_n1_n2) + '_' + str(seed) + '\n')
            parF4_outfile.write('popfilename: sim.popfile_F4stat'+'\n')

            geno_outfile.close()
            snp_outfile.close()
            ind_outfile.close()
            parF4_outfile.close()


#    WRITE VCF FILE FORMAT FOR S* CALCULATIONS , ALONG WITH POPULATION FILE

    elif (options.haplo == "vcf"):
            outdir = options.outdir
            n1_admix_prop = options.n1_admix_prop
            n2_admix_prop = options.n2_admix_prop
            pop = options.pop
            t_n1_n2 = options.t_n1_n2
            length = options.length
            seed = options.seed

            for t, tree_sequence in enumerate(simulation):
                # If a population file does not exist yet, write one
                if os.path.isfile(outdir+'.popfile'):
                    print('popfile already exists', file=sys.stderr)
                else:
                    # Tenn.popfile.gz
                    pop_outfile = open(outdir+'.popfile', 'w')
                    # write header to pop_outfile
                    pop_outfile.write('samp'+'\t'+'pop'+'\t'+'super_pop'+'\n')

                    # For each individual in Tree1
                    for i in range(tree_sequence.get_sample_size()):
                        # get the corresponding population
                        if tree_sequence.get_population(i) == 0:
                            s_pop = 'Neand1'
                        elif tree_sequence.get_population(i) == 1:
                            s_pop = 'Neand2'
                        elif tree_sequence.get_population(i) == 2:
                            s_pop = 'AFR'
                        elif tree_sequence.get_population(i) == 3:
                            s_pop = 'EUR'
                        elif tree_sequence.get_population(i) == 4:
                            s_pop = 'ASN'
                        elif tree_sequence.get_population(i) == 5:
                            s_pop = 'Chimp'
                        elif tree_sequence.get_population(i) == 6:
                            s_pop = 'Deni'
                        elif tree_sequence.get_population(i) == 7:
                            s_pop = 'Split'
                        if i % 2 == 0:
                            # Write, as a string, the sample entry using
                            # the population read from Tree1
                            ind_entry = 'msp_' + str(i/2) + '\t' + s_pop +\
                                        '\t' + s_pop
                            pop_outfile.write(ind_entry + '\n')
                    pop_outfile.close()

                # Tenn_nonAfr_100_n1_0.01_n2_0.01.vcf.gz
                vcf_outfile = gzip.open(outdir + '_' + pop + '_' +
                                        str(seed) + '_n1_' +
                                        str(n1_admix_prop) + '_n2_' +
                                        str(n2_admix_prop) +
                                        '.vcf.gz', 'wb')
                tree_sequence.write_vcf(vcf_outfile, 2)

            vcf_outfile.close()

    print('fin', file=sys.stderr)
