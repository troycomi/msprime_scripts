import sys
import math
import pybedtools
import collections
import gzip
import os.path
from AdmixtureOptionParser import admixture_option_parser
import Demography_Models


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


def get_filename(options, extension):
    return "{}_{}_{}_n1_{}_n2_{}".format(
        options.outdir,
        options.pop,
        options.seed,
        options.n1_admix_prop,
        options.n2_admix_prop) \
        + extension


def get_gz_filename(options, base, include_ext=True):
    name = "{}.{}.n1_{}_n2_{}_t_{}_{}".format(
        options.outdir,
        base,
        options.n1_admix_prop,
        options.n2_admix_prop,
        options.t_n1_n2,
        options.seed)
    if include_ext:
        name += ".gz"
    return name


if __name__ == "__main__":
    parser = admixture_option_parser()
    options = parser.parse_args()

    print("Options")
    for key, item in options.__dict__.items():
        print("{}: {}".format(key, item))

    S_N1 = 2
    S_N2 = 2

    models = {
        "Tenn": Demography_Models.Tenn_demography(S_N1, S_N2, options),
        "Sriram": Demography_Models.Sriram_demography(S_N1, S_N2, options),
        "SplitPop": Demography_Models.SplitPop_demography(S_N1, S_N2, options),
    }

    model = models.get(options.outdir, None)
    if model is None:
        print("unsupported model: {}".format(options.outdir))
        sys.exit(1)

    simulation = model.simulate()
    long_names = model.get_long_name_map()

    # combine the sample indices
    non_trgt = S_N1 + S_N2 + options.AF_sample_size
    non_af = options.EU_sample_size + options.AS_sample_size

    if (options.pop == "EAS"):
        human_samples = range(options.EU_sample_size + non_trgt,
                              non_af + non_trgt)

    elif (options.pop == "EUR"):
        human_samples = range(non_trgt, options.EU_sample_size + non_trgt)

    elif (options.pop == "nonAfr"):
        human_samples = range(non_trgt, non_af + non_trgt)

    # GET HAPLOTYPES FROM SIMULATED TREES #######
    if (options.haplo == "haplo"):
        # Create a .bed file to write to for the simulation
        haplo_outfile = gzip.open(
            get_filename(options, '.bed.merged.gz'), 'wb')
        haplo_entry_list = []

        for t, ts in enumerate(simulation):

            # Defines "node_map" as an empty dictionary of lists :
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
                # Run the introgressed_samples function,
                # using the given tree, the defined Neandertal_mrca/node,
                # the defined tree interval/'segments', and the Neand_samples
                iterator = introgressed_samples_fn(ts, neanderthal_mrca,
                                                   range(0, S_N1 + S_N2),
                                                   segments)

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
        with gzip.open(
                get_gz_filename(options, 'eigenstratgeno'),
                'wb') as geno_outfile,\
            gzip.open(
                get_gz_filename(options, 'snp'),
                'wb') as snp_outfile,\
            gzip.open(
                get_gz_filename(options, 'ind'),
                'wb') as ind_outfile,\
            gzip.open(
                'parfile.F4stat.{}.n1_{}_n2_{}_t_{}_{}.gz'.format(
                    options.outdir,
                    options.n1_admix_prop,
                    options.n2_admix_prop,
                    options.t_n1_n2,
                    options.seed
                ), 'wb') as parF4_outfile:

            rs_num = 0

            # enumerate creates a list of tuples, where each chrom
            # is a tuple with an index and the tree
            # e.g. [(1,Tree1), (2,Tree2)...]
            for t, tree_sequence in enumerate(simulation):
                #  WRITE THE .IND FILE  ###
                # Write .ind file based on output of Tree1
                if t < 1:
                    # Get the sample size from the tree
                    for i in range(tree_sequence.get_sample_size()):
                        pop = long_names[tree_sequence.get_population(i)]

                        ind_outfile.write(
                            str.encode('Sample_{}\tU\t{}\n'.format(i, pop)))

            #  WRITE THE .EIGENSTRATGENO AND .SNP FILES  ###
            chr_num = t+1
            for variant in tree_sequence.variants(as_bytes=True):
                rs_num += 1

                # write genotypes to .eigenstratgeno file
                geno_outfile.write(variant.genotypes+b'\n')

                # write snp_allele info to .snp file
                snp_outfile.write(
                    str.encode(
                        'rs{}\t{}\t{}\t{}\tA\tT\n'.format(
                            rs_num,
                            chr_num,
                            variant.position / options.length,
                            int(variant.position))))

            parF4_outfile.write(
                str.encode('genotypename: {}\n'.format(
                    get_gz_filename(options, 'eigenstratgeno', False))))

            parF4_outfile.write(
                str.encode('snpname: {}\n'.format(
                    get_gz_filename(options, 'snp', False))))

            parF4_outfile.write(
                str.encode('indivname: {}\n'.format(
                    get_gz_filename(options, 'ind', False))))

            parF4_outfile.write(
                str.encode('popfilename: sim.popfile_F4stat'+'\n'))


#    WRITE VCF FILE FORMAT FOR S* CALCULATIONS , ALONG WITH POPULATION FILE
    elif (options.haplo == "vcf"):
        vcf_outfile = open(get_filename(options, '.vcf'), 'w')
        for t, tree_sequence in enumerate(simulation):
            # If a population file does not exist yet, write one
            if os.path.isfile(options.outdir+'.popfile'):
                print('popfile already exists', file=sys.stderr)
            else:
                # Tenn.popfile.gz
                pop_outfile = open(options.outdir+'.popfile', 'w')
                # write header to pop_outfile
                pop_outfile.write('samp'+'\t'+'pop'+'\t'+'super_pop'+'\n')

                # For each individual in Tree1
                for i in range(tree_sequence.get_sample_size()):
                    if i % 2 == 0:
                        pop_outfile.write(
                            'msp_{0}\t{1}\t{1}\n'.format(
                                i//2,
                                long_names[tree_sequence.get_population(i)]))

                pop_outfile.close()
            tree_sequence.write_vcf(vcf_outfile, 2)

        vcf_outfile.close()

        # have to gzip seperately
        with open(get_filename(options, '.vcf'), 'r') as f_in,\
                gzip.open(get_filename(options, '.vcf.gz'), 'wb') as f_out:

            for lin in f_in.readlines():
                f_out.write(str.encode(lin))

        os.remove(get_filename(options, '.vcf'))
