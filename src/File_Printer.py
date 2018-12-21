import sys
import io
import gzip
import os.path
import pybedtools
from Option_Parser import admixture_option_parser


def get_basename(file):
    return os.path.splitext(
        os.path.basename(file))[0]


class file_printer(object):
    '''object to validate and handle the output of the various files
    from the admixture_option_parser'''
    def __init__(self, options: admixture_option_parser):
        '''parse in the provided inputs, validate the set
        and build filenames'''
        self.options = options
        self.files = {}

        self.files['debug'] = options.debug_file
        self.files['haplo'] = options.haplo_file
        self.files['ils'] = options.ils_file
        self.files['option'] = options.option_file
        self.files['vcf'] = options.vcf_file
        self.files['f4dstat'] = options.f4dstat_file
        self.out_dir = options.out_dir

        self.validate_options()

    def __enter__(self):
        self.build_files()
        self.open_writers()
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close_writers()

    def validate_options(self):
        '''test if the provided inputs are valid, raise ValueError
        otherwise'''
        # No values set
        if all([f is None for f in self.files.values()]):
            return

        # failure modes
        numset = sum([f == '*' for f in self.files.values()])
        if numset > 1:
            raise ValueError('Expected at most one output to stdout, '
                             'got {} instead.'.format(numset))

    def build_files(self):
        '''generate all filenames needed for printing outputs'''
        out_dir_set = self.out_dir is not None
        numset = sum([f is not None for f in self.files.values()])
        # when outdir is set without other options, print everything
        # otherwise print only set files in outdir
        print_all = numset == 0 and out_dir_set

        self.build_out_dir()
        self.build_debug(numset, print_all)
        self.build_haplo(print_all)
        self.build_ils(print_all)
        self.build_option(print_all)
        self.build_vcf(print_all)
        self.build_f4dstat(print_all)

    class file_struct:
        def __init__(self, default, fmt="{}"):
            self.default = default
            self.fmt = fmt

        def non_default(self, base):
            return self.fmt.format(base)

    def build_out_dir(self):
        if self.out_dir is not None:
            if not os.path.exists(self.out_dir):
                os.mkdir(self.out_dir)
        else:
            self.out_dir = os.getcwd()

    def build_generic(self, inputKey, defaults, print_all,
                      numset=-1, allow_stdout=False):
        '''helper function to handle building filenames
        defaults is a dict of output keys and file_structs'''

        filename = self.files[inputKey]

        if allow_stdout:
            assert len(defaults) == 1  # only one to stdout
            if numset == 0 and not print_all:  # default to stdout
                for key, value in defaults.items():
                    self.files[key] = sys.stdout
                return

            elif filename == '*':  # print to stdout regardless of out_dir
                for key, value in defaults.items():
                    self.files[key] = sys.stdout
                return

        if filename is None:
            if print_all:  # print default name to out file
                for key, value in defaults.items():
                    self.files[key] = os.path.join(
                        self.out_dir,
                        value.default)

            else:  # not set and not printing all, do nothing
                for key, value in defaults.items():
                    self.files[key] = None

        else:  # set with filename, join with out_dir
            for key, value in defaults.items():
                self.files[key] = os.path.join(
                    self.out_dir,
                    value.non_default(filename))

    def build_debug(self, numset, print_all):
        self.build_generic('debug',
                           {'debug': self.file_struct('debug.txt')},
                           print_all,
                           numset,
                           allow_stdout=True)

    def build_haplo(self, print_all):
        self.build_generic('haplo',
                           {'haplo':
                            self.file_struct(
                                self.get_filename('.bed.merged.gz'))},
                           print_all,
                           allow_stdout=True)

    def build_ils(self, print_all):
        self.build_generic('ils',
                           {'ils':
                            self.file_struct(
                                self.get_filename('.ils.bed.merged.gz'))},
                           print_all,
                           allow_stdout=True)

    def build_option(self, print_all):
        self.build_generic('option',
                           {'option': self.file_struct('options.txt')},
                           print_all,
                           allow_stdout=True)

    def build_vcf(self, print_all):
        self.build_generic('vcf',
                           {'vcf':
                            self.file_struct(
                                self.get_filename('.vcf.gz'),
                                "{}.vcf.gz"),
                            'popfile':
                            self.file_struct(
                                self.get_filename('.popfile'),
                                "{}.popfile")},
                           print_all)

    def build_f4dstat(self, print_all):
        self.build_generic('f4dstat',
                           {'f4dstat':
                            self.file_struct(
                                self.get_f4dstat_filename('parfile.F4stat'),
                                "parfile.F4stat.{}.gz"),
                            'eigen':
                            self.file_struct(
                                self.get_f4dstat_filename('eigenstratgeno'),
                                "eigenstratgeno.{}.gz"),
                            'snp':
                            self.file_struct(
                                self.get_f4dstat_filename('snp'),
                                "snp.{}.gz"),
                            'ind':
                            self.file_struct(
                                self.get_f4dstat_filename('ind'),
                                "ind.{}.gz"),
                            },
                           print_all)

    def get_filename(self,  extension):
        options = self.options
        return "{model}_{pop}_{seed}_n1_{n1}_n2_{n2}".format(
            model=options.model,
            pop=options.pop,
            seed=options.seed,
            n1=options.n1_admix_prop,
            n2=options.n2_admix_prop) \
            + extension

    def get_f4dstat_filename(self, base):
        options = self.options
        return "{model}.{base}.n1_{n1}_n2_{n2}_t_{tn1n2}_{seed}.gz".format(
            model=options.model,
            base=base,
            n1=options.n1_admix_prop,
            n2=options.n2_admix_prop,
            tn1n2=options.t_n1_n2,
            seed=options.seed)

    def open_writers(self):
        self.writers = {}
        for key, filename in self.files.items():
            if filename is None:
                writer = None
            elif filename == sys.stdout:
                writer = sys.stdout
            else:
                if os.path.splitext(filename)[1] == '.gz':
                    writer = gzip.open(filename, 'wb')
                else:
                    writer = open(filename, 'w')

            self.writers[key] = writer

    def close_writers(self):
        for writer in self.writers.values():
            if writer != sys.stdout and writer is not None:
                writer.close()

    def print_options(self):
        writer = self.writers['option']

        if writer is not None:
            writer.write("Options\n")
            for key, item in self.options.__dict__.items():
                if isinstance(item, list):
                    writer.write("{}: [\n".format(key))
                    for el in item:
                        writer.write("\t{},\n".format(el))
                    writer.write("]\n")
                else:
                    writer.write("{}: {}\n".format(key, item))

    def print_debug(self, model):
        writer = self.writers['debug']

        if writer is not None:
            model.print_debug(writer)

    def print_popfile(self, model, tree_sequence):
        writer = self.writers['popfile']

        if writer is None:
            return

        long_names = model.get_long_name_map()
        writer.write('samp\tpop\tsuper_pop\n')

        # For each individual in Tree1
        for i in range(0, tree_sequence.get_sample_size(), 2):
            writer.write(
                'msp_{0}\t{1}\t{1}\n'.format(
                    i//2,
                    long_names[tree_sequence.get_population(i)]))

    def single_simulation_needed(self):
        return self.vcf_needed() or \
            self.haplo_needed() or \
            self.ils_needed()

    def vcf_needed(self):
        return self.writers['vcf'] is not None

    def print_vcf(self, tree_sequence):
        writer = self.writers['vcf']

        if writer is None:
            return

        with io.StringIO() as tempfile:
            tree_sequence.write_vcf(tempfile, 2)
            tempfile.seek(0)

            for line in tempfile.readlines():
                writer.write(str.encode(line))

    def haplo_needed(self):
        return self.writers['haplo'] is not None

    def print_haplo(self, haplo_entry_list):
        self.print_haplo_helper('haplo', haplo_entry_list)

    def ils_needed(self):
        return self.writers['ils'] is not None

    def print_ils(self, haplo_entry_list):
        self.print_haplo_helper('ils', haplo_entry_list)

    def print_haplo_helper(self, writer, haplo_entry_list):
        writer = self.writers[writer]

        if writer is None:
            return

        # NOTE: After printing to a bed file, still need to sort and merge
        haplo_entry_string = '\n'.join(haplo_entry_list)
        pybedtools.set_tempdir('/scratch/')
        bedfile = pybedtools.BedTool(haplo_entry_string, from_string=True)
        bedfile_sorted_merged = pybedtools.BedTool.sort(bedfile).merge()

        # Read the BEDfile line by line, add in the chr#,
        # reorder so that the columns are: chr, strt, end, ind
        for bed_line in bedfile_sorted_merged:
            new_bed_line = '{}\t{}'.format(
                str(bed_line).strip(),
                self.options.seed)
            new_bed_line = new_bed_line.split('\t')
            new_bed_line[0], new_bed_line[3] = new_bed_line[3], new_bed_line[0]
            writer.write(str.encode('\t'.join(new_bed_line)+'\n'))

        pybedtools.cleanup(verbose=False)

    def f4dstat_needed(self):
        return self.writers['f4dstat'] is not None

    def write_to(self, writer, line):
        '''Note that unlike other printers, makes no check that
        writer is open and valid'''

        self.writers[writer].write(line)

    def print_f4dstat(self):
        writer = self.writers['f4dstat']
        if writer is None:
            return

        writer.write(
            str.encode('genotypename: {}\n'.format(
                get_basename(self.files['eigen']))))

        writer.write(
            str.encode('snpname: {}\n'.format(
                get_basename(self.files['snp']))))

        writer.write(
            str.encode('indivname: {}\n'.format(
                get_basename(self.files['ind']))))

        writer.write(
            str.encode('popfilename: sim.popfile_F4stat'+'\n'))
