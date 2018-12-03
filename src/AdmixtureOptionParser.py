import sys
from argparse import ArgumentParser


class admixture_option_parser(ArgumentParser):
    def __init__(self):
        description = "Simulate admixture"
        ArgumentParser.__init__(self, description=description)

        self.add_argument("-p", "--population",
                        dest="pop",
                        default="nonAfr",
                        help="call introg_haplotypes in EUR, EAS, \
                                  or all nonAfr; default=nonAfr")

        self.add_argument("-o", "--outdir",
                        dest="outdir",
                        default="Tenn",
                        help="output directory name, \
                              best to include model name; \
                              default=Tenn")

        self.add_argument("-s", "--seed",
                        type=int,
                        dest="seed",
                        default=1,
                        help="Set random seed for replicate chromosomes; \
                          default=1")

        self.add_argument("-i", "--introgress_pulses",
                        type=int,
                        dest="pulses",
                        default=2,
                        help="Set number of introgression pulses; default=2")

        self.add_argument("-n", "--neand1_admixture_proportion",
                        type=float,
                        dest="n1_admix_prop",
                        default=0.02,
                        help="Set N1 admixture proportion; default=0.02")

        self.add_argument("-d", "--deni_admixture_propotion",
                        type=float,
                        dest="n2_admix_prop",
                        default=0.0,
                        help="Set N2 admixture proportion; default=0.0")

        self.add_argument("-t", "--time_N1_N2_split",
                        type=int,
                        dest="t_n1_n2",
                        default=350,
                        help="Set N1 N2 split time in kya; default=350")

        self.add_argument("-c", "--calls_or_stats",
                        dest="haplo",
                        default="F4Dstat",
                        help="Pick haplotype calls ('haplo') or 'F4Dstat' \
                            or 'vcf' output, or 'debug'; default=F4Dstat")

        self.add_argument("-l", "--length_chrom",
                        type=float,
                        dest="length",
                        default=1e6,
                        help="Define length of simulated chromosome ; \
                        default=1e6")

        self.add_argument("-e", "--european_sample_size",
                        type=int,
                        dest="EU_sample_size",
                        default=1006,
                        help="Set EU haploid sample size ; default=1006")

        self.add_argument("-a", "--asian_sample_size",
                        type=int,
                        dest="AS_sample_size",
                        default=1008,
                        help="Set AS haploid sample size ; default=1008")

        self.add_argument("-r", "--reference_sample_size",
                        type=int,
                        dest="AF_sample_size",
                        default=2,
                        help="Set AF haploid sample size ; default=2")

        self.add_argument("--migration_AF_B",
                        type=float,
                        dest="m_AF_B",
                        default=15e-5,
                        help="Set African <---> ancestral Eurasian \
                              migration rate ; default=15e-5")

        self.add_argument("--migration_AF_AS",
                        type=float,
                        dest="m_AF_AS",
                        default=0.78e-5,
                        help="Set African <-> Asian migration rate ;\
                        default=0.78e-5")

        self.add_argument("--migration_AF_EU",
                        type=float,
                        dest="m_AF_EU",
                        default=2.5e-5,
                        help="Set African <-> European migration rate ;\
                        default=2.5e-5")

        self.add_argument("--migration_EU_AS",
                        type=float,
                        dest="m_EU_AS",
                        default=3.11e-5,
                        help="Set European <-> Asian migration rate ; \
                        default=3.11e-5")


    # taken from https://stackoverflow.com/questions/5943249 for testing
    def _get_action_from_name(self, name):
        """Given a name, get the Action instance registered with this parser.
        If only it were made available in the ArgumentError object. It is 
        passed as it's first arg...
        """
        container = self._actions
        if name is None:
            return None
        for action in container:
            if '/'.join(action.option_strings) == name:
                return action
            elif action.metavar == name:
                return action
            elif action.dest == name:
                return action

    def error(self, message):
        exc = sys.exc_info()[1]
        if exc:
            exc.argument = self._get_action_from_name(exc.argument_name)
            raise exc
        super(ArgumentParser, self).error(message) 
        
