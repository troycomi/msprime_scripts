from optparse import OptionParser

class admixture_option_parser(OptionParser):
    def __init__(self):
        usage = "usage: python <script> -p nonAfr -o Tenn \
                -s 1 -i 2 -n 0.02 -d 0.0 -t 350 -c F4Dstat \
                -l 1e6 -e 1006 -a 1008"
        Base.__init__(usage=usage)

        self.add_option("-p", "--population", \
                          type = "string", \
                          dest = "pop", \
                          default="nonAfr",\
                          help="call introg_haplotypes in EUR, EAS, \
                                  or all nonAfr; default=nonAfr")

        self.add_option("-o", "--outdir", \
                          type = "string", \
                          dest = "outdir", \
                          default="Tenn", \
                          help="output directory name, \
                              best to include model name; \
                              default=Tenn")

        self.add_option("-s", "--seed", \
                          type = "int",\
                          dest = "seed",\
                          default=1, \
                          help="Set random seed for replicate chromosomes; \
                          default=1")

        self.add_option("-i", "--introgress_pulses",\
                          type = "int",\
                          dest = "pulses",\
                          default=2,\
                          help="Set number of introgression pulses; default=2")

        self.add_option("-n", "--neand1_admixture_proportion",\
                          type = "float",\
                          dest = "n1_admix_prop",\
                          default=0.02,\
                          help="Set N1 admixture proportion; default=0.02")

        self.add_option("-d", "--deni_admixture_propotion",\
                          type = "float",\
                          dest = "n2_admix_prop",\
                          default=0.0,\
                          help="Set N2 admixture proportion; default=0.0")

        self.add_option("-t", "--time_N1_N2_split",\
                          type = "int",\
                          dest = "t_n1_n2",\
                          default=350,\
                          help="Set N1 N2 split time in kya; default=350")

        self.add_option("-c", "--calls_or_stats", \
                          type = "string",\
                          dest = "haplo",\
                          default="F4Dstat",\
                          help="Pick haplotype calls ('haplo') or 'F4Dstat' \
                            or 'vcf' output, or 'debug'; default=F4Dstat")

        self.add_option("-l", "--length_chrom", \
                          type= "float",\
                          dest = "length",\
                          default=1e6,\
                          help="Define length of simulated chromosome; default=1e6")

        self.add_option("-e", "--european_sample_size", \
                          type="int",\
                          dest = "EU_sample_size",\
                          default=1006,\
                          help="Set EU haploid sample size ; default=1006")

        self.add_option("-a", "--asian_sample_size", \
                          type="int",\
                          dest="AS_sample_size",\
                          default=1008,\
                          help="Set AS haploid sample size ; default=1008")

        self.add_option("-r", "--reference_sample_size", \
                          type="int",\
                          dest="AF_sample_size",\
                          default=2,\
                          help="Set AF haploid sample size ; default=2")

        self.add_option("--migration_AF_B", \
                          type="float",\
                          dest="m_AF_B",\
                          default=15e-5,\
                          help="Set African <-> ancestral Eurasian migration rate ; \
                              default=15e-5")

        self.add_option("--migration_AF_AS", \
                          type="float",\
                          dest="m_AF_AS",\
                          default=0.78e-5,\
                          help="Set African <-> Asian migration rate ; default=0.78e-5")

        self.add_option("--migration_AF_EU", \
                          type="float",\
                          dest="m_AF_EU",\
                          default=2.5e-5,\
                          help="Set African <-> European migration rate ; default=2.5e-5")

        self.add_option("--migration_EU_AS", \
                          type="float",\
                          dest="m_EU_AS",\
                          default=3.11e-5,\
                          help="Set European <-> Asian migration rate ; default=3.11e-5")
