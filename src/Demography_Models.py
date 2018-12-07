import msprime


class Population(object):
    def __init__(self,
                 name,
                 population_size,
                 growth_rate,
                 samples,
                 generations):
        self.name = name
        self.size = population_size
        self.rate = growth_rate
        self.samples = samples
        self.generations = generations

    def get_sample(self, index):
        return [msprime.Sample(
            population=index,
            time=self.generations)] * self.samples

    def get_configuration(self):
        """get msprime configuration for the population"""
        return msprime.PopulationConfiguration(
            initial_size=self.size,
            growth_rate=self.rate)

    def get_debug_configuration(self, includeRate=False):
        if includeRate:
            return msprime.PopulationConfiguration(
                initial_size=self.size,
                sample_size=2,
                growth_rate=self.rate)

        else:
            return msprime.PopulationConfiguration(
                initial_size=self.size,
                sample_size=2)


class Base_demography(object):
    def __init__(self, S_N1, S_N2, options):
        self.S_N1 = S_N1
        self.S_N2 = S_N2
        self.options = options

        # Set Percent introgression level
        # an instantaneous movement of population from
        # EU_AS to Neand1 (backward in time)
        self.m_PULSE1 = self.options.n1_admix_prop
        self.m_PULSE2 = self.options.n2_admix_prop
        self.set_constants()
        self.set_populations()
        self.set_migrations()
        self.set_demographic_events()

    def simulate(self):
        if(self.options.haplo == "debug"):
            msprime.DemographyDebugger(
                Ne=self.N_A,
                population_configurations=self.get_debug_configuration(),
                migration_matrix=self.get_migration_matrix(),
                demographic_events=self.get_demographic_events()
            ).print_history()

        else:
            replicates = 20 if self.options.haplo == "F4Dstat"\
                            else 1

            return msprime.simulate(
                Ne=self.N_A,
                length=self.length,
                recombination_rate=self.recombination_rate,
                mutation_rate=self.mutation_rate,
                samples=self.get_samples(),
                population_configurations=self.get_population_configuration(),
                migration_matrix=self.get_migration_matrix(),
                demographic_events=self.get_demographic_events(),
                num_replicates=replicates,
                random_seed=self.options.seed
            )

    def set_constants(self):
        self.generation_time = 25  # years/generation
        self.length = self.options.length
        # 2e-8 original from Rajiv, not sure where it is from.
        self.recombination_rate = 1e-8

        # used in Vernot 2014 and 2016
        # matches 0migration and CWZ num_s_star_snps,
        # matches 0migration s_star scores
        # self.mutation_rate = 2.5e-8
        # value best matches CWZ data n_region_ind_snps and s_star scores
        self.mutation_rate = 1.2e-8

        # effective population size used for scaling
        self.N_A = 7310
        # Altai/Vindija lineage effective population size
        self.N_N1 = 1000
        # Eastern Neandertal effective population size
        self.N_N2 = 1000
        # Population size after EUR + ASN join
        self.N_B = 1861
        # Population size of EUR_ASN bottleneck
        # (following initial admixture event)
        self.N_AF0 = 14474
        self.N_EU0 = 1032
        self.N_AS0 = 554
        self.N_CH = 1000
        self.N_DE = 1000

        # Chimp_Human lineage split
        self.T_MH_CH = 7e6 / self.generation_time
        # Neanderthal / modern human split
        self.T_MH_N = 700e3 / self.generation_time
        # Denisovan / Neandertal split time
        self.T_DE_N = 500e3 / self.generation_time
        # Altai/Vindija and Eastern Neandertal popualation split,
        self.T_N1_N2 = ((self.options.t_n1_n2) * 1e3) / self.generation_time
        # African population expansion
        self.T_AF = 148e3 / self.generation_time
        # out of Africa
        self.T_B = 100e3 / self.generation_time
        # Neandertal introgression pulse 1
        self.T_PULSE1 = 55e3 / self.generation_time
        # European / East Asian split
        self.T_EU_AS = 23e3 / self.generation_time
        # Neandertal introgression pulse 2
        self.T_PULSE2 = 18e3 / self.generation_time
        # time of rapid population growth
        self.T_ACL_GRW = 5.1e3 / self.generation_time
        # We need to work out the starting (diploid) population sizes based on
        # the growth rates provided for these two populations
        self.r_EU1 = 0.00307
        self.r_AS1 = 0.0048
        # rapid popualtion growth at 5.1kya
        self.N_EU = 512e3
        self.N_AS = 640e3
        self.N_AF = 424e3

    def set_populations(self):
        """set the self.populations variable
        a list of Populations for the simulation"""
        self.populations = [
            # Altai/Vindija lineage effective population size
            Population(name='N1',
                       population_size=self.N_N1,
                       growth_rate=0.0,
                       samples=self.S_N1,
                       generations=50e3 / self.generation_time),

            # Eastern Neandertal effective population size
            Population(name='N2',
                       population_size=self.N_N2,
                       growth_rate=0.0,
                       samples=self.S_N2,
                       generations=50e3 / self.generation_time),

            Population(name='AF',
                       population_size=self.N_AF,
                       growth_rate=0.0166,
                       samples=self.options.AF_sample_size,
                       generations=0),

            Population(name='EU',
                       population_size=self.N_EU,
                       growth_rate=0.0195,
                       samples=self.options.EU_sample_size,
                       generations=0),

            Population(name='AS',
                       population_size=self.N_AS,
                       growth_rate=0.025,
                       samples=self.options.AS_sample_size,
                       generations=0),

            Population(name='CH',
                       population_size=self.N_CH,
                       growth_rate=0.0,
                       samples=2,
                       generations=0),

            Population(name='DE',
                       population_size=self.N_DE,
                       growth_rate=0.0,
                       samples=2,
                       generations=50e3 / self.generation_time),

        ]

    def set_migrations(self):
        """specify the non-zero migrations in the simulation
        a list of tuples of (name1, name2, value)
        names must match those given in set_populations"""
        self.migrations = [
            ('AF', 'AS', self.options.m_AF_AS),
            ('AF', 'EU', self.options.m_AF_EU),
            ('EU', 'AS', self.options.m_EU_AS)
        ]

    def set_demographic_events(self):
        """set the list of demographic events in self.events"""
        ids = self.get_population_map()

        self.events = []

        self.events += [
            # stop rapid population growth in AF
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                initial_size=self.N_AF0,
                growth_rate=0,
                population_id=ids['AF']),

            # stop rapid population growth in EU
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_EU1,
                population_id=ids['EU']),

            # stop rapid population growth in AS
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_AS1,
                population_id=ids['AS']),

            # set AS popsize to AS0
            msprime.PopulationParametersChange(
                time=self.T_PULSE2-1,
                initial_size=self.N_AS0,
                growth_rate=0,
                population_id=ids['AS']),

            # Neand2 to EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE2,
                source=ids['AS'],
                destination=ids['N2'],
                proportion=self.m_PULSE2),

            # set EU popsize at EU0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_EU0,
                growth_rate=0,
                population_id=ids['EU']),

            # set AS popsize to AS0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_AS0,
                growth_rate=0,
                population_id=ids['AS']),

            # AS merges into EU, now termed "B"
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AS'],
                destination=ids['EU'],
                proportion=1.0)]

        ids['B'] = ids['EU']

        self.events += [
            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=0),

            # migration between "B" and Africa begins
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=self.options.m_AF_B,
                matrix_index=(ids['AF'], ids['B'])),
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=self.options.m_AF_B,
                matrix_index=(ids['B'], ids['AF'])),

            # set parameters of population "B"
            msprime.PopulationParametersChange(
                time=self.T_EU_AS,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # Neand1 to EUR_EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE1,
                source=ids['B'],
                destination=ids['N1'],
                proportion=self.m_PULSE1),

            # Population B merges into Africa at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['B'],
                destination=ids['AF'],
                proportion=1.0),

            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_B,
                rate=0),

            # set parameters of ancestral modern human population
            msprime.PopulationParametersChange(
                time=self.T_AF,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # N_2 merges with N_1 at T_N1_N2
            msprime.MassMigration(
                time=self.T_N1_N2,
                source=ids['N2'],
                destination=ids['N1'],
                proportion=1.0),

            # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(
                time=self.T_N1_N2,
                initial_size=self.N_N1,
                population_id=ids['N1']),

            # DE merges with N1
            msprime.MassMigration(
                time=self.T_DE_N,
                source=ids['DE'],
                destination=ids['N1'],
                proportion=1.0),
            msprime.PopulationParametersChange(
                time=self.T_DE_N,
                initial_size=self.N_N1,
                population_id=ids['N1']),

            # Neandertals merge into modern human lineage at time T_MH_N
            msprime.MassMigration(
                time=self.T_MH_N,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancetral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_N,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # Chimp lineage merges into ancestral hominin
            # population at time T_MH_CH
            msprime.MassMigration(
                time=self.T_MH_CH,
                source=ids['CH'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancestral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_CH,
                initial_size=self.N_A,
                population_id=ids['AF'])
        ]

    def get_population_map(self):
        result = {}
        for i, pop in enumerate(self.populations):
            result[pop.name] = i
        return result

    def get_population_configuration(self):
        """build array of population configurations from self.populations"""
        return [pop.get_configuration()
                for pop in self.populations]

    def get_debug_configuration(self):
        """build array of population configurations for debugging"""
        return [pop.get_debug_configuration()
                for pop in self.populations]

    def get_samples(self):
        result = []
        for i, p in enumerate(self.populations):
            result += p.get_sample(i)
        return result

    def get_migration_matrix(self):
        result = []
        # terrible runtime but I can't think of another general solution
        names = [pop.name for pop in self.populations]
        for n in names:
            row = []
            for m in names:
                found = False
                for entry in self.migrations:
                    if (n == entry[0] and m == entry[1]) or \
                            (n == entry[1] and m == entry[0]):
                        row.append(entry[2])
                        found = True
                        break
                if not found:
                    row.append(0)

            result.append(row)
        return result

    def get_demographic_events(self):
        return sorted(self.events, key=lambda e: e.time)


class Tenn_demography(Base_demography):
    """Small change from base as base was made from Tenn"""
    
    def get_debug_configuration(self):
        return [pop.get_debug_configuration(includeRate=True)
                for pop in self.populations]


class Sriram_demography(Base_demography):
    def set_constants(self):
        Base_demography.set_constants(self)

        # effective population size used for scaling
        self.N_A = 10000
        # Altai/Vindija lineage effective population size
        self.N_N1 = 10000
        # Eastern Neandertal effective population size
        self.N_N2 = 10000
        # Neandertal population bottleneck size
        self.N_N_BN = 100
        # Population size after EUR + ASN join
        self.N_B = 10000
        # Population size of EUR_ASN bottleneck
        # (following initial admixture event)
        self.N_B_BN = 100
        self.N_AF = 10000
        self.N_CH = 10000
        self.N_DE = 10000
        self.N_EU = 10000
        self.N_EU0 = 1000
        self.N_AS0 = 510
        self.N_AS = 10000

        # Neandertal bottleneck            ### DIFFERENT FROM TENNESSEN
        self.T_N_BN = 150e3 / self.generation_time
        # out of Africa                    ### DIFFERENT FROM TENNESSEN
        self.T_B = 75e3 / self.generation_time
        # Neandertal introgression pulse 1 ### DIFFERENT FROM TENNESSEN
        self.T_PULSE1 = 50e3 / self.generation_time
        # A 3rd Bottleneck following shortly after initial admixture event
        # SET AT 1800 GEN
        self.T_B_BN = 45e3 / self.generation_time
        # European / East Asian split
        self.T_EU_AS = 23e3 / self.generation_time

    def set_populations(self):
        """set the self.populations variable
        a list of Populations for the simulation"""

        Base_demography.set_populations(self)

        # only change from tennison is all growth rates are 0
        for pop in self.populations:
            pop.rate = 0.0

    def set_demographic_events(self):
        ids = self.get_population_map()

        self.events = []

        self.events += [
            # Neand2 to EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE2,
                source=ids['AS'],
                destination=ids['N2'],
                proportion=self.m_PULSE2),

            # AS merges into EU, now termed "B"
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AS'],
                destination=ids['EU'],
                proportion=1.0)
        ]

        ids['B'] = ids['EU']

        self.events += [
            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=0),

            # set parameters of population "B"
            msprime.PopulationParametersChange(
                time=self.T_EU_AS,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # population bottleneck of B begins
            msprime.PopulationParametersChange(
                time=self.T_B_BN,
                initial_size=self.N_B_BN,
                growth_rate=0,
                population_id=ids['B']),

            # population bottleneck of B ends shortly before initial admixture
            msprime.PopulationParametersChange(
                time=self.T_B_BN + 20,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # Neand1 to EUR_EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE1,
                source=ids['B'],
                destination=ids['N1'],
                proportion=self.m_PULSE1),

            # Population B merges into Africa at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['B'],
                destination=ids['AF'],
                proportion=1.0),

            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_B,
                rate=0),

            # Neand1 has short bottleneck
            msprime.PopulationParametersChange(
                time=self.T_N_BN,
                initial_size=self.N_N_BN,
                growth_rate=0,
                population_id=ids['N1']),

            # Neand1 bottleneck ends
            msprime.PopulationParametersChange(
                time=self.T_N_BN + 20,
                initial_size=self.N_A,
                growth_rate=0,
                population_id=ids['N1']),

            # N_2 merges with N_1 at T_N1_N2
            msprime.MassMigration(
                time=self.T_N1_N2,
                source=ids['N2'],
                destination=ids['N1'],
                proportion=1.0),

            # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(
                time=self.T_N1_N2,
                initial_size=self.N_A,
                population_id=ids['N1']),

            # DE merges with N1
            msprime.MassMigration(
                time=self.T_DE_N,
                source=ids['DE'],
                destination=ids['N1'],
                proportion=1.0),
            msprime.PopulationParametersChange(
                time=self.T_DE_N,
                initial_size=self.N_A,
                population_id=ids['N1']),

            # Neandertals merge into modern human lineage at time T_MH_N
            msprime.MassMigration(
                time=self.T_MH_N,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancetral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_N,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # Chimp lineage merges into ancestral hominin population
            # at time T_MH_CH
            msprime.MassMigration(
                time=self.T_MH_CH,
                source=ids['CH'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancestral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_CH,
                initial_size=self.N_A,
                population_id=ids['AF'])
        ]


class SplitPop_demography(Base_demography):
    def set_constants(self):
        Base_demography.set_constants(self)

        self.N_SP = 100
        # split population off of B
        self.T_SP = 54e3 / self.generation_time

        # rapid population growth at 5.1kya
        self.r_EU2 = 0.0195
        # rapid population growth at 5.1kya
        self.r_AS2 = 0.025
        # rapid popualtion growth at 5.1kya
        self.r_AF2 = 0.0166

    def set_populations(self):
        self.populations = [
            Population(name='N1',
                       population_size=self.N_N1,
                       growth_rate=0.0,
                       samples=self.S_N1,
                       generations=50e3 / self.generation_time),

            Population(name='N2',
                       population_size=self.N_N2,
                       growth_rate=0.0,
                       samples=self.S_N2,
                       generations=50e3 / self.generation_time),

            Population(name='AF',
                       population_size=self.N_AF,
                       growth_rate=self.r_AF2,
                       samples=self.options.AF_sample_size,
                       generations=0),

            Population(name='EU',
                       population_size=self.N_EU,
                       growth_rate=self.r_EU2,
                       samples=self.options.EU_sample_size,
                       generations=0),

            Population(name='AS',
                       population_size=self.N_AS,
                       growth_rate=self.r_AS2,
                       samples=self.options.AS_sample_size,
                       generations=0),

            Population(name='CH',
                       population_size=self.N_CH,
                       growth_rate=0.0,
                       samples=2,
                       generations=0),

            Population(name='DE',
                       population_size=self.N_DE,
                       growth_rate=0.0,
                       samples=2,
                       generations=50e3 / self.generation_time),

            Population(name='SP',
                       population_size=self.N_SP,
                       growth_rate=0.0,
                       samples=2,
                       generations=50e3 / self.generation_time)

        ]

    def get_samples(self):
        """need to override default to remove last (SP) population"""
        result = []
        for i, p in enumerate(self.populations[:-1]):
            result += p.get_sample(i)
        return result

    def set_demographic_events(self):
        ids = self.get_population_map()

        self.events = []

        self.events += [
            # stop rapid population growth in AF
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                initial_size=self.N_AF0,
                growth_rate=0,
                population_id=ids['AF']),

            # stop rapid population growth in EU
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_EU1,
                population_id=ids['EU']),

            # stop rapid population growth in AS
            msprime.PopulationParametersChange(
                time=self.T_ACL_GRW,
                growth_rate=self.r_AS1,
                population_id=ids['AS']),

            # set AS popsize to AS0
            msprime.PopulationParametersChange(
                time=self.T_PULSE2-1,
                initial_size=self.N_AS0,
                growth_rate=0,
                population_id=ids['AS']),

            # Neand2 to EAS pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE2,
                source=ids['AS'],
                destination=ids['N2'],
                proportion=self.m_PULSE2),

            # set EU popsize at EU0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_EU0,
                growth_rate=0,
                population_id=ids['EU']),

            # set AS popsize to AS0
            msprime.PopulationParametersChange(
                time=self.T_EU_AS-1,
                initial_size=self.N_AS0,
                growth_rate=0,
                population_id=ids['AS']),

            # AS merges into EU, now termed "B"
            msprime.MassMigration(
                time=self.T_EU_AS,
                source=ids['AS'],
                destination=ids['EU'],
                proportion=1.0)
        ]

        ids['B'] = ids['EU']

        self.events += [
            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=0),

            # migration between "B" and Africa begins
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=self.options.m_AF_B,
                matrix_index=(ids['AF'], ids['B'])),
            msprime.MigrationRateChange(
                time=self.T_EU_AS,
                rate=self.options.m_AF_B,
                matrix_index=(ids['B'], ids['AF'])),

            # set parameters of population "B"
            msprime.PopulationParametersChange(
                time=self.T_EU_AS,
                initial_size=self.N_B,
                growth_rate=0,
                population_id=ids['B']),

            # Population splits off of EUR_EAS, 10% migration rate to SplitPop
            msprime.MassMigration(
                time=self.T_SP,
                source=ids['B'],
                destination=ids['SP'],
                proportion=0.1),

            # SplitPop is set at size Ne=100
            msprime.PopulationParametersChange(
                time=self.T_SP,
                initial_size=self.N_SP,
                growth_rate=0,
                population_id=ids['SP']),

            # Neand1 to SplitPop pulse of introgression
            msprime.MassMigration(
                time=self.T_PULSE1,
                source=ids['SP'],
                destination=ids['N1'],
                proportion=self.m_PULSE1),

            # Population B merges into African at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['B'],
                destination=ids['AF'],
                proportion=1.0),

            # SplitPop merges into African at T_B
            msprime.MassMigration(
                time=self.T_B,
                source=ids['SP'],
                destination=ids['AF'],
                proportion=1.0),

            # set all migration rates to zero
            msprime.MigrationRateChange(
                time=self.T_B,
                rate=0),

            # set parameters of ancestral modern human population
            msprime.PopulationParametersChange(
                time=self.T_AF,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # N_2 merges with N_1 at T_N1_N2
            msprime.MassMigration(
                time=self.T_N1_N2,
                source=ids['N2'],
                destination=ids['N1'],
                proportion=1.0),

            # set parameters of ancestral Neandertal population
            msprime.PopulationParametersChange(
                time=self.T_N1_N2,
                initial_size=self.N_N1,
                population_id=ids['N1']),

            # DE merges with N1
            msprime.MassMigration(
                time=self.T_DE_N,
                source=ids['DE'],
                destination=ids['N1'],
                proportion=1.0),
            msprime.PopulationParametersChange(
                time=self.T_DE_N,
                initial_size=self.N_N1,
                population_id=ids['N1']),

            # Neandertals merge into modern human lineage at time T_MH_N
            msprime.MassMigration(
                time=self.T_MH_N,
                source=ids['N1'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancetral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_N,
                initial_size=self.N_A,
                population_id=ids['AF']),

            # Chimp lineage merges into ancestral hominin population
            # at time T_MH_CH
            msprime.MassMigration(
                time=self.T_MH_CH,
                source=ids['CH'],
                destination=ids['AF'],
                proportion=1.0),

            # set parameters of ancestral hominin population
            msprime.PopulationParametersChange(
                time=self.T_MH_CH,
                initial_size=self.N_A,
                population_id=ids['AF'])
        ]

    def get_debug_configuration(self):
        return [pop.get_debug_configuration(includeRate=True)
                for pop in self.populations]
