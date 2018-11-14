import msprime

# simulate 1Mb chromosomes for faster F4 and D calculation
length = 1e6
# 2e-8 original from Rajiv, not sure where it is from.
recombination_rate = 1e-8
# used in Vernot 2014 and 2016
# matches 0migration and CWZ num_s_star_snps,
# matches 0migration s_star scores
mutation_rate = 2.5e-8
# value best matches CWZ data n_region_ind_snps and s_star scores
mutation_rate = 1.2e-8


def Tenn_demography(
        S_N1, S_N2, S_AF, S_EU, S_AS,
        pulses, seed,
        n1_admix_prop, n2_admix_prop,
        outdir, t_n1_n2,
        haplo, length,
        m_AF_B, m_AF_AS, m_AF_EU, m_EU_AS):
    # effective population size used for scaling
    N_A = 7310
    # Altai/Vindija lineage effective population size
    N_N1 = 1000
    # Eastern Neandertal effective population size
    N_N2 = 1000
    # Population size after EUR + ASN join
    N_B = 1861
    # Population size of EUR_ASN bottleneck
    # (following initial admixture event)
    N_AF0 = 14474
    N_EU0 = 1032
    N_AS0 = 554
    N_CH = 1000
    N_DE = 1000
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    # Chimp_Human lineage split
    T_MH_CH = 7e6 / generation_time
    # Neanderthal / modern human split
    T_MH_N = 700e3 / generation_time
    # Denisovan / Neandertal split time
    T_DE_N = 500e3 / generation_time
    # Altai/Vindija and Eastern Neandertal popualation split,
    # default set to 350kya
    T_N1_N2 = ((t_n1_n2) * 1e3) / generation_time
    # African population expansion
    T_AF = 148e3 / generation_time
    # out of Africa
    T_B = 100e3 / generation_time
    # Neandertal introgression pulse 1
    T_PULSE1 = 55e3 / generation_time
    # European / East Asian split
    T_EU_AS = 23e3 / generation_time
    # Neandertal introgression pulse 2
    T_PULSE2 = 18e3 / generation_time
    # time of rapid population growth
    T_ACL_GRW = 5.1e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU1 = 0.00307
    r_AS1 = 0.0048
    # rapid population growth at 5.1kya
    r_EU2 = 0.0195
    # rapid population growth at 5.1kya
    r_AS2 = 0.025
    # rapid popualtion growth at 5.1kya
    r_AF2 = 0.0166
    N_EU = 512e3
    N_AS = 640e3
    N_AF = 424e3

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
    # an instantaneous movement of population from
    # EU_AS to Neand1 (backward in time)
    m_PULSE1 = n1_admix_prop
    m_PULSE2 = n2_admix_prop

    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            initial_size=N_N1),  # population 0
        msprime.PopulationConfiguration(
            initial_size=N_N2),  # population 1
        msprime.PopulationConfiguration(
            initial_size=N_AF, growth_rate=r_AF2),  # population 2
        msprime.PopulationConfiguration(
            initial_size=N_EU, growth_rate=r_EU2),  # population 3
        msprime.PopulationConfiguration(
            initial_size=N_AS, growth_rate=r_AS2),  # population 4
        msprime.PopulationConfiguration(
            initial_size=N_CH),  # popultaion 5
        msprime.PopulationConfiguration(
            initial_size=N_DE)  # population 6
    ]
    # specify desired sample
    N1_sample = [msprime.Sample(
                    population=0,
                    time=50e3 / generation_time
                   )
                 ] * S_N1
    N2_sample = [msprime.Sample(
                    population=1,
                    time=50e3 / generation_time
                   )
                 ] * S_N2
    AF_sample = [msprime.Sample(
                    population=2,
                    time=0
                   )
                 ] * S_AF
    EU_sample = [msprime.Sample(
                    population=3,
                    time=0
                   )
                 ] * S_EU
    AS_sample = [msprime.Sample(
                    population=4,
                    time=0
                   )
                 ] * S_AS
    CH_sample = [msprime.Sample(
                    population=5,
                    time=0
                   )
                 ] * 2
    DE_sample = [msprime.Sample(
                    population=6,
                    time=50e3 / generation_time
                   )
                 ] * 2
    samples = \
        N1_sample + \
        N2_sample + \
        AF_sample + \
        EU_sample + \
        AS_sample + \
        CH_sample + \
        DE_sample

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
        # stop rapid population growth in AF
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            initial_size=N_AF0,
            growth_rate=0,
            population_id=2),
        # stop rapid population growth in EU
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            growth_rate=r_EU1,
            population_id=3),
        # stop rapid population growth in AS
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            growth_rate=r_AS1,
            population_id=4),
        # set AS popsize to AS0
        msprime.PopulationParametersChange(
            time=T_PULSE2-1,
            initial_size=N_AS0,
            growth_rate=0,
            population_id=4),
        # set EU popsize at EU0
        msprime.PopulationParametersChange(
            time=T_EU_AS-1,
            initial_size=N_EU0,
            growth_rate=0,
            population_id=3),
        # set AS popsize to AS0
        msprime.PopulationParametersChange(
            time=T_EU_AS-1,
            initial_size=N_AS0,
            growth_rate=0,
            population_id=4),
        # AS merges into EU, now termed "B"
        msprime.MassMigration(
            time=T_EU_AS,
            source=4,
            destination=3,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=0),
        # migration between "B" and Africa begins
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=m_AF_B,
            matrix_index=(2, 3)),
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=m_AF_B,
            matrix_index=(3, 2)),
        # set parameters of population "B"
        msprime.PopulationParametersChange(
            time=T_EU_AS,
            initial_size=N_B,
            growth_rate=0,
            population_id=3),
        # Neand1 to EUR_EAS pulse of introgression
        msprime.MassMigration(
            time=T_PULSE1,
            source=3,
            destination=0,
            proportion=m_PULSE1),
        # Population B merges into Africa at T_B
        msprime.MassMigration(
            time=T_B,
            source=3,
            destination=2,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_B,
            rate=0),
        # set parameters of ancestral modern human population
        msprime.PopulationParametersChange(
            time=T_AF,
            initial_size=N_A,
            population_id=2),
        # N_2 merges with N_1 at T_N1_N2
        msprime.MassMigration(
            time=T_N1_N2,
            source=1,
            destination=0,
            proportion=1.0),
        # set parameters of ancestral Neandertal population
        msprime.PopulationParametersChange(
            time=T_N1_N2,
            initial_size=N_N1,
            population_id=0),
        # DE merges with N1
        msprime.MassMigration(
            time=T_DE_N,
            source=6,
            destination=0,
            proportion=1.0),
        msprime.PopulationParametersChange(
            time=T_DE_N,
            initial_size=N_N1,
            population_id=0),
        # Neandertals merge into modern human lineage at time T_MH_N
        msprime.MassMigration(
            time=T_MH_N,
            source=0,
            destination=2,
            proportion=1.0),
        # set parameters of ancetral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_N,
            initial_size=N_A,
            population_id=2),
        # Chimp lineage merges into ancestral hominin population
        # at time T_MH_Ch
        msprime.MassMigration(
            time=T_MH_CH,
            source=5,
            destination=2,
            proportion=1.0),
        # set parameters of ancestral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_CH,
            initial_size=N_A,
            population_id=2)
        ]
    two_pulse = [
        # stop rapid population growth in AF
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            initial_size=N_AF0,
            growth_rate=0,
            population_id=2),
        # stop rapid population growth in EU
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            growth_rate=r_EU1,
            population_id=3),
        # stop rapid population growth in AS
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            growth_rate=r_AS1,
            population_id=4),
        # set AS popsize to AS0
        msprime.PopulationParametersChange(
            time=T_PULSE2-1,
            initial_size=N_AS0,
            growth_rate=0,
            population_id=4),
        # Neand2 to EAS pulse of introgression
        msprime.MassMigration(
            time=T_PULSE2,
            source=4,
            destination=1,
            proportion=m_PULSE2),
        # set EU popsize at EU0
        msprime.PopulationParametersChange(
            time=T_EU_AS-1,
            initial_size=N_EU0,
            growth_rate=0,
            population_id=3),
        # set AS popsize to AS0
        msprime.PopulationParametersChange(
            time=T_EU_AS-1,
            initial_size=N_AS0,
            growth_rate=0,
            population_id=4),
        # AS merges into EU, now termed "B"
        msprime.MassMigration(
            time=T_EU_AS,
            source=4,
            destination=3,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=0),
        # migration between "B" and Africa begins
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=m_AF_B,
            matrix_index=(2, 3)),
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=m_AF_B,
            matrix_index=(3, 2)),
        # set parameters of population "B"
        msprime.PopulationParametersChange(
            time=T_EU_AS,
            initial_size=N_B,
            growth_rate=0,
            population_id=3),
        # Neand1 to EUR_EAS pulse of introgression
        msprime.MassMigration(
            time=T_PULSE1,
            source=3,
            destination=0,
            proportion=m_PULSE1),
        # Population B merges into Africa at T_B
        msprime.MassMigration(
            time=T_B,
            source=3,
            destination=2,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_B,
            rate=0),
        # set parameters of ancestral modern human population
        msprime.PopulationParametersChange(
            time=T_AF,
            initial_size=N_A,
            population_id=2),
        # N_2 merges with N_1 at T_N1_N2
        msprime.MassMigration(
            time=T_N1_N2,
            source=1,
            destination=0,
            proportion=1.0),
        # set parameters of ancestral Neandertal population
        msprime.PopulationParametersChange(
            time=T_N1_N2,
            initial_size=N_N1,
            population_id=0),
        # DE merges with N1
        msprime.MassMigration(
            time=T_DE_N,
            source=6,
            destination=0,
            proportion=1.0),
        msprime.PopulationParametersChange(
            time=T_DE_N,
            initial_size=N_N1,
            population_id=0),
        # Neandertals merge into modern human lineage at time T_MH_N
        msprime.MassMigration(
            time=T_MH_N,
            source=0,
            destination=2,
            proportion=1.0),
        # set parameters of ancetral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_N,
            initial_size=N_A,
            population_id=2),
        # Chimp lineage merges into ancestral hominin
        # population at time T_MH_CH
        msprime.MassMigration(
            time=T_MH_CH,
            source=5,
            destination=2,
            proportion=1.0),
        # set parameters of ancestral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_CH,
            initial_size=N_A,
            population_id=2)
    ]

    if (pulses == 1):
        demographic_events = one_pulse
    elif (pulses == 2):
        demographic_events = two_pulse

    # DEBUG OPTION
    if (haplo == "debug"):
        dp = msprime.DemographyDebugger(
            Ne=N_A,
            population_configurations=[
                # population 0
                msprime.PopulationConfiguration(
                    initial_size=N_N1,
                    sample_size=2),
                # population 1
                msprime.PopulationConfiguration(
                    initial_size=N_N2,
                    sample_size=2),
                # population 2
                msprime.PopulationConfiguration(
                    initial_size=N_AF,
                    growth_rate=r_AF2,
                    sample_size=2),
                # population 3
                msprime.PopulationConfiguration(
                    initial_size=N_EU,
                    growth_rate=r_EU2,
                    sample_size=2),
                # population 4
                msprime.PopulationConfiguration(
                    initial_size=N_AS,
                    growth_rate=r_AS2,
                    sample_size=2),
                # popultaion 5
                msprime.PopulationConfiguration(
                    initial_size=N_CH,
                    sample_size=2),
                # population 6
                msprime.PopulationConfiguration(
                    initial_size=N_DE,
                    sample_size=2)
                ],
            migration_matrix=migration_matrix,
            demographic_events=demographic_events
        )
        dp.print_history()

    # HAPLOTYPE SIMULATION FOR F4 and DSTAT CALCULATIONS
    elif (haplo == "F4Dstat"):
        return msprime.simulate(
                    Ne=N_A,
                    length=length,
                    recombination_rate=recombination_rate,
                    mutation_rate=mutation_rate,
                    samples=samples,
                    population_configurations=population_configurations,
                    migration_matrix=migration_matrix,
                    demographic_events=demographic_events,
                    num_replicates=20,
                    random_seed=seed
                )

    # HAPLOTYPE SIMULATION FOR DESERT DISTRIBUTION
    elif (haplo == "haplo" or haplo == "vcf"):
        # Simulate 10Mb chromosomes for looking at desert distribution
        return msprime.simulate(
                    Ne=N_A,
                    length=length,
                    recombination_rate=recombination_rate,
                    mutation_rate=mutation_rate,
                    samples=samples,
                    population_configurations=population_configurations,
                    migration_matrix=migration_matrix,
                    demographic_events=demographic_events,
                    num_replicates=1,
                    random_seed=seed
                )


def Sriram_demography(
        S_N1, S_N2, S_AF, S_EU, S_AS,
        pulses, seed,
        n1_admix_prop, n2_admix_prop,
        outdir, t_n1_n2,
        haplo, length,
        m_AF_B, m_AF_AS, m_AF_EU, m_EU_AS):
    # effective population size used for scaling
    N_A = 10000
    # Altai/Vindija lineage effective population size
    N_N1 = 10000
    # Eastern Neandertal effective population size
    N_N2 = 10000
    # Neandertal population bottleneck size
    N_N_BN = 100
    # Population size after EUR + ASN join
    N_B = 10000
    # Population size of EUR_ASN bottleneck (following initial admixture event)
    N_B_BN = 100
    N_AF = 10000
    N_CH = 10000
    N_DE = 10000
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    # Chimp_Human lineage split (7 MYA)
    T_MH_CH = 7e6 / generation_time
    # Neanderthal / modern human split (700 KYA)
    T_MH_N = 700e3 / generation_time
    # Denisovan / Neandertal split time
    T_DE_N = 500e3 / generation_time
    # Altai/Vindija and Eastern Neandertal popualation split,
    # default set to 350kya
    T_N1_N2 = ((t_n1_n2) * 1e3) / generation_time
    # Neandertal bottleneck            ### DIFFERENT FROM TENNESSEN
    T_N_BN = 150e3 / generation_time
    # out of Africa                    ### DIFFERENT FROM TENNESSEN
    T_B = 75e3 / generation_time
    # Neandertal introgression pulse 1 ### DIFFERENT FROM TENNESSEN
    T_PULSE1 = 50e3 / generation_time
    # A 3rd Bottleneck following shortly after initial admixture event
    # SET AT 1800 GEN
    T_B_BN = 45e3 / generation_time
    # European / East Asian split
    T_EU_AS = 23e3 / generation_time
    # Neandertal introgression pulse 2
    T_PULSE2 = 18e3 / generation_time

    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    N_EU = 10000
    N_AS = 10000

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
    # an instantaneous movement of population
    # from EU_AS to Neand1 (backward in time)
    m_PULSE1 = n1_admix_prop
    m_PULSE2 = n2_admix_prop

    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        # population 0
        msprime.PopulationConfiguration(
            initial_size=N_N1),
        # population 1
        msprime.PopulationConfiguration(
            initial_size=N_N2),
        # population 2
        msprime.PopulationConfiguration(
            initial_size=N_AF),
        # population 3
        msprime.PopulationConfiguration(
            initial_size=N_EU),
        # population 4
        msprime.PopulationConfiguration(
            initial_size=N_AS),
        # popultaion 5
        msprime.PopulationConfiguration(
            initial_size=N_CH),
        # population 6
        msprime.PopulationConfiguration(
            initial_size=N_DE)
    ]
    # specify desired sample
    N1_sample = [msprime.Sample(population=0,
                                time=50e3 / generation_time)] * S_N1
    N2_sample = [msprime.Sample(population=1,
                                time=50e3 / generation_time)] * S_N2
    AF_sample = [msprime.Sample(population=2,
                                time=0)] * S_AF
    EU_sample = [msprime.Sample(population=3,
                                time=0)] * S_EU
    AS_sample = [msprime.Sample(population=4,
                                time=0)] * S_AS
    CH_sample = [msprime.Sample(population=5,
                                time=0)] * 2
    DE_sample = [msprime.Sample(population=6,
                                time=50e3 / generation_time)] * 2
    samples = \
        N1_sample + \
        N2_sample + \
        AF_sample + \
        EU_sample + \
        AS_sample + \
        CH_sample + \
        DE_sample

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
        # AS merges into EU, now termed "B"
        msprime.MassMigration(
            time=T_EU_AS,
            source=4,
            destination=3,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=0),
        # set parameters of population "B"
        msprime.PopulationParametersChange(
            time=T_EU_AS,
            initial_size=N_B,
            growth_rate=0,
            population_id=3),
        # population bottleneck of B begins
        msprime.PopulationParametersChange(
            time=T_B_BN,
            initial_size=N_B_BN,
            growth_rate=0,
            population_id=3),
        # population bottleneck ends shortly before initial admixture
        msprime.PopulationParametersChange(
            time=T_B_BN + 20,
            initial_size=N_B,
            growth_rate=0,
            population_id=3),
        # Neand1 to EUR_EAS pulse of introgression
        msprime.MassMigration(
            time=T_PULSE1,
            source=3,
            destination=0,
            proportion=m_PULSE1),
        # Population B merges into Africa at T_B
        msprime.MassMigration(
            time=T_B,
            source=3,
            destination=2,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_B,
            rate=0),
        # Neand1 has short bottleneck
        msprime.PopulationParametersChange(
            time=T_N_BN,
            initial_size=N_N_BN,
            growth_rate=0,
            population_id=0),
        # Neand1 bottleneck ends
        msprime.PopulationParametersChange(
            time=T_N_BN + 20,
            initial_size=N_N1,
            growth_rate=0,
            population_id=0),
        # N_2 merges with N_1 at T_N1_N2
        msprime.MassMigration(
            time=T_N1_N2,
            source=1,
            destination=0,
            proportion=1.0),
        # set parameters of ancestral Neandertal population
        msprime.PopulationParametersChange(
            time=T_N1_N2,
            initial_size=N_N1,
            population_id=0),
        # DE merges with N1
        msprime.MassMigration(
            time=T_DE_N,
            source=6,
            destination=0,
            proportion=1.0),
        msprime.PopulationParametersChange(
            time=T_DE_N,
            initial_size=N_N1,
            population_id=0),
        # Neandertals merge into modern human lineage at time T_MH_N
        msprime.MassMigration(
            time=T_MH_N,
            source=0,
            destination=2,
            proportion=1.0),
        # set parameters of ancetral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_N,
            initial_size=N_A,
            population_id=2),
        # Chimp lineage merges into ancestral hominin population
        # at time T_MH_Ch
        msprime.MassMigration(
            time=T_MH_CH,
            source=5,
            destination=2,
            proportion=1.0),
        # set parameters of ancestral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_CH,
            initial_size=N_A,
            population_id=2)
    ]
    two_pulse = [
        # Neand2 to EAS pulse of introgression
        msprime.MassMigration(
            time=T_PULSE2,
            source=4,
            destination=1,
            proportion=m_PULSE2),
        # AS merges into EU, now termed "B"
        msprime.MassMigration(
            time=T_EU_AS,
            source=4,
            destination=3,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=0),
        # set parameters of population "B"
        msprime.PopulationParametersChange(
            time=T_EU_AS,
            initial_size=N_B,
            growth_rate=0,
            population_id=3),
        # population bottleneck of B begins
        msprime.PopulationParametersChange(
            time=T_B_BN,
            initial_size=N_B_BN,
            growth_rate=0,
            population_id=3),
        # population bottleneck of B ends shortly before initial admixture
        msprime.PopulationParametersChange(
            time=T_B_BN + 20,
            initial_size=N_B,
            growth_rate=0,
            population_id=3),
        # Neand1 to EUR_EAS pulse of introgression
        msprime.MassMigration(
            time=T_PULSE1,
            source=3,
            destination=0,
            proportion=m_PULSE1),
        # Population B merges into Africa at T_B
        msprime.MassMigration(
            time=T_B,
            source=3,
            destination=2,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_B,
            rate=0),
        # Neand1 has short bottleneck
        msprime.PopulationParametersChange(
            time=T_N_BN,
            initial_size=N_N_BN,
            growth_rate=0,
            population_id=0),
        # Neand1 bottleneck ends
        msprime.PopulationParametersChange(
            time=T_N_BN + 20,
            initial_size=N_A,
            growth_rate=0,
            population_id=0),
        # N_2 merges with N_1 at T_N1_N2
        msprime.MassMigration(
            time=T_N1_N2,
            source=1,
            destination=0,
            proportion=1.0),
        # set parameters of ancestral Neandertal population
        msprime.PopulationParametersChange(
            time=T_N1_N2,
            initial_size=N_A,
            population_id=0),
        # DE merges with N1
        msprime.MassMigration(
            time=T_DE_N,
            source=6,
            destination=0,
            proportion=1.0),
        msprime.PopulationParametersChange(
            time=T_DE_N,
            initial_size=N_A,
            population_id=0),
        # Neandertals merge into modern human lineage at time T_MH_N
        msprime.MassMigration(
            time=T_MH_N,
            source=0,
            destination=2,
            proportion=1.0),
        # set parameters of ancetral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_N,
            initial_size=N_A,
            population_id=2),
        # Chimp lineage merges into ancestral hominin population
        # at time T_MH_CH
        msprime.MassMigration(
            time=T_MH_CH,
            source=5,
            destination=2,
            proportion=1.0),
        # set parameters of ancestral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_CH,
            initial_size=N_A,
            population_id=2)
    ]

    if (pulses == 1):
        demographic_events = one_pulse
    elif (pulses == 2):
        demographic_events = two_pulse

    # DEBUG OPTION
    if (haplo == "debug"):
        dp = msprime.DemographyDebugger(
            Ne=N_A,
            population_configurations=[
                # population 0
                msprime.PopulationConfiguration(
                    initial_size=N_N1,
                    sample_size=2),
                # population 1
                msprime.PopulationConfiguration(
                    initial_size=N_N2,
                    sample_size=2),
                # population 2
                msprime.PopulationConfiguration(
                    initial_size=N_AF,
                    sample_size=2),
                # population 3
                msprime.PopulationConfiguration(
                    initial_size=N_EU,
                    sample_size=2),
                # population 4
                msprime.PopulationConfiguration(
                    initial_size=N_AS,
                    sample_size=2),
                # popultaion 5
                msprime.PopulationConfiguration(
                    initial_size=N_CH,
                    sample_size=2),
                # population 6
                msprime.PopulationConfiguration(
                    initial_size=N_DE,
                    sample_size=2)
            ],
            migration_matrix=migration_matrix,
            demographic_events=demographic_events
        )
        dp.print_history()

    # HAPLOTYPE SIMULATION FOR F4 and DSTAT CALCULATIONS
    elif (haplo == "F4Dstat"):
        return msprime.simulate(
                    Ne=N_A,
                    length=length,
                    recombination_rate=recombination_rate,
                    mutation_rate=mutation_rate,
                    samples=samples,
                    population_configurations=population_configurations,
                    migration_matrix=migration_matrix,
                    demographic_events=demographic_events,
                    num_replicates=20,
                    random_seed=seed
                )

    # HAPLOTYPE SIMULATION FOR DESERT DISTRIBUTION
    elif (haplo == "haplo" or haplo == "vcf"):
        # Simulate 10Mb chromosomes for looking at desert distribution
        return msprime.simulate(
                Ne=N_A,
                length=length,
                recombination_rate=recombination_rate,
                mutation_rate=mutation_rate,
                samples=samples,
                population_configurations=population_configurations,
                migration_matrix=migration_matrix,
                demographic_events=demographic_events,
                num_replicates=1,
                random_seed=seed
                )


def SplitPop_demography(
        S_N1, S_N2, S_AF, S_EU, S_AS,
        pulses, seed,
        n1_admix_prop, n2_admix_prop,
        outdir, t_n1_n2,
        haplo, length,
        m_AF_B, m_AF_AS, m_AF_EU, m_EU_AS):
    # effective population size used for scaling
    N_A = 7310
    # Altai/Vindija lineage effective population size
    N_N1 = 1000
    # Eastern Neandertal effective population size
    N_N2 = 1000
    # Population size after EUR + ASN join
    N_B = 1861
    N_AF0 = 14474
    N_EU0 = 1032
    N_AS0 = 554
    N_CH = 1000
    N_DE = 1000
    N_SP = 100
    # Times are provided in years, so we convert into generations.
    generation_time = 25
    # Chimp_Human lineage split
    T_MH_CH = 7e6 / generation_time
    # Neanderthal / modern human split
    T_MH_N = 700e3 / generation_time
    # Denisovan / Neandertal split time
    T_DE_N = 500e3 / generation_time
    # Altai/Vindija and Eastern Neandertal popualation split,
    # default set to 350kya
    T_N1_N2 = ((t_n1_n2) * 1e3) / generation_time
    # African population expansion
    T_AF = 148e3 / generation_time
    # out of Africa
    T_B = 100e3 / generation_time
    # Neandertal introgression pulse 1
    T_PULSE1 = 55e3 / generation_time
    # split population off of B
    T_SP = 54e3 / generation_time
    # European / East Asian split
    T_EU_AS = 23e3 / generation_time
    # Neandertal introgression pulse 2
    T_PULSE2 = 18e3 / generation_time
    # time of rapid population growth
    T_ACL_GRW = 5.1e3 / generation_time
    # We need to work out the starting (diploid) population sizes based on
    # the growth rates provided for these two populations
    r_EU1 = 0.00307
    r_AS1 = 0.0048
    # rapid population growth at 5.1kya
    r_EU2 = 0.0195
    # rapid population growth at 5.1kya
    r_AS2 = 0.025
    # rapid popualtion growth at 5.1kya
    r_AF2 = 0.0166
    N_EU = 512e3
    N_AS = 640e3
    N_AF = 424e3

    # Migration rates during the various epochs.
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
    m_SP_N1 = 0
    m_SP_N2 = 0
    m_SP_AF = 0
    m_SP_EU = 0
    m_SP_AS = 0
    m_SP_CH = 0
    m_SP_DE = 0

    # Set Percent introgression level
    # an instantaneous movement of population
    # from EU_AS to Neand1 (backward in time)
    m_PULSE1 = n1_admix_prop
    m_PULSE2 = n2_admix_prop

    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=YRI, 1=CEU and 2=CHB
    # initially.
    population_configurations = [
        # population 0
        msprime.PopulationConfiguration(
            initial_size=N_N1),
        # population 1
        msprime.PopulationConfiguration(
            initial_size=N_N2),
        # population 2
        msprime.PopulationConfiguration(
            initial_size=N_AF,
            growth_rate=r_AF2),
        # population 3
        msprime.PopulationConfiguration(
            initial_size=N_EU,
            growth_rate=r_EU2),
        # population 4
        msprime.PopulationConfiguration(
            initial_size=N_AS,
            growth_rate=r_AS2),
        # popultaion 5
        msprime.PopulationConfiguration(
            initial_size=N_CH),
        # population 6
        msprime.PopulationConfiguration(
            initial_size=N_DE),
        # population 7
        msprime.PopulationConfiguration(
            initial_size=N_SP)
        ]

    # specify desired sample
    N1_sample = [msprime.Sample(population=0,
                                time=50e3 / generation_time)] * S_N1
    N2_sample = [msprime.Sample(population=1,
                                time=50e3 / generation_time)] * S_N2
    AF_sample = [msprime.Sample(population=2,
                                time=0)] * S_AF
    EU_sample = [msprime.Sample(population=3,
                                time=0)] * S_EU
    AS_sample = [msprime.Sample(population=4,
                                time=0)] * S_AS
    CH_sample = [msprime.Sample(population=5,
                                time=0)] * 2
    DE_sample = [msprime.Sample(population=6,
                                time=50e3 / generation_time)] * 2
    samples = \
        N1_sample + \
        N2_sample + \
        AF_sample + \
        EU_sample + \
        AS_sample + \
        CH_sample + \
        DE_sample

    # set up the initial migration matrix
    migration_matrix = [
        [      0, m_N1_N2, m_N1_AF, m_N1_EU, m_N1_AS, m_N1_CH, m_DE_N1, m_SP_N1],
        [m_N1_N2,       0, m_N2_AF, m_N2_EU, m_N2_AS, m_N2_CH, m_DE_N2, m_SP_N2],
        [m_N1_AF, m_N2_AF,       0, m_AF_EU, m_AF_AS, m_AF_CH, m_DE_AF, m_SP_AF],
        [m_N1_EU, m_N2_EU, m_AF_EU,       0, m_EU_AS, m_EU_CH, m_DE_EU, m_SP_EU],
        [m_N1_AS, m_N2_AS, m_AF_AS, m_EU_AS,       0, m_AS_CH, m_DE_AS, m_SP_AS],
        [m_N1_CH, m_N2_CH, m_AF_CH, m_EU_CH, m_AS_CH,       0, m_DE_CH, m_SP_CH],
        [m_DE_N1, m_DE_N2, m_DE_AF, m_DE_EU, m_DE_AS, m_DE_CH,       0, m_SP_DE],
        [m_SP_N1, m_SP_N2, m_SP_AF, m_SP_EU, m_SP_AS, m_SP_CH, m_SP_DE,       0]
    ]
    one_pulse = [
        # stop rapid population growth in AF
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            initial_size=N_AF0,
            growth_rate=0,
            population_id=2),
        # stop rapid population growth in EU
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            growth_rate=r_EU1,
            population_id=3),
        # stop rapid population growth in AS
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            growth_rate=r_AS1,
            population_id=4),
        # set AS popsize to AS0
        msprime.PopulationParametersChange(
            time=T_PULSE2-1,
            initial_size=N_AS0,
            growth_rate=0,
            population_id=4),
        # set EU popsize at EU0
        msprime.PopulationParametersChange(
            time=T_EU_AS-1,
            initial_size=N_EU0,
            growth_rate=0,
            population_id=3),
        # set AS popsize to AS0
        msprime.PopulationParametersChange(
            time=T_EU_AS-1,
            initial_size=N_AS0,
            growth_rate=0,
            population_id=4),
        # AS merges into EU, now termed "B"
        msprime.MassMigration(
            time=T_EU_AS,
            source=4,
            destination=3,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=0),
        # migration between "B" and Africa begins
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=m_AF_B,
            matrix_index=(2, 3)),
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=m_AF_B,
            matrix_index=(3, 2)),
        # set parameters of population "B"
        msprime.PopulationParametersChange(
            time=T_EU_AS,
            initial_size=N_B,
            growth_rate=0,
            population_id=3),
        # Population splits off of EUR_EAS, 10% migration rate to SplitPop
        msprime.MassMigration(
            time=T_SP,
            source=3,
            destination=7,
            proportion=0.1),
        # SplitPop is set at size Ne=100
        msprime.PopulationParametersChange(
            time=T_SP,
            initial_size=N_SP,
            growth_rate=0,
            population_id=7),
        # Neand1 to SplitPop pulse of introgression
        msprime.MassMigration(
            time=T_PULSE1,
            source=7,
            destination=0,
            proportion=m_PULSE1),
        # Population B merges into African at T_B
        msprime.MassMigration(
            time=T_B,
            source=3,
            destination=2,
            proportion=1.0),
        # SplitPop merges into African at T_B
        msprime.MassMigration(
            time=T_B,
            source=7,
            destination=2,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_B,
            rate=0),
        # set parameters of ancestral modern human population
        msprime.PopulationParametersChange(
            time=T_AF,
            initial_size=N_A,
            population_id=2),
        # N_2 merges with N_1 at T_N1_N2
        msprime.MassMigration(
            time=T_N1_N2,
            source=1,
            destination=0,
            proportion=1.0),
        # set parameters of ancestral Neandertal population
        msprime.PopulationParametersChange(
            time=T_N1_N2,
            initial_size=N_N1,
            population_id=0),
        # DE merges with N1
        msprime.MassMigration(
            time=T_DE_N,
            source=6,
            destination=0,
            proportion=1.0),
        msprime.PopulationParametersChange(
            time=T_DE_N,
            initial_size=N_N1,
            population_id=0),
        # Neandertals merge into modern human lineage at time T_MH_N
        msprime.MassMigration(
            time=T_MH_N,
            source=0,
            destination=2,
            proportion=1.0),
        # set parameters of ancetral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_N,
            initial_size=N_A,
            population_id=2),
        # Chimp lineage merges into ancestral hominin population
        # at time T_MH_Ch
        msprime.MassMigration(
            time=T_MH_CH,
            source=5,
            destination=2,
            proportion=1.0),
        # set parameters of ancestral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_CH,
            initial_size=N_A,
            population_id=2)
    ]

    two_pulse = [
        # stop rapid population growth in AF
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            initial_size=N_AF0,
            growth_rate=0,
            population_id=2),
        # stop rapid population growth in EU
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            growth_rate=r_EU1,
            population_id=3),
        # stop rapid population growth in AS
        msprime.PopulationParametersChange(
            time=T_ACL_GRW,
            growth_rate=r_AS1,
            population_id=4),
        # set AS popsize to AS0
        msprime.PopulationParametersChange(
            time=T_PULSE2-1,
            initial_size=N_AS0,
            growth_rate=0,
            population_id=4),
        # Neand2 to EAS pulse of introgression
        msprime.MassMigration(
            time=T_PULSE2,
            source=4,
            destination=1,
            proportion=m_PULSE2),
        # set EU popsize at EU0
        msprime.PopulationParametersChange(
            time=T_EU_AS-1,
            initial_size=N_EU0,
            growth_rate=0,
            population_id=3),
        # set AS popsize to AS0
        msprime.PopulationParametersChange(
            time=T_EU_AS-1,
            initial_size=N_AS0,
            growth_rate=0,
            population_id=4),
        # AS merges into EU, now termed "B"
        msprime.MassMigration(
            time=T_EU_AS,
            source=4,
            destination=3,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=0),
        # migration between "B" and Africa begins
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=m_AF_B,
            matrix_index=(2, 3)),
        msprime.MigrationRateChange(
            time=T_EU_AS,
            rate=m_AF_B,
            matrix_index=(3, 2)),
        # set parameters of population "B"
        msprime.PopulationParametersChange(
            time=T_EU_AS,
            initial_size=N_B,
            growth_rate=0,
            population_id=3),
        # Population splits off of EUR_EAS, 10% migration rate to SplitPop
        msprime.MassMigration(
            time=T_SP,
            source=3,
            destination=7,
            proportion=0.1),
        # SplitPop is set at size Ne=100
        msprime.PopulationParametersChange(
            time=T_SP,
            initial_size=N_SP,
            growth_rate=0,
            population_id=7),
        # Neand1 to SplitPop pulse of introgression
        msprime.MassMigration(
            time=T_PULSE1,
            source=7,
            destination=0,
            proportion=m_PULSE1),
        # Population B merges into African at T_B
        msprime.MassMigration(
            time=T_B,
            source=3,
            destination=2,
            proportion=1.0),
        # SplitPop merges into African at T_B
        msprime.MassMigration(
            time=T_B,
            source=7,
            destination=2,
            proportion=1.0),
        # set all migration rates to zero
        msprime.MigrationRateChange(
            time=T_B,
            rate=0),
        # set parameters of ancestral modern human population
        msprime.PopulationParametersChange(
            time=T_AF,
            initial_size=N_A,
            population_id=2),
        # N_2 merges with N_1 at T_N1_N2
        msprime.MassMigration(
            time=T_N1_N2,
            source=1,
            destination=0,
            proportion=1.0),
        # set parameters of ancestral Neandertal population
        msprime.PopulationParametersChange(
            time=T_N1_N2,
            initial_size=N_N1,
            population_id=0),
        # DE merges with N1
        msprime.MassMigration(
            time=T_DE_N,
            source=6,
            destination=0,
            proportion=1.0),
        msprime.PopulationParametersChange(
            time=T_DE_N,
            initial_size=N_N1,
            population_id=0),
        # Neandertals merge into modern human lineage at time T_MH_N
        msprime.MassMigration(
            time=T_MH_N,
            source=0,
            destination=2,
            proportion=1.0),
        # set parameters of ancetral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_N,
            initial_size=N_A,
            population_id=2),
        # Chimp lineage merges into ancestral hominin population
        # at time T_MH_CH
        msprime.MassMigration(
            time=T_MH_CH,
            source=5,
            destination=2,
            proportion=1.0),
        # set parameters of ancestral hominin population
        msprime.PopulationParametersChange(
            time=T_MH_CH,
            initial_size=N_A,
            population_id=2)
    ]

    if (pulses == 1):
        demographic_events = one_pulse
    elif (pulses == 2):
        demographic_events = two_pulse

    # DEBUG OPTION
    if (haplo == "debug"):
        dp = msprime.DemographyDebugger(
                Ne=N_A,
                population_configurations=[
                    # population 0
                    msprime.PopulationConfiguration(
                        initial_size=N_N1,
                        sample_size=2),
                    # population 1
                    msprime.PopulationConfiguration(
                        initial_size=N_N2,
                        sample_size=2),
                    # population 2
                    msprime.PopulationConfiguration(
                        initial_size=N_AF,
                        growth_rate=r_AF2,
                        sample_size=2),
                    # population 3
                    msprime.PopulationConfiguration(
                        initial_size=N_EU,
                        growth_rate=r_EU2,
                        sample_size=2),
                    # population 4
                    msprime.PopulationConfiguration(
                        initial_size=N_AS,
                        growth_rate=r_AS2,
                        sample_size=2),
                    # popultaion 5
                    msprime.PopulationConfiguration(
                        initial_size=N_CH,
                        sample_size=2),
                    # population 6
                    msprime.PopulationConfiguration(
                        initial_size=N_DE,
                        sample_size=2),
                    # population 7
                    msprime.PopulationConfiguration(
                        initial_size=N_SP,
                        sample_size=2)
                ],
                migration_matrix=migration_matrix,
                demographic_events=demographic_events
            )
        dp.print_history()

    # HAPLOTYPE SIMULATION FOR F4 and DSTAT CALCULATIONS
    elif (haplo == "F4Dstat"):
        return msprime.simulate(
                    Ne=N_A,
                    length=length,
                    recombination_rate=recombination_rate,
                    mutation_rate=mutation_rate,
                    samples=samples,
                    population_configurations=population_configurations,
                    migration_matrix=migration_matrix,
                    demographic_events=demographic_events,
                    num_replicates=20,
                    random_seed=seed
                )

    # HAPLOTYPE SIMULATION FOR DESERT DISTRIBUTION
    elif (haplo == "haplo" or haplo == "vcf"):
        # Simulate 10Mb chromosomes for looking at desert distribution
        return msprime.simulate(
                    Ne=N_A,
                    length=length,
                    recombination_rate=recombination_rate,
                    mutation_rate=mutation_rate,
                    samples=samples,
                    population_configurations=population_configurations,
                    migration_matrix=migration_matrix,
                    demographic_events=demographic_events,
                    num_replicates=1,
                    random_seed=seed
                )
