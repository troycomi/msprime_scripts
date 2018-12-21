import pytest
import Option_Parser
import Demography_Models
import msprime


@pytest.fixture
def default():
    parser = Option_Parser.admixture_option_parser()
    options = parser.parse_args([])
    return Demography_Models.Base_demography(options)


@pytest.fixture
def nonDefault():
    parser = Option_Parser.admixture_option_parser()
    options = parser.parse_args(['-s', '4',
                                 '-n', '0.003',
                                 '-d', '0.001',
                                 '--Neand1_sample_size', '3',
                                 '--Neand2_sample_size', '1',
                                 '-t', '450',
                                 '--time_N1_sample', '50',
                                 '--time_N2_sample', '120',
                                 '-l', '1e4',
                                 '-e', '1007',
                                 '-a', '1007',
                                 '-r', '3',
                                 '-g', 'AF_EU_1.5e-5',
                                 '-g', 'EU_AF_1.5e-5',
                                 '-g', 'AF_AS_0.78e-5',
                                 '-g', 'AS_AF_0.79e-5',
                                 '-g', 'EU_AS_3.11e-5',
                                 ])
    return Demography_Models.Base_demography(options)


def test_default_init(default):
    assert default.S_N1 == 2
    assert default.S_N2 == 2
    assert default.m_PULSE1 == 0.02
    assert default.m_PULSE2 == 0.0

    assert default.generation_time == 25
    assert default.length == 1e6
    assert default.recombination_rate == 1e-8
    assert default.mutation_rate == 1.2e-8
    assert default.N_A == 7310
    assert default.N_N1 == 1000
    assert default.N_N2 == 1000
    assert default.N_B == 1861
    assert default.N_AF0 == 14474
    assert default.N_EU0 == 1032
    assert default.N_AS0 == 554
    assert default.N_CH == 1000
    assert default.N_DE == 1000
    assert default.T_MH_CH == 7e6 / 25
    assert default.T_MH_N == 700e3 / 25
    assert default.T_DE_N == 500e3 / 25
    assert default.T_N1_N2 == ((145) * 1e3) / 25
    assert default.T_AF == 148e3 / 25
    assert default.T_B == 100e3 / 25
    assert default.T_PULSE1 == 55e3 / 25
    assert default.T_EU_AS == 23e3 / 25
    assert default.T_PULSE2 == 18e3 / 25
    assert default.T_ACL_GRW == 5.1e3 / 25
    assert default.r_EU1 == 0.00307
    assert default.r_AS1 == 0.0048
    assert default.N_EU == 512e3
    assert default.N_AS == 640e3
    assert default.N_AF == 424e3


def test_nondefault_init(nonDefault):
    assert nonDefault.S_N1 == 3
    assert nonDefault.S_N2 == 1
    assert nonDefault.m_PULSE1 == 0.003
    assert nonDefault.m_PULSE2 == 0.001

    assert nonDefault.length == 1e4
    assert nonDefault.T_N1_N2 == ((450) * 1e3) / 25


def test_population_map(default):
    pop_map = default.get_population_map()

    assert pop_map['N1'] == 0
    assert pop_map['N2'] == 1
    assert pop_map['AF'] == 2
    assert pop_map['EU'] == 3
    assert pop_map['AS'] == 4
    assert pop_map['CH'] == 5
    assert pop_map['DE'] == 6


def test_population_configuration(default):
    config = default.get_population_configuration()
    assert config[0].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            growth_rate=0.0).get_ll_representation()

    assert config[1].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            growth_rate=0.0).get_ll_representation()

    assert config[2].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=424e3,
            growth_rate=0.0166).get_ll_representation()

    assert config[3].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=512e3,
            growth_rate=0.0195).get_ll_representation()

    assert config[4].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=640e3,
            growth_rate=0.025).get_ll_representation()

    assert config[5].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            growth_rate=0.0).get_ll_representation()

    assert config[6].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            growth_rate=0.0).get_ll_representation()


def test_debug_population_configuration(default):
    config = default.get_debug_configuration()
    assert config[0].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2).get_ll_representation()

    assert config[1].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2).get_ll_representation()

    assert config[2].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=424e3,
            sample_size=2).get_ll_representation()

    assert config[3].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=512e3,
            sample_size=2).get_ll_representation()

    assert config[4].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=640e3,
            sample_size=2).get_ll_representation()

    assert config[5].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2).get_ll_representation()

    assert config[6].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2).get_ll_representation()


def test_samples(default, nonDefault):
    samples = default.get_samples()

    curr = 0
    gens = [55e3/25, 125e3/25, 0, 0, 0, 0, 50e3/25]
    numSamps = [2, 2, 2, 1006, 1008, 2, 2]
    for j, n in enumerate(numSamps):
        for i in range(curr, curr + n):
            assert samples[i] == msprime.Sample(population=j, time=gens[j])
        curr += n

    samples = nonDefault.get_samples()
    curr = 0
    gens = [50e3/25, 120e3/25, 0, 0, 0, 0, 50e3/25]
    numSamps = [3, 1, 3, 1007, 1007, 2, 2]
    for j, n in enumerate(numSamps):
        for i in range(curr, curr + n):
            assert samples[i] == msprime.Sample(population=j, time=gens[j])
        curr += n


def test_migration_matrix(default, nonDefault):
    '''
    NOTE ON ORDER
    the names are retrieved from the migration dictionary as
    COLUMN_ROW.  In the default, row 3, column 2 (2.5e-5) has
    a key of AF_EU.  Keys are constant in second underscore
    along a row
    '''
    mm = default.get_migration_matrix()
    assert mm == \
        [[0,   0,  0,  0,  0,  0,  0],
         [0,   0,  0,  0,  0,  0,  0],
         [0,   0,  0,  2.5e-5,  0.78e-5,  0,  0],  # AF
         [0,   0,  2.5e-5,  0,  3.11e-5,  0,  0],  # EU
         [0,   0,  0.78e-5,  3.11e-5,  0,  0,  0],  # AS
         [0,   0,  0,  0,  0,  0,  0],
         [0,   0,  0,  0,  0,  0,  0]]

    mm = nonDefault.get_migration_matrix()
    assert mm == \
        [[0,   0,  0,  0,  0,  0,  0],
         [0,   0,  0,  0,  0,  0,  0],
         [0,   0,  0,  1.5e-5,  0.79e-5,  0,  0],  # AF
         [0,   0,  1.5e-5,  0,  3.11e-5,  0,  0],  # EU
         [0,   0,  0.78e-5,  3.11e-5,  0,  0,  0],  # AS
         [0,   0,  0,  0,  0,  0,  0],
         [0,   0,  0,  0,  0,  0,  0]]


def test_construction():
    parser = Option_Parser.admixture_option_parser()
    options = parser.parse_args([])

    Demography_Models.Tenn_no_modern_migration(options).simulate(1)
    Demography_Models.Tenn_pulsed_migration(options).simulate(1)
    Demography_Models.Sriram_demography(options).simulate(1)
    Demography_Models.SplitPop_demography(options).simulate(1)
    Demography_Models.Out_of_africa_demography(options)
