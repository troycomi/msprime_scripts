from Option_Parser import admixture_option_parser as parser
import Demography_Models
import msprime


def test_Tenn_debug():
    config = Demography_Models.Tenn_demography(parser().parse_args([]))\
        .get_debug_configuration()

    assert config[0].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2,
            growth_rate=0.0).get_ll_representation()

    assert config[1].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2,
            growth_rate=0.0).get_ll_representation()

    assert config[2].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=424e3,
            sample_size=2,
            growth_rate=0.0166).get_ll_representation()

    assert config[3].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=512e3,
            sample_size=2,
            growth_rate=0.0195).get_ll_representation()

    assert config[4].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=640e3,
            sample_size=2,
            growth_rate=0.025).get_ll_representation()

    assert config[5].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2,
            growth_rate=0.0).get_ll_representation()

    assert config[6].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2,
            growth_rate=0.0).get_ll_representation()


def test_SplitPop_debug():
    config = Demography_Models.SplitPop_demography(parser().parse_args([]))\
        .get_debug_configuration()

    assert config[0].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2,
            growth_rate=0.0).get_ll_representation()

    assert config[1].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2,
            growth_rate=0.0).get_ll_representation()

    assert config[2].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=424e3,
            sample_size=2,
            growth_rate=0.0166).get_ll_representation()

    assert config[3].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=512e3,
            sample_size=2,
            growth_rate=0.0195).get_ll_representation()

    assert config[4].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=640e3,
            sample_size=2,
            growth_rate=0.025).get_ll_representation()

    assert config[5].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2,
            growth_rate=0.0).get_ll_representation()

    assert config[6].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2,
            growth_rate=0.0).get_ll_representation()


def test_OOA_debug():
    config = Demography_Models.\
        Out_of_africa_demography(parser().parse_args([]))\
        .get_debug_configuration()

    assert config[0].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2,
            growth_rate=0.0).get_ll_representation()

    assert config[1].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=12300,
            sample_size=2,
            growth_rate=0.0).get_ll_representation()

    ll = config[2].get_ll_representation()
    ll['initial_size'] //= 1
    assert ll == \
        msprime.PopulationConfiguration(
            initial_size=29725,
            sample_size=2,
            growth_rate=0.004).get_ll_representation()

    ll = config[3].get_ll_representation()
    ll['initial_size'] //= 1
    assert ll == \
        msprime.PopulationConfiguration(
            initial_size=54090,
            sample_size=2,
            growth_rate=0.0055).get_ll_representation()

    assert config[4].get_ll_representation() == \
        msprime.PopulationConfiguration(
            initial_size=1000,
            sample_size=2,
            growth_rate=0.0).get_ll_representation()
