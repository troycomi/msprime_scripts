import Demography_Models
import msprime


def test_population():
    pop = Demography_Models.Population('test', 1e3, 0.0, 2, 75)
    assert pop.abbreviation == 'test'
    assert pop.size == 1e3
    assert pop.rate == 0.0
    assert pop.samples == 2
    assert pop.generations == 75


def test_population_sample():
    pop = Demography_Models.Population('test', 1e3, 0.0, 2, 75)
    sample = pop.get_sample(1)
    assert len(sample) == 2
    for i in range(2):
        assert sample[i] == msprime.Sample(population=1, time=75)

    pop = Demography_Models.Population('test', 4e3, 0.1, 3, 85)
    sample = pop.get_sample(2)
    assert len(sample) == 3
    for i in range(3):
        assert sample[i] == msprime.Sample(population=2, time=85)


def test_population_configuration():
    pop = Demography_Models.Population('test', 1e3, 0.0, 2, 75)
    config = pop.get_configuration()
    assert config.get_ll_representation() == msprime.PopulationConfiguration(
        initial_size=1e3, growth_rate=0.0).get_ll_representation()

    pop = Demography_Models.Population('test', 4e3, 0.1, 3, 85)
    config = pop.get_configuration()
    assert config.get_ll_representation() == msprime.PopulationConfiguration(
        initial_size=4e3, growth_rate=0.1).get_ll_representation()


def test_population_debug_configuration():
    pop = Demography_Models.Population('test', 1e3, 0.0, 2, 75)
    config = pop.get_debug_configuration()
    assert config.get_ll_representation() == msprime.PopulationConfiguration(
        initial_size=1e3, sample_size=2).get_ll_representation()
    config = pop.get_debug_configuration(includeRate=True)
    assert config.get_ll_representation() == msprime.PopulationConfiguration(
        initial_size=1e3, sample_size=2,
        growth_rate=0.0).get_ll_representation()

    pop = Demography_Models.Population('test', 4e3, 0.1, 3, 85)
    config = pop.get_debug_configuration()
    assert config.get_ll_representation() == msprime.PopulationConfiguration(
        initial_size=4e3, sample_size=2).get_ll_representation()
    config = pop.get_debug_configuration(includeRate=True)
    assert config.get_ll_representation() == msprime.PopulationConfiguration(
        initial_size=4e3, sample_size=2,
        growth_rate=0.1).get_ll_representation()
