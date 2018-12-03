import pytest
import AdmixtureOptionParser
from argparse import ArgumentError


@pytest.fixture
def parser():
    return AdmixtureOptionParser.admixture_option_parser()


def compare_options(options, nondefault={}):
    defaults = {'pop': 'nonAfr',
                'outdir': 'Tenn',
                'seed': 1,
                'pulses': 2,
                'n1_admix_prop': 0.02,
                'n2_admix_prop': 0.0,
                't_n1_n2': 350,
                'haplo': 'F4Dstat',
                'length': 1e6,
                'EU_sample_size': 1006,
                'AS_sample_size': 1008,
                'AF_sample_size': 2,
                'm_AF_B': 15e-5,
                'm_AF_AS': 0.78e-5,
                'm_AF_EU': 2.5e-5,
                'm_EU_AS': 3.11e-5}

    for k, v in nondefault.items():
        defaults[k] = v

    for k in options.__dict__.keys():
        assert options.__dict__[k] == defaults[k]


def test_defaults(parser):
    options = parser.parse_args([])
    compare_options(options)


def test_single_changes(parser):
    options = parser.parse_args(['-p', 'EUR'])
    compare_options(options, {'pop': 'EUR'})
    options = parser.parse_args(['-o', 'OUTPUT'])
    compare_options(options, {'outdir': 'OUTPUT'})
    options = parser.parse_args(['-s', '0'])
    compare_options(options, {'seed': 0})
    options = parser.parse_args(['-i', '3'])
    compare_options(options, {'pulses': 3})
    options = parser.parse_args(['-n', '0.003'])
    compare_options(options, {'n1_admix_prop': 0.003})
    options = parser.parse_args(['-d', '0.001'])
    compare_options(options, {'n2_admix_prop': 0.001})
    options = parser.parse_args(['-t', '450'])
    compare_options(options, {'t_n1_n2': 450})
    options = parser.parse_args(['-c', 'haplo'])
    compare_options(options, {'haplo': 'haplo'})
    options = parser.parse_args(['-l', '1e7'])
    compare_options(options, {'length': 1e7})
    options = parser.parse_args(['-e', '1007'])
    compare_options(options, {'EU_sample_size': 1007})
    options = parser.parse_args(['-a', '1007'])
    compare_options(options, {'AS_sample_size': 1007})
    options = parser.parse_args(['-r', '3'])
    compare_options(options, {'AF_sample_size': 3})
    options = parser.parse_args(['--migration_AF_B', '13e-5'])
    compare_options(options, {'m_AF_B': 13e-5})
    options = parser.parse_args(['--migration_AF_AS', '0.7e-5'])
    compare_options(options, {'m_AF_AS': 0.7e-5})
    options = parser.parse_args(['--migration_AF_EU', '2e-5'])
    compare_options(options, {'m_AF_EU': 2e-5})
    options = parser.parse_args(['--migration_EU_AS', '3e-5'])
    compare_options(options, {'m_EU_AS': 3e-5})


def test_type_error(parser):
    with pytest.raises(ArgumentError):
        parser.parse_args(['-s', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['-i', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['-n', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['-d', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['-t', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['-l', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['-e', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['-a', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['-r', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['--migration_AF_B', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['--migration_AF_AS', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['--migration_AF_EU', 'string'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['--migration_EU_AS', 'string'])
