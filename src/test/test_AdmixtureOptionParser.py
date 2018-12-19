import pytest
import AdmixtureOptionParser
from argparse import ArgumentError


@pytest.fixture
def parser():
    return AdmixtureOptionParser.admixture_option_parser()


def compare_options(options, nondefault={}):
    defaults = {'pop': 'nonAfr',
                'model': 'Tenn',
                'seed': 1,
                'n1_admix_prop': 0.02,
                'n2_admix_prop': 0.0,
                's_n1': 2,
                's_n2': 2,
                't_n1_n2': 350,
                't_n1_sample': 55,
                't_n2_sample': 125,
                'haplo': 'F4Dstat',
                'length': 1e6,
                'EU_sample_size': 1006,
                'AS_sample_size': 1008,
                'AF_sample_size': 2,
                'initial_migrations': [
                    'AF_EU_2.5e-5', 'EU_AF_2.5e-5',
                    'AF_AS_0.78e-5', 'AS_AF_0.78e-5',
                    'EU_AS_3.11e-5', 'AS_EU_3.11e-5'],
                'later_migrations': [
                    'AF_B_15e-5', 'B_AF_15e-5'],
                'debug_file': None,
                'haplo_file': None,
                'vcf_file': None,
                'f4dstat_file': None,
                'option_file': None,
                'out_dir': None,
                }

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
    options = parser.parse_args(['-m', 'OUTPUT'])
    compare_options(options, {'model': 'OUTPUT'})
    options = parser.parse_args(['-s', '0'])
    compare_options(options, {'seed': 0})
    options = parser.parse_args(['-n', '0.003'])
    compare_options(options, {'n1_admix_prop': 0.003})
    options = parser.parse_args(['-d', '0.001'])
    compare_options(options, {'n2_admix_prop': 0.001})
    options = parser.parse_args(['-t', '450'])
    compare_options(options, {'t_n1_n2': 450})

    options = parser.parse_args(['--debug'])
    compare_options(options, {'debug_file': '*'})
    options = parser.parse_args(['--haplo'])
    compare_options(options, {'haplo_file': '*'})
    options = parser.parse_args(['--option'])
    compare_options(options, {'option_file': '*'})

    options = parser.parse_args(['--debug', 'test debug'])
    compare_options(options, {'debug_file': 'test debug'})
    options = parser.parse_args(['--haplo', 'test haplo'])
    compare_options(options, {'haplo_file': 'test haplo'})
    options = parser.parse_args(['--vcf', 'test vcf'])
    compare_options(options, {'vcf_file': 'test vcf'})
    options = parser.parse_args(['--f4dstat', 'test f4'])
    compare_options(options, {'f4dstat_file': 'test f4'})
    options = parser.parse_args(['--option', 'test option'])
    compare_options(options, {'option_file': 'test option'})
    options = parser.parse_args(['--out-dir', 'test dir'])
    compare_options(options, {'out_dir': 'test dir'})

    options = parser.parse_args(['-l', '1e7'])
    compare_options(options, {'length': 1e7})
    options = parser.parse_args(['-e', '1007'])
    compare_options(options, {'EU_sample_size': 1007})
    options = parser.parse_args(['-a', '1007'])
    compare_options(options, {'AS_sample_size': 1007})
    options = parser.parse_args(['-r', '3'])
    compare_options(options, {'AF_sample_size': 3})
    options = parser.parse_args(['-G', 'AF_B_13e-5'])
    compare_options(options, {'later_migrations': ['AF_B_13e-5']})
    options = parser.parse_args(['-g', 'EU_AS_3e-5'])
    compare_options(options, {'initial_migrations': ['EU_AS_3e-5']})


def test_type_error(parser):
    with pytest.raises(ArgumentError):
        parser.parse_args(['-s', 'string'])
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
        parser.parse_args(['--vcf'])
    with pytest.raises(ArgumentError):
        parser.parse_args(['--f4dstat'])
