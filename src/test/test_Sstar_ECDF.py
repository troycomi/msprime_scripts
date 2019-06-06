import pytest
import Sstar_ECDF
from io import StringIO
import os
import pandas as pd
from pandas.testing import assert_series_equal as ase
from pandas.testing import assert_frame_equal as afe
from numpy.testing import assert_array_equal as aae
from click.testing import CliRunner
import gzip


def test_main_build_null_db():
    runner = CliRunner()
    with runner.isolated_filesystem():
        # build chromosome file
        with open('chroms.txt', 'w') as output:
            output.write(
                '10\n'
                '12\n'
                '16\n'
            )

        # make fake window files
        with gzip.open('10.windowcalc.gz', 'wt') as output:
            output.write(
                'chrom\twinstart\twinend\tn_snps\tn_ind_snps\t'
                'n_region_ind_snps\tind_id\tpop\ts_star\n' +
                '2 0 50000 333 4 161 msp_110 EUR 25211\n'.replace(' ', '\t') +
                '2 0 50000 333 4 162 msp_111 EUR 25211\n'.replace(' ', '\t') +
                '2 0 50000 333 4 162 msp_112 EUR 25212\n'.replace(' ', '\t') +
                '2 0 50000 333 4 162 msp_113 EUR 25212\n'.replace(' ', '\t') +
                '2 0 50000 333 4 162 msp_114 EUR 25212\n'.replace(' ', '\t') +
                '2 0 50000 333 4 162 msp_115 EUR 25212\n'.replace(' ', '\t') +
                '2 0 50000 333 4 162 msp_116 EUR 25212\n'.replace(' ', '\t') +
                '2 0 50000 333 4 162 msp_116 ASN 25212\n'.replace(' ', '\t') +
                ''
            )

        with gzip.open('12.windowcalc.gz', 'wt') as output:
            output.write(
                'chrom\twinstart\twinend\tn_snps\tn_ind_snps\t'
                'n_region_ind_snps\tind_id\tpop\ts_star\n' +
                '3 0 50000 333 4 161 msp_110 ASN 25211\n'.replace(' ', '\t') +
                '3 0 50000 333 4 162 msp_111 EUR 25211\n'.replace(' ', '\t') +
                '3 0 50000 333 4 162 msp_112 ASN 25212\n'.replace(' ', '\t') +
                '3 0 50000 333 4 162 msp_113 EUR 25212\n'.replace(' ', '\t') +
                '3 0 50000 333 4 162 msp_114 ASN 25212\n'.replace(' ', '\t') +
                '3 0 50000 333 4 162 msp_115 EUR 25212\n'.replace(' ', '\t') +
                '3 0 50000 333 4 162 msp_116 EUR 25212\n'.replace(' ', '\t') +
                '3 0 50000 333 4 162 msp_116 ASN 25212\n'.replace(' ', '\t') +
                ''
            )

        with gzip.open('16.windowcalc.gz', 'wt') as output:
            output.write(
                'chrom\twinstart\twinend\tn_snps\tn_ind_snps\t'
                'n_region_ind_snps\tind_id\tpop\ts_star\n' +
                '1 0 50000 333 4 164 msp_110 EUR 25211\n'.replace(' ', '\t') +
                '1 0 50000 333 4 164 msp_111 ASN 25211\n'.replace(' ', '\t') +
                '1 0 50000 333 4 164 msp_112 EUR 25212\n'.replace(' ', '\t') +
                '1 0 50000 333 4 164 msp_113 ASN 25212\n'.replace(' ', '\t') +
                '1 0 50000 333 4 164 msp_114 ASN 25212\n'.replace(' ', '\t') +
                '1 0 50000 333 4 164 msp_115 EUR 25212\n'.replace(' ', '\t') +
                '1 0 50000 333 4 164 msp_116 EUR 25212\n'.replace(' ', '\t') +
                '1 0 50000 333 4 164 msp_116 ASN 25212\n'.replace(' ', '\t') +
                ''
            )

        with gzip.open('18.windowcalc.gz', 'wt') as output:
            output.write(
                'chrom\twinstart\twinend\tn_snps\tn_ind_snps\t'
                'n_region_ind_snps\tind_id\tpop\ts_star\n' +
                '1 0 50000 333 4 161 msp_110 EUR 25211\n'.replace(' ', '\t') +
                'NOT READ\n'
                '1 0 50000 333 4 162 msp_116 ASN 25212\n'.replace(' ', '\t') +
                ''
            )

        # build db
        result = runner.invoke(
            Sstar_ECDF.main,
            ' build-null-db '
            '--chr-list chroms.txt '
            '--outfile test.pkl'
        )

        assert result.exit_code == 0

        # check db
        assert os.path.exists('test.pkl')

        db = Sstar_ECDF.Null_DB()
        db.load('test.pkl')
        index = pd.MultiIndex.from_tuples([
            ('ASN', 161, 25211),
            ('ASN', 162, 25212),
            ('ASN', 164, 25211),
            ('ASN', 164, 25212),
            ('EUR', 161, 25211),
            ('EUR', 162, 25211),
            ('EUR', 162, 25212),
            ('EUR', 164, 25211),
            ('EUR', 164, 25212),
        ], names=['pop', 'n_region_ind_snps', 's_star'])
        ase(db.DB, pd.Series([1, 4, 1, 3, 1, 2, 8, 1, 3], index=index))


def test_main_generate_bed():
    runner = CliRunner()
    with runner.isolated_filesystem():
        # build up null db
        index = pd.MultiIndex.from_tuples([
            ('ASN', 161, 25211),
            ('ASN', 162, 25212),
            ('ASN', 164, 25211),
            ('ASN', 164, 25212),
            ('EUR', 161, 25211),
            ('EUR', 162, 25211),
            ('EUR', 162, 25212),
            ('EUR', 164, 25211),
            ('EUR', 164, 25212),
            ('EUR', 164, 25213),
            ('EUR', 164, 25214),
            ('EUR', 164, 25215),
        ], names=['pop', 'n_region_ind_snps', 's_star'])
        dat = list(range(12))
        null_db = Sstar_ECDF.Null_DB()
        null_db.DB = pd.Series(dat, index=index)
        null_db.save('null.pkl')

        # build up admixed window calc files
        with gzip.open('1.windowcalc.gz', 'wt') as window:
            window.write(
                'chrom\twinstart\twinend\tn_region_ind_snps\tind_id\tpop'
                '\ts_star\tnum_s_star_snps\tn_s_star_snps_hap1\tn_s_star_snps_hap2\n'
                '1\t0\t50000\t179\tmsp_110\tEUR\t85315\t10\t2\t8\n'
            )

        with gzip.open('2.windowcalc.gz', 'wt') as window:
            window.write(
                'chrom\twinstart\twinend\tn_region_ind_snps\tind_id\tpop'
                '\ts_star\tnum_s_star_snps\tn_s_star_snps_hap1\tn_s_star_snps_hap2\n'
                '2\t0\t50000\t179\tmsp_110\tEUR\t85315\t10\t2\t8\n'
            )

        # build tsv files
        with gzip.open('1.tsv.gz', 'wt') as tsv:
            tsv.write(
                'chr\tstart\tend\tisc\thaplotype\tpopulation'
                '\tmatch_pct\tpvalue\tmatching_windows\n'
                '1\t0\t50000\t199\tmsp_110:1\tEUR\t0.8\t0.03\t2217\n'
                '1\t0\t50000\t199\tmsp_110:2\tEUR\t0.8\t0.01\t2217\n'
            )

        with gzip.open('2.tsv.gz', 'wt') as tsv:
            tsv.write(
                'chr\tstart\tend\tisc\thaplotype\tpopulation'
                '\tmatch_pct\tpvalue\tmatching_windows\n'
                '2\t0\t50000\t199\tmsp_110:1\tEUR\t0.8\t0.01\t2217\n'
            )

        # chrom list
        with open('chroms.txt', 'w') as chrom:
            chrom.write('1\n2\n')

        # run main
        result = runner.invoke(
            Sstar_ECDF.main,
            'generate-bed '
            '--chr-list chroms.txt '
            '--outfile {chrom}.test.bed '
            '--null-db null.pkl '
        )
        assert result.exit_code == 0

        # check output files
        out1 = open('1.test.bed').readlines()
        assert out1 == [
            'msp_ID\tstart\tend\n',
            'msp_110:2_1\t0\t50000\n'
        ]
        out2 = open('2.test.bed').readlines()
        assert out2 == [
            'msp_ID\tstart\tend\n'
        ]


def test_get_chromosomes(mocker):
    data = StringIO(
        '1063938750\n'
        '1222356006\n'
        '1634154403\n'
    )
    mock_file = mocker.patch('Sstar_ECDF.open',
                             return_value=data)

    result = [e for e in Sstar_ECDF.get_chromosomes('mock')]

    mock_file.assert_called_once_with('mock', 'r')

    assert result == [
        '1063938750',
        '1222356006',
        '1634154403'
    ]

    mock_file.reset_mock()

    data = StringIO('')
    mock_file = mocker.patch('Sstar_ECDF.open',
                             return_value=data)
    assert [e for e in Sstar_ECDF.get_chromosomes('mock_empty')] == []

    mock_file.assert_called_once_with('mock_empty', 'r')


def test_region_files():
    result = [c for c in Sstar_ECDF.region_files('base', 'a b c'.split())]

    assert result == [
        os.path.join('base', 'a.windowcalc.gz'),
        os.path.join('base', 'b.windowcalc.gz'),
        os.path.join('base', 'c.windowcalc.gz'),
    ]

    result = [c for c in Sstar_ECDF.region_files('',
                                                 '1 2 3'.split(),
                                                 '{chrom}.txt')]

    assert result == [
        '1.txt',
        '2.txt',
        '3.txt',
    ]

    result = [c for c in Sstar_ECDF.region_files('base', [])]

    assert result == []


def test_process_windowcalc():
    infile = (
        'chrom\twinstart\twinend\tn_region_ind_snps\tind_id\tpop'
        '\ts_star\tnum_s_star_snps\tn_s_star_snps_hap1\tn_s_star_snps_hap2\n'
        '363104291\t0\t50000\t179\tmsp_110\tEUR\t85315\t10\t2\t8\n'
        '363104291\t0\t50000\t173\tmsp_111\tEUR\t47145\t4\t2\t2\n'
        '363104291\t0\t50000\t174\tmsp_112\tASN\t22028\t4\t2\t2\n'
        '363104291\t0\t50000\t180\tmsp_113\tEUR\t82882\t12\t7\t4\n'
        '363104291\t0\t50000\t184\tmsp_114\tEUR\t111004\t15\t8\t7\n'
        '363104291\t0\t50000\t174\tmsp_116\tEUR\t50976\t5\t4\t1\n'
        '363104291\t0\t50000\t176\tmsp_117\tASN\t70315\t7\t7\t0\n'
        '363104292\t0\t50000\t175\tmsp_118\tEUR\t67973\t3\t0\t3\n'
        '363104291\t0\t50000\t171\tmsp_115\tASN\t0\t0\t0\t0\n'
    )
    result = Sstar_ECDF.process_windowcalc(StringIO(infile))

    aae(result.columns, ['start', 'end', 'n_region_ind_snps',
                         'pop', 's_star', 'haplotype', 'msp_ID'])
    aae(result.start.values, [0, 0, 0, 0])
    aae(result.end.values, [50000, 50000, 50000, 50000])
    aae(result.n_region_ind_snps.values, [179, 174, 176, 175])
    aae(result['pop'].values, ['EUR', 'EUR', 'ASN', 'EUR'])
    aae(result.s_star.values, [85315, 50976, 70315, 67973])
    aae(result.haplotype.values, ['msp_110:2', 'msp_116:1',
                                  'msp_117:1', 'msp_118:2'])
    aae(result.msp_ID.values, ['msp_110:2_363104291', 'msp_116:1_363104291',
                               'msp_117:1_363104291', 'msp_118:2_363104292'])

    result = Sstar_ECDF.process_windowcalc(StringIO(infile), 1)

    aae(result.columns, ['start', 'end', 'n_region_ind_snps',
                         'pop', 's_star', 'haplotype', 'msp_ID'])
    aae(result.start.values, [0, 0])
    aae(result.end.values, [50000, 50000])
    aae(result.n_region_ind_snps.values, [176, 175])
    aae(result['pop'].values, ['ASN', 'EUR'])
    aae(result.s_star.values, [70315, 67973])
    aae(result.haplotype.values, ['msp_117:1', 'msp_118:2'])
    aae(result.msp_ID.values, ['msp_117:1_363104291', 'msp_118:2_363104292'])


def test_process_windowcalc_s_snps():
    infile = (
        'chrom\twinstart\twinend\tn_region_ind_snps\tind_id\tpop'
        '\ts_star\tnum_s_star_snps\tn_s_star_snps_hap1\tn_s_star_snps_hap2\n'
        '363104291\t0\t50000\t179\tmsp_110\tEUR\t85315\t8\t2\t8\n'
        '363104291\t0\t50000\t173\tmsp_111\tEUR\t47145\t2\t2\t2\n'
        '363104291\t0\t50000\t174\tmsp_112\tASN\t22028\t4\t2\t2\n'
        '363104291\t0\t50000\t180\tmsp_113\tEUR\t82882\t12\t7\t4\n'
        '363104291\t0\t50000\t184\tmsp_114\tEUR\t111004\t15\t8\t7\n'
        '363104291\t0\t50000\t174\tmsp_116\tEUR\t50976\t15\t4\t1\n'
        '363104291\t0\t50000\t176\tmsp_117\tASN\t70315\t7\t7\t0\n'
        '363104292\t0\t50000\t175\tmsp_118\tEUR\t67973\t3\t0\t3\n'
        '363104291\t0\t50000\t171\tmsp_115\tASN\t0\t0\t0\t0\n'
    )
    result = Sstar_ECDF.process_windowcalc(StringIO(infile))

    # 116 no longer passes as the num is too high
    # 111 should also not pass as it matches both (high overlap)
    aae(result.columns, ['start', 'end', 'n_region_ind_snps',
                         'pop', 's_star', 'haplotype', 'msp_ID'])
    aae(result.start.values, [0, 0, 0])
    aae(result.end.values, [50000, 50000, 50000])
    aae(result.n_region_ind_snps.values, [179, 176, 175])
    aae(result['pop'].values, ['EUR', 'ASN', 'EUR'])
    aae(result.s_star.values, [85315, 70315, 67973])
    aae(result.haplotype.values, ['msp_110:2',
                                  'msp_117:1', 'msp_118:2'])
    aae(result.msp_ID.values, ['msp_110:2_363104291',
                               'msp_117:1_363104291', 'msp_118:2_363104292'])


def test_filter_by_match():
    infile = (
        'chrom\twinstart\twinend\tn_region_ind_snps\tind_id\tpop'
        '\ts_star\tnum_s_star_snps\tn_s_star_snps_hap1\tn_s_star_snps_hap2\n'
        '363104291\t0\t50000\t179\tmsp_110\tEUR\t85315\t10\t2\t8\n'
        '363104291\t0\t50000\t173\tmsp_111\tEUR\t47145\t4\t2\t2\n'
        '363104291\t0\t50000\t174\tmsp_112\tASN\t22028\t4\t2\t2\n'
        '363104291\t0\t50000\t180\tmsp_113\tEUR\t82882\t12\t7\t4\n'
        '363104291\t0\t50000\t184\tmsp_114\tEUR\t111004\t15\t8\t7\n'
        '363104291\t0\t50000\t174\tmsp_116\tEUR\t50976\t5\t4\t1\n'
        '363104291\t0\t50000\t176\tmsp_117\tASN\t70315\t7\t7\t0\n'
        '363104292\t0\t50000\t175\tmsp_118\tEUR\t67973\t3\t0\t3\n'
        '363104291\t0\t50000\t171\tmsp_115\tASN\t0\t0\t0\t0\n'
    )
    window = Sstar_ECDF.process_windowcalc(StringIO(infile), 0.5)

    match_file = (
        'chr\tstart\tend\tisc\thaplotype\tpopulation'
        '\tmatch_pct\tpvalue\tmatching_windows\n'
        '363104291\t0\t50000\t199\tmsp_110:1\tEUR\t0.8\t0.01\t2217\n'
        '363104291\t0\t50000\t199\tmsp_110:2\tEUR\t0.8\t0.02\t2217\n'
    )

    result = Sstar_ECDF.filter_by_match(window, StringIO(match_file), 0.02)
    afe(result, pd.DataFrame({
        'start': [0],
        'end': [50000],
        'n_region_ind_snps': [179],
        'pop': ['EUR'],
        's_star': [85315],
        'msp_ID': ['msp_110:2_363104291'],
    }))

    # result empty as pvalue too low
    result = Sstar_ECDF.filter_by_match(window, StringIO(match_file), 0.01)
    afe(result, pd.DataFrame({
        'start': [],
        'end': [],
        'n_region_ind_snps': [],
        'pop': [],
        's_star': [],
        'msp_ID': [],
    }), check_index_type=False, check_dtype=False)


def test_filter_by_sstar():
    df = pd.DataFrame({
        'start': [0, 0, 0],
        'end': [50000, 50000, 50000],
        'n_region_ind_snps': [179, 179, 178],
        'pop': ['EUR', 'EUR', 'EUR'],
        's_star': [85315, 85314, 85313],
        'msp_ID': ['msp_110:2_363104292', 'msp_110:2_363104293',
                   'msp_110:2_363104291'],
    })

    null_db = Sstar_ECDF.Null_DB()
    index = pd.MultiIndex.from_tuples([
        ('EUR', 179, 25211),
        ('EUR', 179, 85313),
        ('EUR', 179, 85314),
        ('EUR', 179, 85315),
        ('EUR', 179, 85316),
        ('EUR', 178, 25211),
    ], names=['pop', 'n_region_ind_snps', 's_star'])
    null_db.DB = pd.Series([2, 2, 2, 2, 2, 2], index=index)

    # nothing
    result = Sstar_ECDF.filter_by_sstar(df, null_db, -1)
    afe(result, pd.DataFrame({
        'start': [],
        'end': [],
        'msp_ID': []
    }), check_like=True, check_dtype=False)

    # last
    result = Sstar_ECDF.filter_by_sstar(df, null_db, 0)
    afe(result, pd.DataFrame({
        'start': [0],
        'end': [50000],
        'msp_ID': ['msp_110:2_363104291']
    }), check_like=True, check_dtype=False)

    result = Sstar_ECDF.filter_by_sstar(df, null_db, 0.2)
    afe(result, pd.DataFrame({
        'start': [0, 0],
        'end': [50000, 50000],
        'msp_ID': ['msp_110:2_363104291', 'msp_110:2_363104292']
    }), check_like=True, check_dtype=False)

    # Note, groupby sorts by pop, then n_region
    result = Sstar_ECDF.filter_by_sstar(df, null_db, 0.4)
    afe(result, pd.DataFrame({
        'start': [0, 0, 0],
        'end': [50000, 50000, 50000],
        'msp_ID': ['msp_110:2_363104291', 'msp_110:2_363104292',
                   'msp_110:2_363104293']
    }), check_like=True, check_dtype=False)

    # empty window input
    df = pd.DataFrame({
        'start': [],
        'end': [],
        'n_region_ind_snps': [],
        'pop': [],
        's_star': [],
        'msp_ID': [],
    })
    result = Sstar_ECDF.filter_by_sstar(df, null_db, 0.05)
    afe(result, pd.DataFrame({
        'start': [],
        'end': [],
        'msp_ID': []
    }), check_like=True, check_dtype=False)


@pytest.fixture
def null_db():
    return Sstar_ECDF.Null_DB()


def test_null_db_init(null_db):
    ase(null_db.DB, pd.Series(dtype='int64'))


def test_null_db_read_windowcalc(null_db):
    infile = StringIO(
        'chrom\twinstart\twinend\tn_snps\tn_ind_snps\t'
        'n_region_ind_snps\tind_id\tpop\ts_star\n' +
        '1 0 50000 333 4 161 msp_110 EUR 25211\n'.replace(' ', '\t') +
        '1 0 50000 333 4 162 msp_111 EUR 25211\n'.replace(' ', '\t') +
        '1 0 50000 333 4 162 msp_112 EUR 25212\n'.replace(' ', '\t') +
        '1 0 50000 333 4 162 msp_113 EUR 25212\n'.replace(' ', '\t') +
        '1 0 50000 333 4 162 msp_114 EUR 25212\n'.replace(' ', '\t') +
        '1 0 50000 333 4 162 msp_115 EUR 25212\n'.replace(' ', '\t') +
        '1 0 50000 333 4 162 msp_116 EUR 25212\n'.replace(' ', '\t') +
        '1 0 50000 333 4 162 msp_116 ASN 25212\n'.replace(' ', '\t') +
        '1 0 50000 333 4 162 msp_116 EUR 0\n'.replace(' ', '\t') +  # skip
        ''
    )
    null_db.read_windowcalc(infile)

    index = pd.MultiIndex.from_tuples([
        ('ASN', 162, 25212),
        ('EUR', 161, 25211),
        ('EUR', 162, 25211),
        ('EUR', 162, 25212),
    ], names=['pop', 'n_region_ind_snps', 's_star'])
    ase(null_db.DB, pd.Series([1, 1, 1, 5], index=index))

    infile = StringIO(
        'chrom\twinstart\twinend\tn_snps\tn_ind_snps\t'
        'n_region_ind_snps\tind_id\tpop\ts_star\n' +
        '1 0 50000 333 4 161 msp_110 EUR 25211\n'.replace(' ', '\t') +
        '1 0 50000 333 4 162 msp_111 EUR 25211\n'.replace(' ', '\t') +
        ''
    )
    null_db.read_windowcalc(infile)

    index = pd.MultiIndex.from_tuples([
        ('ASN', 162, 25212),
        ('EUR', 161, 25211),
        ('EUR', 162, 25211),
        ('EUR', 162, 25212),
    ], names=['pop', 'n_region_ind_snps', 's_star'])
    ase(null_db.DB, pd.Series([1, 2, 2, 5], index=index))


def test_null_db_ecdf(null_db):
    index = pd.MultiIndex.from_tuples([
        ('ASN', 162, 25212),
        ('EUR', 161, 25211),
        ('EUR', 162, 25211),
        ('EUR', 162, 25212),
    ], names=['pop', 'n_region_ind_snps', 's_star'])
    null_db.DB = pd.Series([1, 2, 2, 6], index=index)

    ecdf = null_db.ecdf('EUR', 162, [25210, 25211, 25212, 25213])
    aae(ecdf, [0, 0.25, 1, 1])
    ecdf = null_db.ecdf('EUR', 161, [25210, 25211, 25212, 25213])
    aae(ecdf, [0, 1, 1, 1])
    ecdf = null_db.ecdf('ASN', 162, [25210, 25211, 25212, 25213])
    aae(ecdf, [0, 0, 1, 1])


def test_null_db_ecdf_missing_n(null_db):
    index = pd.MultiIndex.from_tuples([
        ('EUR', 161, 25208),
        ('EUR', 161, 25211),
        ('EUR', 165, 25210),
        ('EUR', 165, 25211),
        ('EUR', 165, 25212),
    ], names=['pop', 'n_region_ind_snps', 's_star'])
    null_db.DB = pd.Series([1, 1, 2, 2, 6], index=index)

    ecdf = null_db.ecdf('EUR', 50, [25210, 25211, 25212, 25213])
    aae(ecdf, [0.5, 1, 1, 1])

    ecdf = null_db.ecdf('EUR', 250, [25210, 25211, 25212, 25213])
    aae(ecdf, [0.2, 0.4, 1, 1])

    # nearer 161
    ecdf = null_db.ecdf('EUR', 162, [25209, 25210, 25211, 25212, 25213])
    aae(ecdf, [3/16, 5/16, 10/16, 1, 1])

    ecdf = null_db.ecdf('EUR', 163, [25209, 25210, 25211, 25212, 25213])
    aae(ecdf, [1/12, 3/12, 6/12, 1, 1])

    ecdf = null_db.ecdf('EUR', 164, [25209, 25210, 25211, 25212, 25213])
    aae(ecdf, [1/32, 7/32, 14/32, 1, 1])


def test_null_db_get_sstar(null_db):
    index = pd.MultiIndex.from_tuples([
        ('EUR', 161, 25208),
        ('EUR', 161, 25211),
        ('EUR', 165, 25210),
        ('EUR', 165, 25211),
        ('EUR', 165, 25212),
    ], names=['pop', 'n_region_ind_snps', 's_star'])
    null_db.DB = pd.Series([1, 1, 2, 2, 6], index=index)

    # not a valid pop
    with pytest.raises(ValueError) as e:
        null_db.get_sstar('ASN', 161)
    assert 'Population "ASN" not found in null database' in str(e)

    # exists, get values
    ase(null_db.get_sstar('EUR', 161),
        pd.Series([1, 1],
                  index=pd.Index([25208, 25211],
                                 name='s_star')))

    # low, get min
    ase(null_db.get_sstar('EUR', 157),
        pd.Series([1, 1],
                  index=pd.Index([25208, 25211],
                                 name='s_star')))

    # high, get max
    ase(null_db.get_sstar('EUR', 167),
        pd.Series([2, 2, 6],
                  index=pd.Index([25210, 25211, 25212],
                                 name='s_star')))

    # interpolate, nearer 161
    ase(null_db.get_sstar('EUR', 162),
        pd.Series([3/4, 2/4, 5/4, 6/4],
                  index=pd.Index([25208, 25210, 25211, 25212],
                                 name='s_star')))

    # interpolate, halfway
    ase(null_db.get_sstar('EUR', 163),
        pd.Series([1/2, 1, 3/2, 3],
                  index=pd.Index([25208, 25210, 25211, 25212],
                                 name='s_star')))

    # interpolate, nearer 165
    ase(null_db.get_sstar('EUR', 164),
        pd.Series([1/4, 6/4, 7/4, 18/4],
                  index=pd.Index([25208, 25210, 25211, 25212],
                                 name='s_star')))


def test_null_db_save_load(null_db, tmp_path):
    index = pd.MultiIndex.from_tuples([
        ('ASN', 162, 25212),
        ('EUR', 162, 25212),
        ('EUR', 162, 25211),
        ('EUR', 161, 25211),
    ], names=['pop', 'n_region_ind_snps', 's_star'])
    null_db.DB = pd.Series([1, 2, 2, 6], index=index)

    ase(null_db.DB, pd.Series([1, 2, 2, 6], index=index))

    pkl = tmp_path / 'db.pkl'
    null_db.save(pkl)

    assert os.path.exists(pkl)
    null_db.DB = None

    null_db.load(pkl)
    # now sorted
    index = pd.MultiIndex.from_tuples([
        ('ASN', 162, 25212),
        ('EUR', 161, 25211),
        ('EUR', 162, 25211),
        ('EUR', 162, 25212),
    ], names=['pop', 'n_region_ind_snps', 's_star'])
    ase(null_db.DB, pd.Series([1, 6, 2, 2], index=index))
