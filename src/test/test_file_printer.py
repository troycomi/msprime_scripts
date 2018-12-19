from AdmixtureOptionParser import admixture_option_parser
from File_Printer import file_printer
import pytest
import os
import sys


def test_defaults():
    opts = admixture_option_parser().parse_args([])
    fp = file_printer(opts)
    for f in fp.files.values():
        assert f is None
    assert fp.out_dir is None


def test_options():
    # single flags set
    opts = admixture_option_parser().parse_args(['--haplo'])
    fp = file_printer(opts)
    assert fp.files['haplo'] == '*'
    assert fp.files['debug'] is None

    opts = admixture_option_parser().parse_args(['--option'])
    fp = file_printer(opts)
    assert fp.files['option'] == '*'
    assert fp.files['debug'] is None

    # flag set with filename on others
    opts = admixture_option_parser().parse_args(
        ['--debug', '--vcf', 'vcf.out'])
    fp = file_printer(opts)
    assert fp.files['option'] is None
    assert fp.files['debug'] == '*'
    assert fp.files['vcf'] == 'vcf.out'

    opts = admixture_option_parser().parse_args(
        ['--option', '--f4dstat', 'f4'])
    fp = file_printer(opts)
    assert fp.files['haplo'] is None
    assert fp.files['f4dstat'] == 'f4'
    assert fp.files['option'] == '*'

    # outdir set ok
    opts = admixture_option_parser().parse_args(
        ['--out-dir', 'outd'])
    fp = file_printer(opts)
    assert fp.files['haplo'] is None
    assert fp.out_dir == 'outd'

    # note this is valid but silly, no output will go to outdir
    opts = admixture_option_parser().parse_args(
        ['--out-dir', 'outd', '--haplo'])
    fp = file_printer(opts)
    assert fp.out_dir == 'outd'
    assert fp.files['haplo'] == '*'

    opts = admixture_option_parser().parse_args(
        ['--out-dir', 'outd', '--haplo', 'hap'])
    fp = file_printer(opts)
    assert fp.out_dir == 'outd'
    assert fp.files['haplo'] == 'hap'


def test_exeception_options():
    # too many flags set
    with pytest.raises(ValueError) as e:
        opts = admixture_option_parser().parse_args(
            ['--option', '--haplo'])
        file_printer(opts)
    assert 'Expected at most one output to stdout, got 2 instead.' in str(e)

    with pytest.raises(ValueError) as e:
        opts = admixture_option_parser().parse_args(
            ['--vcf', 'vc', '--haplo', '--option', '--debug'])
        file_printer(opts)
    assert 'Expected at most one output to stdout, got 3 instead.' in str(e)

    with pytest.raises(ValueError) as e:
        opts = admixture_option_parser().parse_args(
            ['--haplo', '--option', '--debug'])
        file_printer(opts)
    assert 'Expected at most one output to stdout, got 3 instead.' in str(e)


def test_build_out_dir(tmp_path):
    # default set to wd
    opts = admixture_option_parser().parse_args(
        ['--haplo', 'hap'])
    fp = file_printer(opts)
    fp.build_out_dir()
    assert fp.out_dir == os.getcwd()

    # new dir (pytest starts with existing dir)
    os.rmdir(tmp_path)
    assert not os.path.exists(tmp_path)
    opts = admixture_option_parser().parse_args(
        ['--out-dir', str(tmp_path), '--haplo', 'hap'])
    fp = file_printer(opts)
    fp.build_out_dir()
    assert os.path.exists(tmp_path)

    # redoing with existing dir should not fail
    opts = admixture_option_parser().parse_args(
        ['--out-dir', str(tmp_path), '--haplo', 'hap'])
    fp = file_printer(opts)
    fp.build_out_dir()


def test_build_debug(tmp_path):
    opts = admixture_option_parser().parse_args([])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['debug'] == sys.stdout

    opts = admixture_option_parser().parse_args(['--debug'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['debug'] == sys.stdout

    opts = admixture_option_parser().parse_args(
        ['--debug', '--vcf', 'vcfbase'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['debug'] == sys.stdout

    opts = admixture_option_parser().parse_args(
        ['--debug', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['debug'] == sys.stdout

    opts = admixture_option_parser().parse_args(
        ['--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['debug'] == os.path.join(str(tmp_path), 'debug.txt')

    opts = admixture_option_parser().parse_args(
        ['--haplo'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['debug'] is None

    opts = admixture_option_parser().parse_args(
        ['--debug', 'test.txt'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['debug'] == os.path.join(os.getcwd(), 'test.txt')

    opts = admixture_option_parser().parse_args(
        ['--debug', 'test.txt', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()


def test_build_haplo(tmp_path):
    opts = admixture_option_parser().parse_args([])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['haplo'] is None

    opts = admixture_option_parser().parse_args(['--haplo'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['haplo'] == sys.stdout

    opts = admixture_option_parser().parse_args(
        ['--haplo', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['haplo'] == sys.stdout

    opts = admixture_option_parser().parse_args(
        ['--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['haplo'] == os.path.join(
        str(tmp_path),
        fp.get_filename('.bed.merged.gz'))

    opts = admixture_option_parser().parse_args(
        ['--haplo', 'test.txt'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['haplo'] == os.path.join(os.getcwd(), 'test.txt')

    opts = admixture_option_parser().parse_args(
        ['--haplo', 'test.txt', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['haplo'] == os.path.join(str(tmp_path), 'test.txt')


def test_build_option(tmp_path):
    opts = admixture_option_parser().parse_args([])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['option'] is None

    opts = admixture_option_parser().parse_args(['--option'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['option'] == sys.stdout

    opts = admixture_option_parser().parse_args(
        ['--option', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['option'] == sys.stdout

    opts = admixture_option_parser().parse_args(
        ['--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['option'] == os.path.join(
        str(tmp_path),
        'options.txt')

    opts = admixture_option_parser().parse_args(
        ['--option', 'test.txt'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['option'] == os.path.join(os.getcwd(), 'test.txt')

    opts = admixture_option_parser().parse_args(
        ['--option', 'test.txt', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['option'] == os.path.join(str(tmp_path), 'test.txt')


def test_build_vcf(tmp_path):
    opts = admixture_option_parser().parse_args([])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['vcf'] is None
    assert fp.files['popfile'] is None

    opts = admixture_option_parser().parse_args(
        ['--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['vcf'] == os.path.join(
        str(tmp_path), fp.get_filename('.vcf.gz'))
    assert fp.files['popfile'] == os.path.join(
        str(tmp_path), fp.get_filename('.popfile'))

    opts = admixture_option_parser().parse_args(
        ['--vcf', 'test'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['vcf'] == \
        os.path.join(os.getcwd(), 'test.vcf.gz')
    assert fp.files['popfile'] == \
        os.path.join(os.getcwd(), 'test.popfile')

    opts = admixture_option_parser().parse_args(
        ['--vcf', 'test', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['vcf'] == \
        os.path.join(str(tmp_path), 'test.vcf.gz')
    assert fp.files['popfile'] == \
        os.path.join(str(tmp_path), 'test.popfile')


def test_build_f4dstat(tmp_path):
    opts = admixture_option_parser().parse_args([])
    fp = file_printer(opts)
    fp.build_files()
    files = ['f4dstat', 'eigen', 'snp', 'ind']
    for f in files:
        assert fp.files[f] is None

    opts = admixture_option_parser().parse_args(
        ['--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    outbases = ['parfile.F4stat', 'eigenstratgeno', 'snp', 'ind']
    for f, ob in zip(files, outbases):
        assert fp.files[f] == os.path.join(
            str(tmp_path), fp.get_f4dstat_filename(ob))

    opts = admixture_option_parser().parse_args(
        ['--f4dstat', 'test'])
    fp = file_printer(opts)
    fp.build_files()
    outbases = ['parfile.F4stat.{}.gz', 'eigenstratgeno.{}.gz',
                'snp.{}.gz', 'ind.{}.gz']
    for f, ob in zip(files, outbases):
        assert fp.files[f] == \
            os.path.join(os.getcwd(), ob.format('test'))

    opts = admixture_option_parser().parse_args(
        ['--f4dstat', 'test', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    for f, ob in zip(files, outbases):
        assert fp.files[f] == \
            os.path.join(str(tmp_path), ob.format('test'))


def test_open_writers(tmp_path):
    opts = admixture_option_parser().parse_args([])
    with file_printer(opts) as fp:
        for k, v in fp.writers.items():
            if v is not None:
                assert k == 'debug'
        assert fp.writers['debug'] == sys.stdout
        assert fp.haplo_needed() is False
        assert fp.f4dstat_needed() is False

    singleopts = ['debug', 'option', 'haplo']
    for o in singleopts:
        opts = admixture_option_parser().parse_args(['--' + o])
        with file_printer(opts) as fp:
            for k, v in fp.writers.items():
                if v is not None:
                    assert k == o
            assert fp.writers[o] == sys.stdout
            assert fp.haplo_needed() is (o == 'haplo')
            assert fp.f4dstat_needed() is False

    opts = admixture_option_parser().parse_args(['--vcf', 'test',
                                                 '--out-dir', str(tmp_path)])
    with file_printer(opts) as fp:
        for k, v in fp.writers.items():
            if v is not None:
                assert k == 'vcf' or k == 'popfile'

    opts = admixture_option_parser().parse_args(['--f4dstat', 'test',
                                                 '--out-dir', str(tmp_path)])
    with file_printer(opts) as fp:
        for k, v in fp.writers.items():
            if v is not None:
                assert k == 'f4dstat' \
                    or k == 'eigen' \
                    or k == 'snp' \
                    or k == 'ind'
        assert fp.haplo_needed() is False
        assert fp.f4dstat_needed() is True
