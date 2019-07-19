from Option_Parser import admixture_option_parser
from File_Printer import file_printer, get_basename
import Demography_Models
import pytest
import os
import sys
import io


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

    opts = admixture_option_parser().parse_args(['--options'])
    fp = file_printer(opts)
    assert fp.files['options'] == '*'
    assert fp.files['debug'] is None

    # flag set with filename on others
    opts = admixture_option_parser().parse_args(
        ['--debug', '--vcf', 'vcf.out'])
    fp = file_printer(opts)
    assert fp.files['options'] is None
    assert fp.files['debug'] == '*'
    assert fp.files['vcf'] == 'vcf.out'

    opts = admixture_option_parser().parse_args(
        ['--options', '--f4dstat', 'f4'])
    fp = file_printer(opts)
    assert fp.files['haplo'] is None
    assert fp.files['f4dstat'] == 'f4'
    assert fp.files['options'] == '*'

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
            ['--options', '--haplo'])
        file_printer(opts)
    assert 'Expected at most one output to stdout, got 2 instead.' in str(e)

    with pytest.raises(ValueError) as e:
        opts = admixture_option_parser().parse_args(
            ['--vcf', 'vc', '--haplo', '--options', '--debug'])
        file_printer(opts)
    assert 'Expected at most one output to stdout, got 3 instead.' in str(e)

    with pytest.raises(ValueError) as e:
        opts = admixture_option_parser().parse_args(
            ['--haplo', '--options', '--debug', '--pi'])
        file_printer(opts)
    assert 'Expected at most one output to stdout, got 4 instead.' in str(e)


def test_basename():
    assert get_basename('/test/one.txt') == 'one'
    assert get_basename('/test/two.txt.txt') == 'two.txt'
    assert get_basename('bare.txt') == 'bare'
    assert get_basename('bare') == 'bare'


def test_build_file_struct():
    fs = file_printer.file_struct("", "{}.txt")
    assert fs.non_default("test") == "test.txt"
    assert fs.non_default("outdir/test") == "outdir/test.txt"
    assert fs.non_default("test.txt") == "test.txt"
    assert fs.non_default("outdir/test.txt") == "outdir/test.txt"


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


def test_build_pi(tmp_path):
    opts = admixture_option_parser().parse_args([])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['pi'] is None

    opts = admixture_option_parser().parse_args(['--pi'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['pi'] == sys.stdout

    opts = admixture_option_parser().parse_args(
        ['--pi', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['pi'] == sys.stdout

    opts = admixture_option_parser().parse_args(
        ['--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['pi'] == os.path.join(
        str(tmp_path), 'pi.txt')

    opts = admixture_option_parser().parse_args(
        ['--pi', 'test.txt'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['pi'] == os.path.join(os.getcwd(), 'test.txt')

    opts = admixture_option_parser().parse_args(
        ['--pi', 'test.txt', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['pi'] == os.path.join(str(tmp_path), 'test.txt')


def test_build_option(tmp_path):
    opts = admixture_option_parser().parse_args([])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['options'] is None

    opts = admixture_option_parser().parse_args(['--options'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['options'] == sys.stdout

    opts = admixture_option_parser().parse_args(
        ['--options', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['options'] == sys.stdout

    opts = admixture_option_parser().parse_args(
        ['--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['options'] == os.path.join(
        str(tmp_path),
        'options.txt')

    opts = admixture_option_parser().parse_args(
        ['--options', 'test.txt'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['options'] == os.path.join(os.getcwd(), 'test.txt')

    opts = admixture_option_parser().parse_args(
        ['--options', 'test.txt', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['options'] == os.path.join(str(tmp_path), 'test.txt')


def test_build_vcf(tmp_path):
    opts = admixture_option_parser().parse_args([])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['vcf'] is None

    opts = admixture_option_parser().parse_args(
        ['--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['vcf'] == os.path.join(
        str(tmp_path), fp.get_filename('.vcf.gz'))

    opts = admixture_option_parser().parse_args(
        ['--vcf', 'test'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['vcf'] == \
        os.path.join(os.getcwd(), 'test.vcf.gz')

    opts = admixture_option_parser().parse_args(
        ['--vcf', 'test', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['vcf'] == \
        os.path.join(str(tmp_path), 'test.vcf.gz')


def test_build_popfile(tmp_path):
    opts = admixture_option_parser().parse_args([])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['popfile'] is None

    opts = admixture_option_parser().parse_args(
        ['--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['popfile'] == os.path.join(
        str(tmp_path), fp.get_filename('.popfile'))

    opts = admixture_option_parser().parse_args(
        ['--popfile', 'test'])
    fp = file_printer(opts)
    fp.build_files()
    assert fp.files['popfile'] == \
        os.path.join(os.getcwd(), 'test.popfile')

    opts = admixture_option_parser().parse_args(
        ['--popfile', 'test', '--out-dir', str(tmp_path)])
    fp = file_printer(opts)
    fp.build_files()
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

    singleopts = ['debug', 'options', 'haplo']
    for o in singleopts:
        opts = admixture_option_parser().parse_args(['--' + o])
        with file_printer(opts) as fp:
            for k, v in fp.writers.items():
                if v is not None:
                    assert k == o
            assert fp.writers[o] == sys.stdout
            assert fp.single_simulation_needed() is (o == 'haplo')
            assert fp.haplo_needed() is (o == 'haplo')
            assert fp.f4dstat_needed() is False

    opts = admixture_option_parser().parse_args(['--haplo', 'test.gz',
                                                 '--out-dir', str(tmp_path)])
    opts = admixture_option_parser().parse_args(['--haplo', 'test.vcf',
                                                 '--out-dir', str(tmp_path)])

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


def test_print_options():
    output = io.StringIO()
    opts = admixture_option_parser().parse_args([])
    with file_printer(opts) as fp:
        fp.writers['options'] = output
        fp.print_options()

        assert 'pop: nonAfr' in output.getvalue()
        assert 's_n1: 2' in output.getvalue()

    output = io.StringIO()
    opts = admixture_option_parser().parse_args(
        '-p AFR --Neand1_sample_size 3'.split())
    with file_printer(opts) as fp:
        fp.writers['options'] = output
        fp.print_options()

        assert 'pop: AFR' in output.getvalue()
        assert 's_n1: 3' in output.getvalue()


def test_print_debug():
    output = io.StringIO()
    opts = admixture_option_parser().parse_args([])
    model = Demography_Models.Base_demography(opts)
    with file_printer(opts) as fp:
        fp.writers['debug'] = output
        fp.print_debug(model)


def test_print_popfile():
    output = io.StringIO()
    opts = admixture_option_parser().parse_args([])
    with file_printer(opts) as fp:
        # when popfile writer is still none
        fp.print_popfile(None, None)
        model = Demography_Models.Base_demography(opts)
        ts = next(model.simulate(1))
        fp.writers['popfile'] = output
        fp.print_popfile(model, ts)

        output = output.getvalue().split('\n')
        assert "samp\tpop\tsuper_pop" == output[0]
        assert "msp_0\tNeand1\tNeand1" == output[1]
        assert "msp_1011\tDeni\tDeni" == output[-2]


def test_print_vcf():
    output = io.StringIO()
    opts = admixture_option_parser().parse_args([])
    with file_printer(opts) as fp:
        # when writer is still none
        fp.print_vcf(None)
        model = Demography_Models.Base_demography(opts)
        ts = next(model.simulate(1))
        fp.writers['vcf'] = output
        fp.print_vcf(ts)

        output = output.getvalue().split('\n')
        assert "##source=msprime 0.6.1" == output[1]
        assert len(output[5].split('\t')) == 1021


def test_print_pi():
    output = io.StringIO()
    opts = admixture_option_parser().parse_args([])
    with file_printer(opts) as fp:
        # when writer is still none
        fp.print_pi(None, None, None)
        model = Demography_Models.Base_demography(opts)
        ts = next(model.simulate(1))
        fp.writers['pi'] = output
        fp.print_pi(ts, model.get_sample_indices(),
                    model.get_population_map())

        lines = output.getvalue().split('\n')
        assert lines[0].split() == 'N1 N2 AF EU AS CH DE'.split()
        assert lines[1].split()[0] == '5.3e-05'
        assert lines[1].split()[1] == '7.5e-05'


def test_print_ils():
    output = io.StringIO()
    opts = admixture_option_parser().parse_args([])
    with file_printer(opts) as fp:
        # when writer is still none
        fp.print_ils(None)
        assert fp.ils_needed() is False

        haplos = {2: [[100, 300], [200, 400]]}
        fp.writers['ils'] = output
        assert fp.ils_needed() is True
        fp.print_ils(haplos)

        output = output.getvalue().split('\n')
        assert output[0] == '1\t100\t200\t2'
        assert output[1] == '1\t300\t400\t2'


def test_print_haplo():
    output = io.StringIO()
    opts = admixture_option_parser().parse_args([])
    with file_printer(opts) as fp:
        # when writer is still none
        fp.print_haplo(None)
        assert fp.haplo_needed() is False

        haplos = {2: [[100, 300], [200, 400]]}
        fp.writers['haplo'] = output
        assert fp.haplo_needed() is True
        fp.print_haplo(haplos)

        output = output.getvalue().split('\n')
        assert output[0] == '1\t100\t200\t2'
        assert output[1] == '1\t300\t400\t2'


def test_write_to():
    output = io.StringIO()
    opts = admixture_option_parser().parse_args([])
    with file_printer(opts) as fp:
        fp.writers['test'] = output
        fp.write_to('test', 'testing')
        assert output.getvalue() == 'testing'


def test_print_f4dstat():
    output = io.StringIO()
    opts = admixture_option_parser().parse_args([])
    with file_printer(opts) as fp:
        fp.print_f4dstat()

        fp.writers['f4dstat'] = output
        fp.files['eigen'] = 'eigen'
        fp.files['snp'] = 'snp'
        fp.files['ind'] = 'ind'
        fp.print_f4dstat()

        output = output.getvalue().split('\n')
        assert 'genotypename: eigen' == output[0]
        assert 'snpname: snp' == output[1]
        assert 'indivname: ind' == output[2]
        assert 'popfilename: sim.popfile_F4stat' == output[3]
