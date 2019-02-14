import pytest
import collections
from Admixture_Simulation import merge_dict
import Admixture_Simulation
import Option_Parser
import Demography_Models
import File_Printer
import io


def test_merge_dict_build():
    d = collections.OrderedDict()
    merge_dict(d, ('10', 50, 100))
    assert d['10'] == [[50], [100]]
    merge_dict(d, ('10', 150, 200))
    assert d['10'] == [[50, 150], [100, 200]]
    merge_dict(d, ('10', 25, 30))
    assert d['10'] == [[25, 50, 150], [30, 100, 200]]

    d = collections.OrderedDict()
    merge_dict(d, ('10', 50, 100))
    merge_dict(d, ('10', 25, 30))
    merge_dict(d, ('10', 150, 200))
    assert d['10'] == [[25, 50, 150], [30, 100, 200]]

    d = collections.OrderedDict()
    merge_dict(d, ('10', 150, 200))
    merge_dict(d, ('10', 50, 100))
    merge_dict(d, ('10', 25, 30))
    assert d['10'] == [[25, 50, 150], [30, 100, 200]]


@pytest.fixture
def three_elem():
    d = collections.OrderedDict()
    merge_dict(d, ('10', 50, 100))
    merge_dict(d, ('10', 150, 200))
    merge_dict(d, ('10', 25, 30))
    return d


def test_merge_dict_edges(three_elem):
    d = three_elem
    merge_dict(d, ('10', 20, 27))
    assert d['10'] == [[20, 50, 150], [30, 100, 200]]
    merge_dict(d, ('10', 27, 35))
    assert d['10'] == [[20, 50, 150], [35, 100, 200]]
    merge_dict(d, ('10', 15, 40))
    assert d['10'] == [[15, 50, 150], [40, 100, 200]]

    merge_dict(d, ('10', 49, 99))
    assert d['10'] == [[15, 49, 150], [40, 100, 200]]
    merge_dict(d, ('10', 50, 101))
    assert d['10'] == [[15, 49, 150], [40, 101, 200]]
    merge_dict(d, ('10', 48, 102))
    assert d['10'] == [[15, 48, 150], [40, 102, 200]]

    merge_dict(d, ('10', 149, 199))
    assert d['10'] == [[15, 48, 149], [40, 102, 200]]
    merge_dict(d, ('10', 150, 201))
    assert d['10'] == [[15, 48, 149], [40, 102, 201]]
    merge_dict(d, ('10', 148, 202))
    assert d['10'] == [[15, 48, 148], [40, 102, 202]]


def test_merge_dict_overlap(three_elem):
    d = three_elem
    merge_dict(d, ('10', 20, 50))
    assert d['10'] == [[20, 150], [100, 200]]

    merge_dict(d, ('10', 400, 450))
    merge_dict(d, ('10', 300, 350))
    assert d['10'] == [[20, 150, 300, 400], [100, 200, 350, 450]]

    merge_dict(d, ('10', 15, 201))
    assert d['10'] == [[15, 300, 400], [201, 350, 450]]

    merge_dict(d, ('10', 15, 451))
    assert d['10'] == [[15], [451]]


def test_get_model():
    parser = Option_Parser.admixture_option_parser()
    opts = parser.parse_args([])

    assert type(Admixture_Simulation.get_model(opts)) == \
        Demography_Models.Tenn_demography
    opts.model = "SplitPop"
    assert type(Admixture_Simulation.get_model(opts)) == \
        Demography_Models.SplitPop_demography

    with pytest.raises(ValueError) as e:
        opts.model = "NONE"
        Admixture_Simulation.get_model(opts)

    assert "unsupported model: NONE" in str(e)


def test_tree_running():
    parser = Option_Parser.admixture_option_parser()
    opts = parser.parse_args([])
    sim = Demography_Models.Tenn_demography(opts).simulate(1)
    entries = Admixture_Simulation.get_haplo_entries(next(sim), opts)
    assert entries == {}


def test_get_samples():
    parser = Option_Parser.admixture_option_parser()
    opts = parser.parse_args([])

    samp = Admixture_Simulation.get_human_samples(opts)
    assert samp == range(6, 2020)

    opts.pop = "AFR"
    samp = Admixture_Simulation.get_human_samples(opts)
    assert samp == range(4, 6)

    with pytest.raises(ValueError) as e:
        opts.pop = "NONE"
        samp = Admixture_Simulation.get_human_samples(opts)

    assert "unknown human sample: NONE" in str(e)


def test_write_f4dstats():
    parser = Option_Parser.admixture_option_parser()
    opts = parser.parse_args([])

    with File_Printer.file_printer(opts) as fp:
        ind = io.StringIO()
        fp.writers['ind'] = ind
        eigen = io.BytesIO()
        fp.writers['eigen'] = eigen
        snp = io.StringIO()
        fp.writers['snp'] = snp

        model = Demography_Models.Tenn_demography(opts)
        sim = model.simulate(1)
        Admixture_Simulation.write_f4dstats(sim, fp, model)

        ind = ind.getvalue().split('\n')
        assert ind[0] == 'Sample_0\tU\tNeand1'
        assert ind[9] == 'Sample_9\tU\tEUR'
        snp = snp.getvalue().split('\n')
        snp0 = snp[0].split('\t')
        assert snp0[0] == 'rs1'
        assert snp0[1] == '1'
        assert snp0[3] == '129'
        assert snp0[4] == 'A'
        assert snp0[5] == 'T'

        snp9 = snp[9].split('\t')
        assert snp9[0] == 'rs10'
        assert snp9[1] == '1'
        assert snp9[3] == '965'
        assert snp9[4] == 'A'
        assert snp9[5] == 'T'
