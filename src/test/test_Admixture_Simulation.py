import pytest
import collections
from Admixture_Simulation import merge_dict


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
