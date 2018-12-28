import Parameter_Sweeper
import pytest


def test_get_arguments_valid():
    out = Parameter_Sweeper.get_arguments(['-p', 'a;0'])[0]
    assert out['a'] == [0.0]

    out = Parameter_Sweeper.get_arguments(
        ['-p', 'a;0', '-p', 'a;1'])[0]
    assert out['a'] == [0.0, 1.0]

    out = Parameter_Sweeper.get_arguments(
        ['-p', 'a;1:3'])[0]
    assert out['a'] == [1, 2, 3]

    out = Parameter_Sweeper.get_arguments(
        ['-p', 'a;1,3,4'])[0]
    assert out['a'] == [1, 3, 4]

    out = Parameter_Sweeper.get_arguments(
        ['-p', 'a;1,3,4', '-p', 'b;0,4:6,8:9'])[0]
    assert out['a'] == [1, 3, 4]
    assert out['b'] == [0, 4, 5, 6, 8, 9]

    out = Parameter_Sweeper.get_arguments(
        ['-p', 'a;3.1e-5:0.1e-5:3.5e-5'])[0]
    assert out['a'] == [3.1e-5, 3.2e-5, 3.3e-5, 3.4e-5, 3.5e-5]


def test_get_arguments_invalid():
    with pytest.raises(ValueError) as e:
        Parameter_Sweeper.get_arguments([])
    assert 'No arguments provided' in str(e)

    with pytest.raises(ValueError) as e:
        Parameter_Sweeper.get_arguments(['-p' 'a;b'])
    assert 'Ill formed number in a;b' in str(e)

    with pytest.raises(SyntaxError) as e:
        Parameter_Sweeper.get_arguments(['-p' 'a;0:1:20:1'])

    with pytest.raises(ValueError) as e:
        Parameter_Sweeper.get_arguments(['-p' 'a;0:b'])
    assert 'Ill formed number in a;0:b' in str(e)


def test_get_permutations():
    args, fmts = Parameter_Sweeper.get_arguments(
        ['-p', 'a;0:1:2'])
    args = Parameter_Sweeper.get_permutations(args, fmts)
    assert args == ['-a 0.0', '-a 1.0', '-a 2.0']

    args, fmts = Parameter_Sweeper.get_arguments(
        ['-p', 'test;0:1:2'])
    args = Parameter_Sweeper.get_permutations(args, fmts)
    assert args == ['--test 0.0', '--test 1.0', '--test 2.0']

    args, fmts = Parameter_Sweeper.get_arguments(
        ['-p', 'a;0:1:2', '-p', 'b;1'])
    args = Parameter_Sweeper.get_permutations(args, fmts)
    assert args == ['-a 0.0 -b 1.0', '-a 1.0 -b 1.0', '-a 2.0 -b 1.0']

    args, fmts = Parameter_Sweeper.get_arguments(
        ['-p', 'a;0:1:2', '-p', 'b;1', '-p' 'c;2,3'])
    args = Parameter_Sweeper.get_permutations(args, fmts)
    assert args == [
        '-a 0.0 -b 1.0 -c 2.0',
        '-a 0.0 -b 1.0 -c 3.0',
        '-a 1.0 -b 1.0 -c 2.0',
        '-a 1.0 -b 1.0 -c 3.0',
        '-a 2.0 -b 1.0 -c 2.0',
        '-a 2.0 -b 1.0 -c 3.0']
