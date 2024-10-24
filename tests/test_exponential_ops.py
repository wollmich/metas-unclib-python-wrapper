"""
Use this test file to run unit tests against the code.
To use, install pytest via pip and run `pytest` in the root folder of the repo.
"""

from metas_unclib import *
import numpy as np


def test_power():
    a = ufloat(-2)
    result25 = (
        1.7319121124709869e-15 + 5.656854249492381j
    )  # should be equal to value of a**2.5

    assert isinstance(a**2.5, ucomplex)
    assert np.isclose(get_value(a**2.5), result25)

    assert isinstance(a**2, ufloat)
    assert get_value(a**2) == 4

    assert isinstance(a ** (2 + 1j), ucomplex)
    assert np.isclose(
        get_value(a ** (2 + 1j)), (0.13296730803542664 + 0.11044808147333202j)
    )

    assert isinstance((a ** np.array([2, 4]))[0], ufloat)
    assert get_value((a ** np.array([2, 4]))[0]) == 4
    assert isinstance((a ** np.array([2.5, 4.0]))[0], ucomplex)
    assert np.isclose(get_value((a ** np.array([2.5, 4.0]))[0]), result25)


def test_sqrt():
    assert isinstance(ufloat(4).sqrt(), ufloat)
    assert get_value(ufloat(4).sqrt()) == 2.0

    assert isinstance(ufloat(-1).sqrt(), ucomplex)
    assert np.isclose(get_value(ufloat(-1).sqrt()), 1.0j)


def test_log():
    assert isinstance(ufloat(42).log(), ufloat)
    assert np.isclose(get_value(ufloat(42).log()), 3.7376696182833684)

    assert isinstance(ufloat(-1).log(), ucomplex)
    assert np.isclose(get_value(ufloat(-1).log()), (1j*np.pi))


def test_log10():
    assert isinstance(ufloat(100).log10(), ufloat)
    assert get_value(ufloat(100).log10()) == 2.0

    assert isinstance(ufloat(-1).log10(), ucomplex)
    assert np.isclose(get_value(ufloat(-1).log10()), 1.3643763538418412j)
