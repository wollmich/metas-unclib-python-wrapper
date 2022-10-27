"""Test basic usage"""
from metas_unclib import ufloat, get_value, get_stdunc


def test_ufloat():
    a = ufloat(3.0, 0.3, desc="a")
    assert get_value(a) == 3.0
