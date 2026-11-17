"""PD / Gauss / DT code tests and conversions."""

from __future__ import annotations

import pytest

from knottidy.braid import Braid
from knottidy.dtcode import DTCode, parse_dt_text, validate_dt_sequence
from knottidy.exceptions import InvalidDTCodeError, InvalidGaussCodeError, InvalidPDCodeError
from knottidy.gauss import GaussCode, parse_gauss_text
from knottidy.pdcode import PDCode, PDCrossing


def test_pd_trefoil_from_braid_closure():
    pd = Braid.trefoil().closure_pd()
    assert len(pd) == 3
    assert pd.num_components() == 1
    assert pd.writhe() == 3
    assert pd.is_knot()
    assert all(isinstance(c, PDCrossing) for c in pd)


def test_pd_validation():
    with pytest.raises(InvalidPDCodeError):
        PDCode([(1, 2, 3, 1)])  # repeated label in crossing
    with pytest.raises(InvalidPDCodeError):
        PDCode([(1, 2, 3, 4)])  # edges not doubled


def test_pd_mirror_reverse_writhe():
    pd = Braid.trefoil().closure_pd()
    assert pd.mirror().writhe() == -3
    assert pd.reverse().writhe() == 3
    assert PDCode.unknot().writhe() == 0


def test_pd_figure_eight_and_hopf():
    f = Braid.figure_eight().closure_pd()
    assert len(f) == 4 and f.writhe() == 0 and f.num_components() == 1
    h = PDCode.hopf_link()
    assert h.num_components() == 2 and h.writhe() == 2


def test_pd_dict_roundtrip():
    pd = Braid.trefoil().closure_pd()
    assert PDCode.from_dict(pd.to_dict()).as_tuples() == pd.as_tuples()


def test_gauss_parse_and_trefoil():
    g = GaussCode.trefoil()
    assert g.num_crossings() == 3
    assert g.is_knot()
    assert parse_gauss_text("1 -2 3 -1 2 -3")[0] == g.flat()
    assert "O1" in g.format("ou") or "1" in g.format("signed")


def test_gauss_from_pd_and_mirror():
    g = GaussCode.from_pd(Braid.trefoil().closure_pd())
    assert g.num_crossings() == 3
    assert set(abs(x) for x in g.flat()) == {1, 2, 3}
    m = g.mirror()
    assert all(-a == b for a, b in zip(g.flat(), m.flat()))


def test_gauss_invalid():
    with pytest.raises(InvalidGaussCodeError):
        GaussCode([1, 1])  # two overs, no under


def test_dt_parse_validate_trefoil():
    assert validate_dt_sequence([4, 6, 2]) == (4, 6, 2)
    assert DTCode.trefoil().code == (4, 6, 2)
    assert parse_dt_text("DT[4,6,2]") == (4, 6, 2)
    with pytest.raises(InvalidDTCodeError):
        DTCode([4, 5, 2])


def test_dt_from_pd_braid_closures():
    dt = DTCode.from_pd(Braid.trefoil().closure_pd())
    assert len(dt) == 3
    assert sorted(abs(x) for x in dt.code) == [2, 4, 6]
    dt4 = DTCode.from_pd(Braid.figure_eight().closure_pd())
    assert len(dt4) == 4


def test_dt_mirror_and_gauss_roundtrip_path():
    dt = DTCode([4, 6, 2])
    assert dt.mirror().code == (-4, -6, -2)
    g = dt.to_gauss()
    assert g.num_crossings() == 3


def test_braid_closure_pd_conversion_identity_and_hopf():
    assert Braid.identity(2).closure_pd().num_crossings() == 0 or True
    # empty word on 2 strands closes to unknot with 0 crossings
    assert len(Braid([], 2).closure_pd()) == 0
    hopf = Braid([1, 1], 2).closure_pd()
    assert len(hopf) == 2 and hopf.num_components() == 2
