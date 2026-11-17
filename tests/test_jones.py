"""Jones polynomial known-value tests (unknot, trefoil, figure-eight)."""

from __future__ import annotations

from knottidy import Braid, PDCode, jones_polynomial, jones_trefoil, jones_unknot
from knottidy.jones import (
    LaurentPolynomial,
    is_unknot_jones,
    jones_figure_eight,
    jones_left_trefoil,
    kauffman_bracket,
    verify_trefoil_jones,
)


def test_unknot_jones():
    assert jones_polynomial(PDCode.unknot()) == jones_unknot()
    assert jones_polynomial(Braid.identity(1)) == jones_unknot()
    assert jones_polynomial(Braid([1, -1], 2)) == jones_unknot()
    assert is_unknot_jones(jones_unknot())


def test_trefoil_jones_right_and_left():
    rh = jones_polynomial(Braid.trefoil())
    assert rh == jones_trefoil()
    assert rh == LaurentPolynomial({1: 1, 3: 1, 4: -1}, variable="t")
    assert verify_trefoil_jones(rh)
    lh = jones_polynomial(Braid.trefoil().mirror())
    assert lh == jones_left_trefoil()
    assert verify_trefoil_jones(lh)
    assert jones_polynomial(PDCode.trefoil()) == jones_trefoil()


def test_figure_eight_jones():
    v = jones_polynomial(Braid.figure_eight())
    assert v == jones_figure_eight()
    assert v == jones_polynomial(PDCode.figure_eight())


def test_kauffman_bracket_unknot_and_trefoil():
    assert kauffman_bracket(PDCode.unknot()) == LaurentPolynomial.one("A")
    br = kauffman_bracket(Braid.trefoil().closure_pd())
    assert br.coefficient(0) == 0
    assert br  # nonempty


def test_laurent_arithmetic():
    p = LaurentPolynomial({1: 2, -1: 3}, variable="t")
    q = LaurentPolynomial({1: 1}, variable="t")
    assert (p + q).coefficient(1) == 3
    assert (p * q).coefficient(2) == 2
    assert (-p).coefficient(1) == -2
    assert p.shift(2).coefficient(3) == 2
    assert LaurentPolynomial.one("t") ** 5 == LaurentPolynomial.one("t")
    assert str(jones_trefoil())
    assert p.to_dict()["variable"] == "t"


def test_cinquefoil_differs_from_trefoil():
    c = jones_polynomial(Braid([1, 1, 1, 1, 1], 2))
    assert c != jones_trefoil()
    assert c != jones_unknot()
