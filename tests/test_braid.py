"""Tests for braid words, reduction, Markov moves, and closures."""

from __future__ import annotations

import pytest

from knottidy.braid import Braid, format_braid_word, free_reduce, parse_braid_word
from knottidy.exceptions import InvalidBraidError


def test_free_reduce_cancels_inverses():
    assert free_reduce([1, -1, 2]) == (2,)
    assert free_reduce([1, 1, -1, -1]) == ()
    assert free_reduce([1, 2, -2, -1]) == ()
    assert free_reduce([2, 1, -2]) == (2, 1, -2)


def test_braid_init_reduces():
    b = Braid([1, -1, 1], strands=2)
    assert b.word == (1,)
    assert Braid.identity(3).is_identity()


def test_parse_and_format():
    word, n = parse_braid_word("1 -2 1")
    assert word == (1, -2, 1) and n == 3
    b = Braid.from_string("s1 s2^-1 s1")
    assert b.word == (1, -2, 1)
    assert "1" in format_braid_word(b.word, "signed")
    assert "s1" in format_braid_word(b.word, "sigma")
    assert format_braid_word((), "signed") == "e"


def test_algebra_inverse_multiply_pow():
    b = Braid([1, 2], 3)
    assert (b * b.inverse()).is_identity()
    assert (b**3).word == (1, 2, 1, 2, 1, 2)
    assert (b**0).is_identity()
    assert (~b).word == (-2, -1)


def test_writhe_permutation_components():
    t = Braid.trefoil()
    assert t.writhe() == 3
    assert t.is_knot_closure()
    assert t.num_components() == 1
    assert Braid([1, 1], 2).num_components() == 2
    assert Braid([1, -1], 2).is_pure()


def test_mirror_reverse_stabilize():
    t = Braid.trefoil()
    assert t.mirror().word == (-1, -1, -1)
    assert t.reverse().word == (1, 1, 1)
    s = t.stabilize(positive=True, side="right")
    assert s.strands == 3 and s.word[-1] == 2
    assert s.can_destabilize()
    assert s.destabilize().word == t.word


def test_conjugate_and_add_strands():
    b = Braid([1], 2)
    c = b.conjugate(Braid([2], 3))
    assert c.strands == 3
    assert b.add_strands(2, "right").strands == 4


def test_invalid_generators():
    with pytest.raises(InvalidBraidError):
        Braid([0], 2)
    with pytest.raises(InvalidBraidError):
        Braid([3], 2)
    with pytest.raises(InvalidBraidError):
        parse_braid_word("not_a_braid")


def test_closure_pd_crossing_count():
    assert len(Braid.trefoil().closure_pd()) == 3
    assert len(Braid.figure_eight().closure_pd()) == 4
    assert len(Braid.identity(1).closure_pd()) == 0
    assert Braid([1, -1], 2).closure_pd().num_components() == 1


def test_named_constructors_and_dict():
    assert Braid.generator(1, 3, power=2).word == (1, 1)
    d = Braid.trefoil().to_dict()
    assert Braid.from_dict(d).word == (1, 1, 1)
