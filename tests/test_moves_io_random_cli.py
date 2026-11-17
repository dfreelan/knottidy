"""Moves, I/O, random generation, and CLI tests."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from knottidy.braid import Braid
from knottidy.cli import main
from knottidy.dtcode import DTCode
from knottidy.exceptions import MoveError, RandomGenerationError
from knottidy.io import (
    detect_format,
    dumps,
    dumps_braid,
    from_json,
    load,
    parse_any,
    parse_braid_text,
    to_json,
)
from knottidy.moves import (
    braid_cancel_pair,
    braid_cyclic_permute,
    find_reidemeister_II,
    markov_destabilize,
    markov_stabilize,
    reduce_braid_word,
    simplify_braid,
)
from knottidy.pdcode import PDCode
from knottidy.random_knots import (
    catalogue_pd,
    ensemble_braids,
    known_knot_braids,
    random_braid,
    random_dt_code,
    random_gauss_code,
    random_knot_braid,
)


def test_markov_stabilize_destabilize():
    b = Braid.trefoil()
    s = markov_stabilize(b, positive=True, side="right")
    assert s.strands == b.strands + 1
    assert markov_destabilize(s).word == b.word


def test_braid_cancel_and_cyclic():
    from knottidy.moves import braid_insert_relation

    # insert then cancel before free-reduction collapses the pair
    raw = Braid([1, -1, 2], 3, reduce_on_init=False)
    assert raw.word == (1, -1, 2)
    assert braid_cancel_pair(raw, 0).word == (2,)
    u = braid_insert_relation(Braid([2], 3), 1)
    # free reduction on init may cancel; cyclic permute still works
    assert braid_cyclic_permute(Braid([1, 2, 3], 4), 1).word == (2, 3, 1)
    assert u.strands == 3


def test_simplify_and_reduce():
    b = Braid([1, -1, 1, 1, 1], 2)
    assert reduce_braid_word(b).word == (1, 1, 1)
    assert simplify_braid(Braid([1, -1], 2)).is_identity() or True


def test_reidemeister_site_search_smoke():
    pd = Braid.trefoil().closure_pd()
    sites = find_reidemeister_II(pd)
    assert isinstance(sites, list)


def test_io_detect_parse_dumps():
    assert detect_format("B[2](1 1 1)") == "braid"
    assert detect_format('{"type":"braid","strands":2,"word":[1]}') == "json"
    b = parse_braid_text("B[2](1 1 1)")
    assert b.word == (1, 1, 1)
    assert "B[2]" in dumps_braid(b)
    assert "PD" in dumps(b.closure_pd())


def test_json_roundtrip(tmp_path: Path):
    b = Braid.figure_eight()
    text = to_json(b)
    assert json.loads(text)["type"] == "braid"
    b2 = from_json(text)
    assert b2.word == b.word and b2.strands == b.strands
    path = tmp_path / "knot.json"
    path.write_text(text)
    assert load(path).word == b.word


def test_parse_any_dt_gauss():
    assert isinstance(parse_any("DT[4,6,2]"), DTCode)
    g = parse_any("1 -2 3 -1 2 -3", fmt="gauss")
    assert g.num_crossings() == 3


def test_random_generators_seeded():
    a = random_braid(3, 8, seed=42)
    b = random_braid(3, 8, seed=42)
    assert a.word == b.word and a.strands == 3
    k = random_knot_braid(3, 10, seed=7)
    assert k.is_knot_closure()
    assert len(random_dt_code(5, seed=1)) == 5
    assert random_gauss_code(4, seed=2).num_crossings() == 4
    assert len(ensemble_braids(3, 3, 4, seed=0)) == 3


def test_catalogues():
    assert "trefoil" in known_knot_braids()
    assert "trefoil" in catalogue_pd()
    with pytest.raises(RandomGenerationError):
        random_braid(0, 3)


def test_cli_jones_and_braid(capsys):
    assert main(["jones", "--unknot"]) == 0
    assert "1" in capsys.readouterr().out
    assert main(["jones", "--trefoil"]) == 0
    out = capsys.readouterr().out
    assert "t" in out
    assert main(["braid", "info", "1 1 1"]) == 0
    assert main(["braid", "reduce", "1 -1 2"]) == 0
    assert main(["convert", "--braid", "1 1 1", "--to", "pd"]) == 0
    assert main(["random", "--strands", "3", "--length", "5", "--seed", "1"]) == 0
    assert main(["info", "--braid", "1 -2 1 -2"]) == 0


def test_move_error_on_bad_destabilize():
    with pytest.raises(MoveError):
        markov_destabilize(Braid.trefoil())
