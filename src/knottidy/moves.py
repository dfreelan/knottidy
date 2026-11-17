"""Reidemeister-adjacent and Markov moves on diagrams and braids."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Optional

from knottidy.braid import Braid
from knottidy.exceptions import InvalidBraidError, InvalidPDCodeError, MoveError
from knottidy.pdcode import PDCode, PDCrossing

def markov_conjugate(braid: Braid, conjugator: Braid) -> Braid:
    """Markov conjugation: ``β ↦ γ β γ^{-1}``."""
    if not isinstance(braid, Braid) or not isinstance(conjugator, Braid):
        raise MoveError("markov_conjugate expects Braid arguments")
    return braid.conjugate(conjugator)

def markov_stabilize(braid: Braid, positive: bool = True, side: str = "right") -> Braid:
    """Markov stabilization (±)."""
    if not isinstance(braid, Braid):
        raise MoveError("markov_stabilize expects a Braid")
    return braid.stabilize(positive=positive, side=side)

def markov_destabilize(braid: Braid) -> Braid:
    """Markov destabilization when applicable."""
    if not isinstance(braid, Braid):
        raise MoveError("markov_destabilize expects a Braid")
    try:
        return braid.destabilize()
    except InvalidBraidError as exc:
        raise MoveError(str(exc)) from exc

def markov_moves_applicable(braid: Braid) -> dict[str, bool]:
    """Report which elementary Markov moves are immediately applicable."""
    if not isinstance(braid, Braid):
        raise MoveError("expected a Braid")
    return {
        "stabilize_pos_right": True,
        "stabilize_neg_right": True,
        "stabilize_pos_left": True,
        "stabilize_neg_left": True,
        "destabilize": braid.can_destabilize(),
        "conjugate": True,
    }

def braid_insert_relation(braid: Braid, position: int) -> Braid:
    """Insert a free-cancellation pair ``σ_i σ_i^{-1}`` at ``position`` (Reidemeister-II-like on braids).

    ``position`` is an index in ``0..len(word)`` where the pair is spliced in.
    The generator index ``i`` defaults to 1 when strands allow, else the only
    available generator.
    """
    if not isinstance(braid, Braid):
        raise MoveError("expected a Braid")
    if braid.strands < 2:
        raise MoveError("need at least 2 strands to insert a generator pair")
    if position < 0 or position > len(braid.word):
        raise MoveError(f"position {position} out of range")
    i = 1
    word = list(braid.word)
    word[position:position] = [i, -i]
    return Braid(word, braid.strands)

def braid_cancel_pair(braid: Braid, position: int) -> Braid:
    """Cancel a free pair at ``position`` if ``word[position:position+2]`` are inverses."""
    if not isinstance(braid, Braid):
        raise MoveError("expected a Braid")
    w = braid.word
    if position < 0 or position + 1 >= len(w):
        raise MoveError("position out of range for cancel")
    if w[position] != -w[position + 1]:
        raise MoveError("no inverse pair at the given position")
    new_word = w[:position] + w[position + 2 :]
    return Braid(new_word, braid.strands)

def braid_cyclic_permute(braid: Braid, shift: int = 1) -> Braid:
    """Cyclic conjugation by rotating the word (special conjugacy)."""
    if not isinstance(braid, Braid):
        raise MoveError("expected a Braid")
    w = braid.word
    if not w:
        return braid
    shift %= len(w)
    return Braid(w[shift:] + w[:shift], braid.strands)

def reduce_braid_word(braid: Braid) -> Braid:
    """Free-reduce and apply heuristic Artin shorteners."""
    if not isinstance(braid, Braid):
        raise MoveError("expected a Braid")
    return braid.free_reduce().commute_form()

@dataclass(frozen=True)
class ReidemeisterSite:
    """Description of a location where a Reidemeister move applies."""

    move: str  # "I", "II", "III"
    crossings: tuple[int, ...]
    detail: str = ""

def find_reidemeister_I(pd: PDCode) -> list[ReidemeisterSite]:
    """Find combinatorial Reidemeister I curl sites on a PD code."""
    if not isinstance(pd, PDCode):
        raise MoveError("expected PDCode")
    sites: list[ReidemeisterSite] = []
    if len(pd) == 0:
        return sites
    succ = pd.successor_map()
    for i, cr in enumerate(pd.crossings):
        if succ[cr.c] == cr.b:
            sites.append(ReidemeisterSite("I", (i,), "under-to-over curl"))
        elif succ[cr.d] == cr.a:
            sites.append(ReidemeisterSite("I", (i,), "over-to-under curl"))
        elif succ[cr.c] == cr.a and succ.get(cr.d) == cr.b:
            sites.append(ReidemeisterSite("I", (i,), "double-loop"))
        elif succ[cr.c] == cr.b and succ[cr.d] == cr.a:
            sites.append(ReidemeisterSite("I", (i,), "curl-swap"))
    return sites

def find_reidemeister_II(pd: PDCode) -> list[ReidemeisterSite]:
    """Find Reidemeister II digons between two crossings.

    Two crossings form an RII site when two parallel strands connect them
    with opposite crossing signs (or cancelling smoothings), i.e. they share
    two connecting arcs and the crossings can cancel.
    """
    if not isinstance(pd, PDCode):
        raise MoveError("expected PDCode")
    sites: list[ReidemeisterSite] = []
    n = len(pd)
    if n < 2:
        return sites
    succ = pd.successor_map()
    inc: dict[int, tuple[int, str]] = {}
    for i, cr in enumerate(pd.crossings):
        inc[cr.a] = (i, "U")
        inc[cr.b] = (i, "O")

    for i in range(n):
        for j in range(i + 1, n):
            ci, cj = pd.crossings[i], pd.crossings[j]
            outs_i = {ci.c, ci.d}
            outs_j = {cj.c, cj.d}
            ins_i = {ci.a, ci.b}
            ins_j = {cj.a, cj.b}
            i_to_j = outs_i & ins_j
            j_to_i = outs_j & ins_i
            if len(i_to_j) == 2 and len(j_to_i) == 0:
                pass
            if len(i_to_j) >= 1 and len(j_to_i) >= 1:
                if len(i_to_j) + len(j_to_i) >= 2:
                    sites.append(
                        ReidemeisterSite(
                            "II",
                            (i, j),
                            f"connections i→j={sorted(i_to_j)} j→i={sorted(j_to_i)}",
                        )
                    )
    uniq: dict[frozenset[int], ReidemeisterSite] = {}
    for s in sites:
        uniq[frozenset(s.crossings)] = s
    return list(uniq.values())

def find_reidemeister_III(pd: PDCode) -> list[ReidemeisterSite]:
    """Find triangles of three crossings that may admit Reidemeister III.

    Heuristic: three crossings each pair-connected by one arc forming a
    triangle.
    """
    if not isinstance(pd, PDCode):
        raise MoveError("expected PDCode")
    sites: list[ReidemeisterSite] = []
    n = len(pd)
    if n < 3:
        return sites
    adj: dict[int, set[int]] = {i: set() for i in range(n)}
    edge_owner_out: dict[int, int] = {}
    edge_owner_in: dict[int, int] = {}
    for i, cr in enumerate(pd.crossings):
        edge_owner_out[cr.c] = i
        edge_owner_out[cr.d] = i
        edge_owner_in[cr.a] = i
        edge_owner_in[cr.b] = i
    for e, oi in edge_owner_out.items():
        if e in edge_owner_in:
            oj = edge_owner_in[e]
            if oi != oj:
                adj[oi].add(oj)
                adj[oj].add(oi)
    for i in range(n):
        for j in adj[i]:
            if j <= i:
                continue
            for k in adj[i] & adj[j]:
                if k <= j:
                    continue
                sites.append(ReidemeisterSite("III", (i, j, k), "triangle"))
    return sites

def apply_reidemeister_I_remove(pd: PDCode, crossing_index: int) -> PDCode:
    """Remove a Reidemeister I curl when a detected site exists."""
    if not isinstance(pd, PDCode):
        raise MoveError("expected PDCode")
    if crossing_index < 0 or crossing_index >= len(pd):
        raise MoveError("crossing index out of range")
    sites = [s for s in find_reidemeister_I(pd) if crossing_index in s.crossings]
    if not sites:
        raise MoveError(f"no Reidemeister I site at crossing {crossing_index}")
    if len(pd) == 1:
        return PDCode.empty()
    raise MoveError("RI removal for multi-crossing diagrams is not automated yet")

def apply_reidemeister_II_remove(pd: PDCode, i: int, j: int) -> PDCode:
    """Remove a cancelling Reidemeister II pair when the diagram is only that pair."""
    if not isinstance(pd, PDCode):
        raise MoveError("expected PDCode")
    if i == j or min(i, j) < 0 or max(i, j) >= len(pd):
        raise MoveError("invalid RII indices")
    sites = [s for s in find_reidemeister_II(pd) if set(s.crossings) == {i, j}]
    if not sites:
        raise MoveError(f"no RII site at crossings {i}, {j}")
    if len(pd) == 2:
        return PDCode.empty()
    raise MoveError("RII removal for larger diagrams is not automated yet")

def simplify_pd(pd: PDCode, max_passes: int = 32) -> PDCode:
    """Greedily apply RI / RII removals to reduce crossing number."""
    if not isinstance(pd, PDCode):
        raise MoveError("expected PDCode")
    current = pd
    for _ in range(max_passes):
        progressed = False
        for site in find_reidemeister_I(current):
            try:
                current = apply_reidemeister_I_remove(current, site.crossings[0])
                progressed = True
                break
            except MoveError:
                continue
        if progressed:
            continue
        for site in find_reidemeister_II(current):
            if len(site.crossings) != 2:
                continue
            a, b = site.crossings
            try:
                current = apply_reidemeister_II_remove(current, a, b)
                progressed = True
                break
            except MoveError:
                continue
        if not progressed:
            break
    return current

def all_reidemeister_sites(pd: PDCode) -> list[ReidemeisterSite]:
    """Collect RI, RII, and RIII sites."""
    return find_reidemeister_I(pd) + find_reidemeister_II(pd) + find_reidemeister_III(pd)

def simplify_braid(braid: Braid, max_passes: int = 16) -> Braid:
    """Try free reduction, Artin shorten, and destabilization."""
    current = reduce_braid_word(braid)
    for _ in range(max_passes):
        if current.can_destabilize():
            try:
                current = current.destabilize()
                current = reduce_braid_word(current)
                continue
            except InvalidBraidError:
                pass
        broke = False
        w = current.word
        for i in range(len(w) - 1):
            if w[i] == -w[i + 1]:
                current = braid_cancel_pair(current, i)
                broke = True
                break
        if not broke:
            break
    return current

__all__ = [
    "ReidemeisterSite",
    "markov_conjugate",
    "markov_stabilize",
    "markov_destabilize",
    "markov_moves_applicable",
    "braid_insert_relation",
    "braid_cancel_pair",
    "braid_cyclic_permute",
    "reduce_braid_word",
    "find_reidemeister_I",
    "find_reidemeister_II",
    "find_reidemeister_III",
    "apply_reidemeister_I_remove",
    "apply_reidemeister_II_remove",
    "simplify_pd",
    "all_reidemeister_sites",
    "simplify_braid",
]
