"""Random generators for braids, PD sketches, and combinatorial knot data.

All randomness goes through :class:`random.Random` so callers can pass an
explicit ``seed`` for reproducibility. Generation is combinatorial (random
braid words, random DT-like sequences) rather than geometrically uniform in
the sense of a probability measure on knot space.
"""

from __future__ import annotations

import random
from collections.abc import Sequence
from typing import Optional

from knottidy.braid import Braid
from knottidy.dtcode import DTCode
from knottidy.exceptions import RandomGenerationError
from knottidy.gauss import GaussCode
from knottidy.pdcode import PDCode

def _rng(seed: Optional[int] = None) -> random.Random:
    return random.Random(seed)

def random_braid(
    strands: int,
    length: int,
    *,
    seed: Optional[int] = None,
    allow_inverses: bool = True,
    reduce: bool = True,
) -> Braid:
    """Sample a random braid word in ``B_strands`` of the given free length."""
    if strands < 1:
        raise RandomGenerationError("strands must be >= 1")
    if length < 0:
        raise RandomGenerationError("length must be >= 0")
    if strands == 1:
        return Braid.identity(1)
    rng = _rng(seed)
    max_gen = strands - 1
    word: list[int] = []
    for _ in range(length):
        idx = rng.randint(1, max_gen)
        if allow_inverses and rng.random() < 0.5:
            idx = -idx
        word.append(idx)
    return Braid(word, strands, reduce_on_init=reduce)

def random_positive_braid(
    strands: int,
    length: int,
    *,
    seed: Optional[int] = None,
) -> Braid:
    """Random positive braid (no inverse generators)."""
    return random_braid(strands, length, seed=seed, allow_inverses=False, reduce=False)

def random_knot_braid(
    strands: int,
    length: int,
    *,
    seed: Optional[int] = None,
    max_tries: int = 200,
) -> Braid:
    """Sample braids whose closure has a single component."""
    if strands < 2:
        raise RandomGenerationError("need strands >= 2 for a knot braid")
    rng = _rng(seed)
    for _ in range(max_tries):
        b = random_braid(strands, length, seed=rng.randint(0, 2**30))
        if b.is_knot_closure():
            return b
    raise RandomGenerationError(
        f"failed to sample a knot-closure braid in {max_tries} tries"
    )

def random_dt_code(
    crossings: int,
    *,
    seed: Optional[int] = None,
    signed: bool = True,
) -> DTCode:
    """Random DT sequence: random permutation of ``{2,4,...,2n}`` with random signs."""
    if crossings < 0:
        raise RandomGenerationError("crossings must be >= 0")
    if crossings == 0:
        return DTCode.empty()
    rng = _rng(seed)
    evens = list(range(2, 2 * crossings + 1, 2))
    rng.shuffle(evens)
    if signed:
        code = [e if rng.random() < 0.5 else -e for e in evens]
    else:
        code = evens
    return DTCode(code)

def random_gauss_code(
    crossings: int,
    *,
    seed: Optional[int] = None,
) -> GaussCode:
    """Random interlaced Gauss code for a knot with ``crossings`` crossings.

    Builds a random perfect pairing placement on a circle and assigns random
    over/under consistently.
    """
    if crossings < 0:
        raise RandomGenerationError("crossings must be >= 0")
    if crossings == 0:
        return GaussCode.empty()
    rng = _rng(seed)
    slots = []
    for k in range(1, crossings + 1):
        slots.extend([k, k])
    rng.shuffle(slots)
    seen: set[int] = set()
    seq: list[int] = []
    first_is_over: dict[int, bool] = {
        k: rng.random() < 0.5 for k in range(1, crossings + 1)
    }
    for k in slots:
        if k not in seen:
            seen.add(k)
            seq.append(k if first_is_over[k] else -k)
        else:
            seq.append(-k if first_is_over[k] else k)
    return GaussCode(seq)

def random_braid_closure_pd(
    strands: int,
    length: int,
    *,
    seed: Optional[int] = None,
) -> PDCode:
    """Random braid → PD via closure."""
    b = random_braid(strands, length, seed=seed)
    return b.closure_pd()

def ensemble_braids(
    count: int,
    strands: int,
    length: int,
    *,
    seed: Optional[int] = None,
    knot_closures_only: bool = False,
) -> list[Braid]:
    """Generate an ensemble of random braids."""
    if count < 0:
        raise RandomGenerationError("count must be >= 0")
    rng = _rng(seed)
    out: list[Braid] = []
    tries = 0
    max_tries = max(1000, count * 50)
    while len(out) < count and tries < max_tries:
        tries += 1
        b = random_braid(strands, length, seed=rng.randint(0, 2**30))
        if knot_closures_only and not b.is_knot_closure():
            continue
        out.append(b)
    if len(out) < count:
        raise RandomGenerationError(
            f"only sampled {len(out)}/{count} braids matching constraints"
        )
    return out

def known_knot_braids() -> dict[str, Braid]:
    """Dictionary of named example braids."""
    return {
        "unknot": Braid.identity(1),
        "unknot_2": Braid([1, -1], 2),
        "trefoil": Braid.trefoil(),
        "left_trefoil": Braid.trefoil().mirror(),
        "figure_eight": Braid.figure_eight(),
        "hopf_link": Braid([1, 1], 2),
        "cinquefoil": Braid([1, 1, 1, 1, 1], 2),
        "3_1_3strand": Braid([1, 1, 1, 2, -1, 2], 3),
        "solomon_seal": Braid([1, 1, 1, 1, 1], 2),  # 5_1
        "stevedore_cand": Braid([1, 1, 1, -2, -1, 2, -1, 2], 3),
    }

def catalogue_pd() -> dict[str, PDCode]:
    """Dictionary of named PD diagrams."""
    return {
        "unknot": PDCode.unknot(),
        "trefoil": PDCode.trefoil(),
        "figure_eight": PDCode.figure_eight(),
        "hopf_link": PDCode.hopf_link(),
        "trefoil_braid_closure": Braid.trefoil().closure_pd(),
        "figure_eight_braid_closure": Braid.figure_eight().closure_pd(),
    }

__all__ = [
    "random_braid",
    "random_positive_braid",
    "random_knot_braid",
    "random_dt_code",
    "random_gauss_code",
    "random_braid_closure_pd",
    "ensemble_braids",
    "known_knot_braids",
    "catalogue_pd",
]
