# knottidy

**knottidy** is a pure-Python library for combinatorial knot theory: Artin braid
words, planar diagram (PD) codes, Dowker–Thistlethwaite (DT) codes, Gauss codes,
and the Jones polynomial via the Kauffman bracket.

It targets Python 3.10+ and depends only on the standard library.

## Features

- **Braid groups** — words in Artin generators, free reduction, mirrors,
  reverses, strand flips, Markov stabilization / destabilization, conjugacy
  helpers, strand permutation and component count for closures
- **PD codes** — validation, components, writhe, mirror / reverse, conversion
  from braid closures
- **DT & Gauss codes** — parsing, validation, conversion to/from PD
- **Jones polynomial** — exact integer Laurent polynomials; state-sum Kauffman
  bracket with writhe normalisation (practical for small crossing number)
- **Moves** — Markov moves on braids; Reidemeister I/II/III site detection and
  greedy simplification heuristics on PD codes
- **Random generators** — seeded random braids, DT/Gauss samples, ensembles
- **CLI** — `knottidy` for Jones, convert, random, simplify, catalogue

## Install

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -e ".[dev]"
```

Or with pip:

```bash
pip install -e ".[dev]"
```

## Quick start

```python
from knottidy import Braid, PDCode, jones_polynomial, jones_trefoil, jones_unknot

# Unknot
assert jones_polynomial(PDCode.unknot()) == jones_unknot()

# Right-handed trefoil as σ₁³ in B₂
trefoil = Braid.trefoil()
V = jones_polynomial(trefoil)
assert V == jones_trefoil()          # t + t³ − t⁴
print(V)                             # t^4 terms...

# Figure-eight braid closure
b41 = Braid.figure_eight()           # σ₁ σ₂⁻¹ σ₁ σ₂⁻¹
print(b41.closure_pd())
print(b41.num_components())          # 1

# Mirror ↔ invert Jones exponents (left trefoil)
print(jones_polynomial(trefoil.mirror()))
```

### CLI

```bash
knottidy jones --trefoil
knottidy jones --braid "1 1 1"
knottidy braid info "1 -2 1 -2"
knottidy convert --braid "1 1 1" --to pd
knottidy random --strands 3 --length 10 --seed 42 --knot
```

## Representations

| Object | Meaning |
|--------|---------|
| `Braid(word, strands)` | Word in ``B_n``; letters are nonzero ints ``±i`` for ``σ_i^{±1}`` |
| `PDCode` | List of crossings ``(a,b,c,d)`` — CCW from under-in; under ``a→c``, over ``b→d`` |
| `DTCode` | Even signed integers, absolutes ``{2,4,…,2n}`` |
| `GaussCode` | Sequence of ``±k`` with each ``k`` once over and once under |
| `LaurentPolynomial` | Exact sparse Laurent series used for bracket / Jones |

## Jones normalisation

The Kauffman bracket is evaluated by state sum with

- ``⟨○⟩ = 1`` (via ``δ^{c−1}``, ``δ = −A² − A^{−2}``),
- ``⟨crossing⟩ = A⟨A-smooth⟩ + A^{−1}⟨B-smooth⟩``.

Then ``V(L) = (−A³)^{−w(L)} ⟨L⟩`` and ``A = t^{−1/4}``. For the right-handed
trefoil one obtains ``V = t + t³ − t⁴``; for the unknot, ``V = 1``.

## Tests

```bash
.venv/bin/pytest -q
```

## License

MIT — see [LICENSE](LICENSE).
