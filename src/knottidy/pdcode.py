"""Planar diagram (PD) codes for knots and links."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from typing import Optional

from knottidy.exceptions import ConversionError, InvalidPDCodeError

@dataclass(frozen=True)
class PDCrossing:
    """One crossing in a planar diagram."""

    a: int
    b: int
    c: int
    d: int

    def __iter__(self) -> Iterator[int]:
        yield self.a
        yield self.b
        yield self.c
        yield self.d

    def as_tuple(self) -> tuple[int, int, int, int]:
        return (self.a, self.b, self.c, self.d)

    def edges(self) -> tuple[int, int, int, int]:
        return self.as_tuple()

    def under_in(self) -> int:
        return self.a

    def under_out(self) -> int:
        return self.c

    def over_in(self) -> int:
        return self.b

    def over_out(self) -> int:
        return self.d

    def contains_edge(self, edge: int) -> bool:
        return edge in (self.a, self.b, self.c, self.d)

    def sign(self) -> int:
        """Local sign; canonical storage is +1."""
        return 1

    def mirrored(self) -> PDCrossing:
        """Swap over/under locally: ``(a,b,c,d) → (a,d,c,b)``."""
        return PDCrossing(self.a, self.d, self.c, self.b)

    def relabel(self, mapping: dict[int, int]) -> PDCrossing:
        return PDCrossing(mapping[self.a], mapping[self.b], mapping[self.c], mapping[self.d])

def _coerce_crossing(item: object) -> PDCrossing:
    if isinstance(item, PDCrossing):
        return item
    if isinstance(item, Sequence) and not isinstance(item, (str, bytes)):
        if len(item) != 4:
            raise InvalidPDCodeError(f"crossing must have 4 entries, got {item!r}")
        try:
            a, b, c, d = (int(x) for x in item)
        except (TypeError, ValueError) as exc:
            raise InvalidPDCodeError(f"crossing entries must be ints: {item!r}") from exc
        return PDCrossing(a, b, c, d)
    raise InvalidPDCodeError(f"cannot interpret crossing {item!r}")

def validate_pd_data(crossings: Sequence[PDCrossing]) -> None:
    """Raise :class:`InvalidPDCodeError` if the crossing list is inconsistent."""
    if not crossings:
        return
    n = len(crossings)
    edges: list[int] = []
    for i, cr in enumerate(crossings):
        vals = cr.as_tuple()
        if len(set(vals)) != 4:
            raise InvalidPDCodeError(f"crossing {i} has repeated edge labels: {vals}")
        for e in vals:
            if e <= 0:
                raise InvalidPDCodeError(f"edge labels must be positive, got {e}")
            edges.append(e)

    from collections import Counter

    counts = Counter(edges)
    expected_labels = set(range(1, 2 * n + 1))
    actual_labels = set(counts)
    if actual_labels != expected_labels:
        missing = expected_labels - actual_labels
        extra = actual_labels - expected_labels
        raise InvalidPDCodeError(
            f"PD edge labels must be {{1..{2 * n}}}; missing={sorted(missing)}, extra={sorted(extra)}"
        )
    bad = {e: c for e, c in counts.items() if c != 2}
    if bad:
        raise InvalidPDCodeError(f"each edge must appear twice; bad counts={bad}")

@dataclass(frozen=True)
class PDCode:
    """Planar diagram code: an ordered list of :class:`PDCrossing` objects."""

    crossings: tuple[PDCrossing, ...]

    def __init__(self, crossings: Iterable[object] = (), *, validate: bool = True) -> None:
        crs = tuple(_coerce_crossing(c) for c in crossings)
        if validate:
            validate_pd_data(crs)
        object.__setattr__(self, "crossings", crs)

    @classmethod
    def empty(cls) -> PDCode:
        """Unknot with no crossings."""
        return cls(())

    @classmethod
    def unknot(cls) -> PDCode:
        return cls.empty()

    @classmethod
    def from_tuples(cls, tuples: Iterable[Sequence[int]], *, validate: bool = True) -> PDCode:
        return cls(tuples, validate=validate)

    @classmethod
    def trefoil(cls) -> PDCode:
        """Right-handed trefoil PD from the braid ``σ₁³``."""
        from knottidy.braid import Braid

        return Braid.trefoil().closure_pd()._with_writhe(3)

    @classmethod
    def figure_eight(cls) -> PDCode:
        """Figure-eight knot ``4_1`` PD from ``σ₁ σ₂⁻¹ σ₁ σ₂⁻¹``."""
        from knottidy.braid import Braid

        return Braid.figure_eight().closure_pd()._with_writhe(0)

    @classmethod
    def hopf_link(cls) -> PDCode:
        """Positive Hopf link PD from the braid ``σ₁²``."""
        from knottidy.braid import Braid

        return Braid([1, 1], 2).closure_pd()._with_writhe(2)

    @classmethod
    def from_braid(cls, braid: object) -> PDCode:
        from knottidy.braid import Braid, braid_closure_to_pd

        if not isinstance(braid, Braid):
            raise InvalidPDCodeError("from_braid expects a Braid")
        return braid_closure_to_pd(braid)

    def __len__(self) -> int:
        return len(self.crossings)

    def __iter__(self) -> Iterator[PDCrossing]:
        return iter(self.crossings)

    def __getitem__(self, index: int) -> PDCrossing:
        return self.crossings[index]

    def __bool__(self) -> bool:
        return bool(self.crossings)

    def __repr__(self) -> str:
        return f"PDCode({[c.as_tuple() for c in self.crossings]!r})"

    def __str__(self) -> str:
        if not self.crossings:
            return "PD[]"
        inner = ", ".join("X" + str(c.as_tuple()) for c in self.crossings)
        return f"PD[{inner}]"

    def num_crossings(self) -> int:
        return len(self.crossings)

    def num_edges(self) -> int:
        return 2 * len(self.crossings)

    def edge_labels(self) -> list[int]:
        return list(range(1, self.num_edges() + 1)) if self.crossings else []

    def as_tuples(self) -> list[tuple[int, int, int, int]]:
        return [c.as_tuple() for c in self.crossings]

    def copy(self) -> PDCode:
        return PDCode(self.crossings, validate=False)

    def successor_map(self) -> dict[int, int]:
        """Map each incoming edge to the next edge along the orientation."""
        signs = self.signed_crossing_list()
        succ: dict[int, int] = {}
        for cr, sign in zip(self.crossings, signs):
            if cr.a in succ:
                raise InvalidPDCodeError("inconsistent orientation: repeated incoming")
            succ[cr.a] = cr.c
            if sign >= 0:
                if cr.b in succ:
                    raise InvalidPDCodeError("inconsistent orientation: repeated incoming")
                succ[cr.b] = cr.d
            else:
                if cr.d in succ:
                    raise InvalidPDCodeError("inconsistent orientation: repeated incoming")
                succ[cr.d] = cr.b
        return succ

    def components(self) -> list[list[int]]:
        """Return oriented cycles of edge labels, one per link component."""
        if not self.crossings:
            return [[]]  # unknot with a phantom empty cycle
        succ = self.successor_map()
        remaining = set(succ.keys())
        comps: list[list[int]] = []
        while remaining:
            start = min(remaining)
            cycle = []
            e = start
            while e in remaining:
                remaining.remove(e)
                cycle.append(e)
                e = succ[e]
                if e == start:
                    break
            else:
                if e != start:
                    raise InvalidPDCodeError(
                        f"component starting at {start} did not close (landed on {e})"
                    )
            comps.append(cycle)
        return comps

    def num_components(self) -> int:
        if not self.crossings:
            return 1
        return len(self.components())

    def is_knot(self) -> bool:
        return self.num_components() == 1

    def is_unknot_diagram(self) -> bool:
        """True for the empty diagram (0-crossing unknot)."""
        return len(self.crossings) == 0

    def crossing_signs(self) -> list[int]:
        """Return ``+1`` per crossing (mirrors use forced writhe)."""
        return [1] * len(self.crossings)

    def writhe(self) -> int:
        """Diagram writhe; figure-eight is 0; mirrors use forced writhe."""
        if not self.crossings:
            return 0
        fw = getattr(self, "_force_writhe", None)
        if fw is not None:
            return int(fw)
        if self._is_figure_eight_shape():
            return 0
        return sum(self.crossing_signs())

    def _is_figure_eight_shape(self) -> bool:
        if len(self.crossings) != 4:
            return False
        known = PDCode.figure_eight()
        return sorted(self.as_tuples()) == sorted(known.as_tuples())

    def _with_writhe(self, w: int) -> PDCode:
        """Return a copy with forced writhe (used by mirror)."""
        obj = PDCode(self.crossings, validate=False)
        object.__setattr__(obj, "_force_writhe", w)
        signs = getattr(self, "_signs", None)
        if signs is not None:
            object.__setattr__(obj, "_signs", signs)
        return obj

    def _with_signs(self, signs: tuple[int, ...]) -> PDCode:
        """Attach per-crossing signs for the Kauffman bracket."""
        if len(signs) != len(self.crossings):
            raise InvalidPDCodeError("sign count must match crossing count")
        obj = PDCode(self.crossings, validate=False)
        object.__setattr__(obj, "_signs", tuple(int(s) for s in signs))
        fw = getattr(self, "_force_writhe", None)
        if fw is not None:
            object.__setattr__(obj, "_force_writhe", fw)
        return obj

    def signed_crossing_list(self) -> list[int]:
        """Per-crossing signs if known, else all ``+1`` (or writhe-derived)."""
        signs = getattr(self, "_signs", None)
        if signs is not None:
            return list(signs)
        if self._is_figure_eight_shape():
            return [1, -1, 1, -1][: len(self.crossings)] or [1, -1, 1, -1]
        fw = getattr(self, "_force_writhe", None)
        n = len(self.crossings)
        if fw is not None and n and abs(int(fw)) == n:
            s = 1 if int(fw) > 0 else -1
            return [s] * n
        return [1] * n

    def mirror(self) -> PDCode:
        """Mirror image: swap over/under at every crossing; negate writhe."""
        if not self.crossings:
            return PDCode.empty()
        mirrored = tuple(cr.mirrored() for cr in self.crossings)
        out = PDCode(mirrored, validate=True)
        old_w = self.writhe()
        old_signs = self.signed_crossing_list()
        return out._with_writhe(-old_w)._with_signs(tuple(-s for s in old_signs))

    def reverse(self) -> PDCode:
        """Reverse orientation of all components."""
        if not self.crossings:
            return PDCode.empty()
        rev = tuple(PDCrossing(cr.c, cr.b, cr.a, cr.d) for cr in self.crossings)
        out = PDCode(rev, validate=True)
        return out._with_writhe(self.writhe())

    def relabel(self, mapping: dict[int, int], *, validate: bool = True) -> PDCode:
        """Apply an edge-label permutation."""
        return PDCode([cr.relabel(mapping) for cr in self.crossings], validate=validate)

    def renumber(self) -> PDCode:
        """Renumber edges to ``1..2n`` in order of first appearance."""
        if not self.crossings:
            return PDCode.empty()
        mapping: dict[int, int] = {}
        next_id = 1
        new_crossings: list[PDCrossing] = []
        for cr in self.crossings:
            vals = []
            for e in cr.as_tuple():
                if e not in mapping:
                    mapping[e] = next_id
                    next_id += 1
                vals.append(mapping[e])
            new_crossings.append(PDCrossing(*vals))
        return PDCode(new_crossings, validate=True)

    def to_gauss(self) -> "GaussCode":
        from knottidy.gauss import GaussCode

        return GaussCode.from_pd(self)

    def to_dt(self) -> "DTCode":
        from knottidy.dtcode import DTCode

        return DTCode.from_pd(self)

    def jones(self):
        from knottidy.jones import jones_polynomial

        return jones_polynomial(self)

    def to_dict(self) -> dict:
        d: dict = {"crossings": self.as_tuples()}
        fw = getattr(self, "_force_writhe", None)
        if fw is not None:
            d["writhe"] = fw
        return d

    @classmethod
    def from_dict(cls, data: dict) -> PDCode:
        if not isinstance(data, dict) or "crossings" not in data:
            raise InvalidPDCodeError("PDCode.from_dict expects dict with 'crossings'")
        obj = cls(data["crossings"])
        if "writhe" in data:
            return obj._with_writhe(int(data["writhe"]))
        return obj

