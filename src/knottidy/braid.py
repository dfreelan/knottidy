"""Braid group words and related operations."""

from __future__ import annotations

import re
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from typing import Optional

from knottidy.exceptions import ConversionError, InvalidBraidError

_TOKEN_RE = re.compile(
    r"""
    (?:s|sigma|σ)?
    _?
    (-?\d+)
    (?:\^(-?\d+))?
    """,
    re.IGNORECASE | re.VERBOSE,
)


def _validate_strands(strands: int) -> int:
    if not isinstance(strands, int) or isinstance(strands, bool):
        raise InvalidBraidError(f"strands must be a positive int, got {strands!r}")
    if strands < 1:
        raise InvalidBraidError(f"strands must be >= 1, got {strands}")
    return strands


def _normalize_word(word: Iterable[int], strands: int) -> tuple[int, ...]:
    """Validate and freeze a braid word for ``B_strands``."""
    strands = _validate_strands(strands)
    max_gen = strands - 1
    out: list[int] = []
    for raw in word:
        if not isinstance(raw, int) or isinstance(raw, bool):
            raise InvalidBraidError(f"braid generators must be ints, got {raw!r}")
        if raw == 0:
            raise InvalidBraidError("generator index 0 is not allowed")
        if abs(raw) > max_gen:
            raise InvalidBraidError(
                f"generator |{raw}| exceeds max index {max_gen} for B_{strands}"
            )
        out.append(int(raw))
    return tuple(out)

def free_reduce(word: Sequence[int]) -> tuple[int, ...]:
    """Cancel adjacent inverse pairs ``σ_i σ_i^{-1}`` until none remain.

    This is the free-group reduction on the Artin generators; it does not
    apply the braid relation ``σ_i σ_{i+1} σ_i = σ_{i+1} σ_i σ_{i+1}``.
    """
    stack: list[int] = []
    for g in word:
        if not isinstance(g, int) or isinstance(g, bool) or g == 0:
            raise InvalidBraidError(f"invalid generator {g!r}")
        if stack and stack[-1] == -g:
            stack.pop()
        else:
            stack.append(int(g))
    return tuple(stack)

def parse_braid_word(text: str, strands: Optional[int] = None) -> tuple[tuple[int, ...], int]:
    """Parse a human-readable braid word."""
    if not isinstance(text, str):
        raise InvalidBraidError("braid text must be a string")
    cleaned = text.strip()
    if not cleaned or cleaned in ("1", "e", "id", "ε", "emptyset"):
        n = strands if strands is not None else 1
        return (), _validate_strands(n)

    parts = re.split(r"[\s,;]+", cleaned)
    gens: list[int] = []
    for part in parts:
        if not part:
            continue
        m = _TOKEN_RE.fullmatch(part)
        if not m:
            raise InvalidBraidError(f"cannot parse braid token {part!r}")
        idx = int(m.group(1))
        if idx == 0:
            raise InvalidBraidError("generator index 0 is not allowed")
        power = int(m.group(2)) if m.group(2) is not None else 1
        if power == 0:
            continue
        sign = 1 if idx > 0 else -1
        abs_idx = abs(idx)
        unit = sign * abs_idx
        if power > 0:
            gens.extend([unit] * power)
        else:
            gens.extend([-unit] * (-power))

    if not gens:
        n = strands if strands is not None else 1
        return (), _validate_strands(n)

    inferred = 1 + max(abs(g) for g in gens)
    n = strands if strands is not None else inferred
    if strands is not None and inferred > strands:
        raise InvalidBraidError(
            f"word requires at least {inferred} strands, got strands={strands}"
        )
    return _normalize_word(gens, n), n

def format_braid_word(word: Sequence[int], style: str = "signed") -> str:
    """Format a braid word for display."""
    if style == "signed":
        return " ".join(str(g) for g in word) if word else "e"
    if style == "sigma":
        parts = []
        for g in word:
            if g > 0:
                parts.append(f"s{g}")
            else:
                parts.append(f"s{-g}^-1")
        return " ".join(parts) if parts else "e"
    if style == "compact":
        if not word:
            return "e"
        parts: list[str] = []
        i = 0
        while i < len(word):
            g = word[i]
            j = i + 1
            while j < len(word) and word[j] == g:
                j += 1
            count = j - i
            idx = abs(g)
            if g > 0:
                parts.append(f"s{idx}" if count == 1 else f"s{idx}^{count}")
            else:
                parts.append(f"s{idx}^-{count}")
            i = j
        return " ".join(parts)
    raise InvalidBraidError(f"unknown format style {style!r}")

@dataclass(frozen=True)
class Braid:
    """Element of the Artin braid group ``B_n``."""

    word: tuple[int, ...]
    strands: int
    _reduced: bool = field(default=True, repr=False, compare=False)

    def __init__(
        self,
        word: Iterable[int] = (),
        strands: Optional[int] = None,
        *,
        reduce_on_init: bool = True,
    ) -> None:
        w_list = list(word)
        if strands is None:
            strands = 1 + max((abs(g) for g in w_list), default=0)
            strands = max(strands, 1)
        strands = _validate_strands(strands)
        frozen = _normalize_word(w_list, strands)
        if reduce_on_init:
            frozen = free_reduce(frozen)
        object.__setattr__(self, "word", frozen)
        object.__setattr__(self, "strands", strands)
        object.__setattr__(self, "_reduced", True)

    @classmethod
    def identity(cls, strands: int) -> Braid:
        """Return the identity braid on ``strands`` strands."""
        return cls((), strands, reduce_on_init=False)

    @classmethod
    def generator(cls, index: int, strands: int, power: int = 1) -> Braid:
        """Return ``σ_index^power`` in ``B_strands``."""
        strands = _validate_strands(strands)
        if index < 1 or index >= strands:
            raise InvalidBraidError(f"generator index {index} invalid for B_{strands}")
        if power == 0:
            return cls.identity(strands)
        unit = index if power > 0 else -index
        return cls([unit] * abs(power), strands)

    @classmethod
    def from_string(cls, text: str, strands: Optional[int] = None) -> Braid:
        """Parse and construct a braid from text (see :func:`parse_braid_word`)."""
        word, n = parse_braid_word(text, strands)
        return cls(word, n)

    @classmethod
    def trefoil(cls) -> Braid:
        """Standard positive trefoil braid ``σ_1^3`` in ``B_2``."""
        return cls([1, 1, 1], 2)

    @classmethod
    def figure_eight(cls) -> Braid:
        """Figure-eight knot as braid ``σ_1 σ_2^{-1} σ_1 σ_2^{-1}`` in ``B_3``."""
        return cls([1, -2, 1, -2], 3)

    def __len__(self) -> int:
        return len(self.word)

    def __iter__(self) -> Iterator[int]:
        return iter(self.word)

    def __bool__(self) -> bool:
        return bool(self.word)

    def __repr__(self) -> str:
        return f"Braid({list(self.word)!r}, strands={self.strands})"

    def __str__(self) -> str:
        return f"B_{self.strands}<{format_braid_word(self.word)}>"

    def is_identity(self) -> bool:
        return len(self.word) == 0

    def copy(self) -> Braid:
        return Braid(self.word, self.strands, reduce_on_init=False)

    def __mul__(self, other: object) -> Braid:
        if not isinstance(other, Braid):
            return NotImplemented
        strands = max(self.strands, other.strands)
        return Braid(self.word + other.word, strands)

    def __pow__(self, exponent: int) -> Braid:
        if not isinstance(exponent, int) or isinstance(exponent, bool):
            return NotImplemented
        if exponent == 0:
            return Braid.identity(self.strands)
        if exponent > 0:
            word = self.word * exponent
            return Braid(word, self.strands)
        inv = self.inverse()
        return Braid(inv.word * (-exponent), self.strands)

    def inverse(self) -> Braid:
        """Return the group inverse."""
        return Braid([-g for g in reversed(self.word)], self.strands, reduce_on_init=False)

    def __invert__(self) -> Braid:
        return self.inverse()

    def conjugate(self, conjugator: Braid) -> Braid:
        """Return ``conjugator * self * conjugator^{-1}``."""
        if not isinstance(conjugator, Braid):
            raise InvalidBraidError("conjugator must be a Braid")
        strands = max(self.strands, conjugator.strands)
        return Braid(
            conjugator.word + self.word + conjugator.inverse().word,
            strands,
        )

    def commute_form(self) -> "Braid":
        """Heuristic shortener (free reduction only)."""
        return self.free_reduce()

    def free_reduce(self) -> Braid:
        """Return a freely reduced copy (usually already reduced)."""
        return Braid(free_reduce(self.word), self.strands, reduce_on_init=False)

    def mirror(self) -> Braid:
        """Mirror image: invert every generator sign (``σ_i ↔ σ_i^{-1}``)."""
        return Braid([-g for g in self.word], self.strands, reduce_on_init=False)

    def reverse(self) -> Braid:
        """Reverse the word (time-reversal of the braid)."""
        return Braid(tuple(reversed(self.word)), self.strands, reduce_on_init=False)

    def add_strands(self, extra: int = 1, side: str = "right") -> Braid:
        """Stabilize by adding trivial strands.

        ``side='right'`` appends strands on the right (generators unchanged).
        ``side='left'`` prepends strands and shifts all generator indices up.
        """
        if extra < 0:
            raise InvalidBraidError("extra strands must be non-negative")
        if side == "right":
            return Braid(self.word, self.strands + extra, reduce_on_init=False)
        if side == "left":
            shifted = [g + extra if g > 0 else g - extra for g in self.word]
            return Braid(shifted, self.strands + extra, reduce_on_init=False)
        raise InvalidBraidError(f"side must be 'left' or 'right', got {side!r}")

    def stabilize(self, positive: bool = True, side: str = "right") -> Braid:
        """Markov stabilization: add a strand and a new extremal generator.

        Positive right stabilization of ``β ∈ B_n`` yields ``β σ_n ∈ B_{n+1}``.
        """
        if side not in ("left", "right"):
            raise InvalidBraidError(f"side must be 'left' or 'right', got {side!r}")
        if side == "right":
            new_strands = self.strands + 1
            gen = self.strands if positive else -self.strands
            return Braid(self.word + (gen,), new_strands)
        shifted = [g + 1 if g > 0 else g - 1 for g in self.word]
        gen = 1 if positive else -1
        return Braid((gen,) + tuple(shifted), self.strands + 1)

    def destabilize(self) -> Braid:
        """Attempt Markov destabilization when the last strand is locked by one generator.

        Requires that the highest-index generator ``±(n-1)`` appears exactly once
        and that removing it leaves a braid on ``n-1`` strands.
        """
        if self.strands < 2:
            raise InvalidBraidError("cannot destabilize B_1")
        top = self.strands - 1
        count = sum(1 for g in self.word if abs(g) == top)
        if count != 1:
            raise InvalidBraidError(
                f"destabilization requires exactly one ±σ_{top}, found {count}"
            )
        new_word = tuple(g for g in self.word if abs(g) != top)
        if new_word and max(abs(g) for g in new_word) >= top:
            raise InvalidBraidError("remaining word not in B_{n-1}")
        return Braid(new_word, self.strands - 1, reduce_on_init=True)

    def can_destabilize(self) -> bool:
        """Return whether :meth:`destabilize` would succeed."""
        try:
            self.destabilize()
            return True
        except InvalidBraidError:
            return False

    def writhe(self) -> int:
        """Algebraic writhe: ``#positive generators − #negative generators``."""
        return sum(1 if g > 0 else -1 for g in self.word)

    def exponent_sum(self) -> int:
        """Alias for :meth:`writhe`."""
        return self.writhe()

    def permutation(self) -> tuple[int, ...]:
        """Induced permutation of strands as a tuple ``π`` with ``π[i]`` = image of ``i``.

        Strands are labelled ``0 .. n-1``. Crossing ``σ_k`` swaps strands
        currently in positions ``k-1`` and ``k``.
        """
        n = self.strands
        pos = list(range(n))
        for g in self.word:
            k = abs(g) - 1
            pos[k], pos[k + 1] = pos[k + 1], pos[k]
        pi = [0] * n
        for position, original in enumerate(pos):
            pi[original] = position
        return tuple(pi)

    def is_pure(self) -> bool:
        """Return true if the braid induces the identity permutation."""
        return self.permutation() == tuple(range(self.strands))

    def strand_permutation_cycles(self) -> list[tuple[int, ...]]:
        """Return cycles of the induced strand permutation (0-based labels)."""
        pi = self.permutation()
        n = len(pi)
        seen = [False] * n
        cycles: list[tuple[int, ...]] = []
        for i in range(n):
            if seen[i]:
                continue
            cycle = []
            j = i
            while not seen[j]:
                seen[j] = True
                cycle.append(j)
                j = pi[j]
            cycles.append(tuple(cycle))
        return cycles

    def num_components(self) -> int:
        """Number of link components in the braid closure."""
        return len(self.strand_permutation_cycles())

    def is_knot_closure(self) -> bool:
        """True iff the closure is a single-component knot."""
        return self.num_components() == 1

    def closure_pd(self) -> "PDCode":
        """Return a planar diagram code for the braid closure."""
        from knottidy.pdcode import PDCode

        return braid_closure_to_pd(self)

    def jones(self):
        """Jones polynomial of the braid closure (Laurent polynomial in ``t``)."""
        from knottidy.jones import jones_polynomial

        return jones_polynomial(self.closure_pd())

    def to_dict(self) -> dict:
        """Serialise to a plain dictionary."""
        return {"strands": self.strands, "word": list(self.word)}

    @classmethod
    def from_dict(cls, data: dict) -> Braid:
        """Deserialise from :meth:`to_dict` output."""
        if not isinstance(data, dict):
            raise InvalidBraidError("Braid.from_dict expects a dict")
        try:
            return cls(data["word"], int(data["strands"]))
        except (KeyError, TypeError, ValueError) as exc:
            raise InvalidBraidError(f"invalid braid dict: {exc}") from exc

    def format(self, style: str = "signed") -> str:
        return format_braid_word(self.word, style=style)

def braid_closure_to_pd(braid: Braid) -> "PDCode":
    """Convert a braid to a PD code via geometric closure."""
    from knottidy.pdcode import PDCode, PDCrossing

    n = braid.strands
    word = braid.word
    if n == 1 and not word:
        return PDCode([])

    next_edge = 1
    top_edges = [0] * n
    active = [0] * n
    for i in range(n):
        top_edges[i] = next_edge
        active[i] = next_edge
        next_edge += 1

    crossings: list[PDCrossing] = []

    for g in word:
        k = abs(g) - 1  # left position index of the crossing
        left_in = active[k]
        right_in = active[k + 1]
        left_out = next_edge
        next_edge += 1
        right_out = next_edge
        next_edge += 1

        if g > 0:
            # (under_in, over_in, under_out, over_out) for a positive crossing.
            a, b, c, d = right_in, left_in, left_out, right_out
        else:
            # Mirrored cyclic order of the positive pattern (A/B regions swapped).
            # Equivalent to PDCrossing(*positive).mirrored().
            a, b, c, d = right_in, right_out, left_out, left_in

        crossings.append(PDCrossing(a, b, c, d))
        active[k] = left_out
        active[k + 1] = right_out

    parent: dict[int, int] = {}

    def find(x: int) -> int:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            if ra < rb:
                parent[rb] = ra
            else:
                parent[ra] = rb

    for i in range(n):
        union(active[i], top_edges[i])

    remapped_crossings: list[tuple[int, int, int, int]] = []
    used: set[int] = set()
    for cr in crossings:
        vals = (find(cr.a), find(cr.b), find(cr.c), find(cr.d))
        remapped_crossings.append(vals)
        used.update(vals)

    if not remapped_crossings:
        return PDCode([])

    sorted_ids = sorted(used)
    rename = {old: new for new, old in enumerate(sorted_ids, start=1)}
    final = [
        PDCrossing(rename[a], rename[b], rename[c], rename[d])
        for a, b, c, d in remapped_crossings
    ]
    # Cyclic order encodes handedness for the bracket; signs help successor/Gauss.
    signs = tuple(1 if g > 0 else -1 for g in word)
    return (
        PDCode(final, validate=True)
        ._with_writhe(braid.writhe())
        ._with_signs(signs)
    )

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from knottidy.pdcode import PDCode

__all__ = [
    "Braid",
    "free_reduce",
            "parse_braid_word",
    "format_braid_word",
    "braid_closure_to_pd",
]
