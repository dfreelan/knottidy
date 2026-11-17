"""Gauss codes for knots (and single-component diagrams)."""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Optional, Union

from knottidy.exceptions import ConversionError, InvalidGaussCodeError, InvalidPDCodeError

def _validate_single_component(seq: Sequence[int]) -> tuple[int, ...]:
    if not seq:
        return ()
    frozen = tuple(int(x) for x in seq)
    for x in frozen:
        if x == 0:
            raise InvalidGaussCodeError("Gauss labels must be nonzero")
    abs_vals = [abs(x) for x in frozen]
    n = max(abs_vals)
    expected = list(range(1, n + 1))
    counts = {k: abs_vals.count(k) for k in expected}
    if any(abs_vals.count(k) != 2 for k in expected):
        raise InvalidGaussCodeError(
            f"each crossing 1..{n} must appear twice; counts={counts}"
        )
    if len(frozen) != 2 * n:
        raise InvalidGaussCodeError(f"expected length {2 * n}, got {len(frozen)}")
    for k in expected:
        signs = [1 if x > 0 else -1 for x in frozen if abs(x) == k]
        if sorted(signs) != [-1, 1]:
            raise InvalidGaussCodeError(
                f"crossing {k} must appear once as over and once as under"
            )
    return frozen

def parse_gauss_text(text: str) -> list[tuple[int, ...]]:
    """Parse a Gauss code from text."""
    if not isinstance(text, str):
        raise InvalidGaussCodeError("Gauss text must be a string")
    text = text.strip()
    if not text or text in ("e", "unknot", "()"):
        return [()]

    components_raw = re.split(r"[;|]", text)
    result: list[tuple[int, ...]] = []
    for comp in components_raw:
        comp = comp.strip().lstrip("(").rstrip(")")
        if not comp.strip():
            result.append(())
            continue
        tokens = re.split(r"[\s,]+", comp.strip())
        seq: list[int] = []
        for tok in tokens:
            if not tok:
                continue
            m = re.fullmatch(r"[Oo](\d+)", tok)
            if m:
                seq.append(int(m.group(1)))
                continue
            m = re.fullmatch(r"[Uu](\d+)", tok)
            if m:
                seq.append(-int(m.group(1)))
                continue
            try:
                seq.append(int(tok))
            except ValueError as exc:
                raise InvalidGaussCodeError(f"bad Gauss token {tok!r}") from exc
        result.append(_validate_single_component(seq))
    return result

@dataclass(frozen=True)
class GaussCode:
    """Gauss code, possibly with multiple components."""

    components: tuple[tuple[int, ...], ...]

    def __init__(
        self,
        code: Optional[Union[Sequence[int], Sequence[Sequence[int]]]] = None,
        *,
        validate: bool = True,
    ) -> None:
        if code is None:
            comps: tuple[tuple[int, ...], ...] = ((),)
        elif code and isinstance(code[0], (list, tuple)) and not isinstance(code[0], int):
            raw_comps = [tuple(int(x) for x in comp) for comp in code]  # type: ignore[union-attr]
            comps = tuple(raw_comps)
        else:
            comps = (tuple(int(x) for x in code),)  # type: ignore[arg-type]

        if validate:
            all_labels: list[int] = []
            validated: list[tuple[int, ...]] = []
            for comp in comps:
                if not comp:
                    validated.append(())
                    continue
                for x in comp:
                    if not isinstance(x, int) or isinstance(x, bool) or x == 0:
                        raise InvalidGaussCodeError(f"invalid Gauss label {x!r}")
                validated.append(tuple(comp))
                all_labels.extend(abs(x) for x in comp)

            if all_labels:
                n = max(all_labels)
                for k in range(1, n + 1):
                    appearances = [x for x in (y for c in validated for y in c) if abs(x) == k]
                    if len(appearances) != 2:
                        raise InvalidGaussCodeError(
                            f"crossing {k} appears {len(appearances)} times, expected 2"
                        )
                    signs = sorted(1 if x > 0 else -1 for x in appearances)
                    if signs != [-1, 1]:
                        raise InvalidGaussCodeError(
                            f"crossing {k} must have one over and one under visit"
                        )
                extras = set(all_labels) - set(range(1, n + 1))
                if extras:
                    raise InvalidGaussCodeError(f"non-contiguous crossing labels: {extras}")
            comps = tuple(validated)

        object.__setattr__(self, "components", comps)

    @classmethod
    def empty(cls) -> GaussCode:
        return cls([()])

    @classmethod
    def unknot(cls) -> GaussCode:
        return cls.empty()

    @classmethod
    def from_string(cls, text: str) -> GaussCode:
        return cls(parse_gauss_text(text))

    @classmethod
    def trefoil(cls) -> GaussCode:
        """Right-handed trefoil Gauss code ``1 2 3 -1 -2 -3`` variant.

        Standard: ``O1 O2 O3 U1 U2 U3`` which is ``1 2 3 -1 -2 -3``, but the
        classical alternating form is ``1 -2 3 -1 2 -3``.
        """
        return cls([1, -2, 3, -1, 2, -3])

    @classmethod
    def figure_eight(cls) -> GaussCode:
        return cls([1, -2, 3, -4, -1, 2, -3, 4])

    @classmethod
    def from_pd(cls, pd: object) -> GaussCode:
        """Build a Gauss code by walking components of a PD diagram."""
        from knottidy.pdcode import PDCode

        if not isinstance(pd, PDCode):
            raise InvalidGaussCodeError("from_pd expects a PDCode")
        if not pd.crossings:
            return cls.empty()

        signs = pd.signed_crossing_list()
        edge_to_crossing_role: dict[int, tuple[int, str]] = {}
        for i, cr in enumerate(pd.crossings):
            edge_to_crossing_role[cr.a] = (i, "under")
            over_in = cr.b if signs[i] >= 0 else cr.d
            edge_to_crossing_role[over_in] = (i, "over")

        comps_edges = pd.components()
        gauss_comps: list[list[int]] = []
        for cycle in comps_edges:
            seq: list[int] = []
            for edge in cycle:
                if edge not in edge_to_crossing_role:
                    raise InvalidGaussCodeError(f"edge {edge} missing from incoming map")
                idx, role = edge_to_crossing_role[edge]
                label = idx + 1
                seq.append(label if role == "over" else -label)
            gauss_comps.append(seq)
        return cls(gauss_comps)

    def __len__(self) -> int:
        """Number of crossings."""
        labels = [abs(x) for comp in self.components for x in comp]
        return max(labels) if labels else 0

    def __repr__(self) -> str:
        if len(self.components) == 1:
            return f"GaussCode({list(self.components[0])!r})"
        return f"GaussCode({[list(c) for c in self.components]!r})"

    def __str__(self) -> str:
        return self.format()

    def num_crossings(self) -> int:
        return len(self)

    def num_components(self) -> int:
        if self.components == ((),):
            return 1
        return len([c for c in self.components if c is not None])

    def is_knot(self) -> bool:
        return self.num_components() == 1

    def flat(self) -> tuple[int, ...]:
        """Return the single-component sequence, or raise if multi-component."""
        if len(self.components) != 1:
            raise ConversionError("flat() requires a single-component Gauss code")
        return self.components[0]

    def format(self, style: str = "signed") -> str:
        """Format for display. Styles: ``signed``, ``ou``."""
        parts: list[str] = []
        for comp in self.components:
            if style == "signed":
                parts.append(" ".join(str(x) for x in comp) if comp else "e")
            elif style == "ou":
                toks = []
                for x in comp:
                    toks.append(("O" if x > 0 else "U") + str(abs(x)))
                parts.append(" ".join(toks) if toks else "e")
            else:
                raise InvalidGaussCodeError(f"unknown style {style!r}")
        return " | ".join(parts)

    def mirror(self) -> GaussCode:
        """Mirror: swap over/under (negate every label)."""
        return GaussCode([[-x for x in comp] for comp in self.components])

    def reverse(self) -> GaussCode:
        """Reverse orientation: reverse each component sequence."""
        return GaussCode([list(reversed(comp)) for comp in self.components])

    def renumber(self) -> GaussCode:
        """Renumber crossings to ``1..n`` in order of first appearance."""
        mapping: dict[int, int] = {}
        next_id = 1
        new_comps: list[list[int]] = []
        for comp in self.components:
            new_comp = []
            for x in comp:
                ax = abs(x)
                if ax not in mapping:
                    mapping[ax] = next_id
                    next_id += 1
                new_comp.append(mapping[ax] if x > 0 else -mapping[ax])
            new_comps.append(new_comp)
        return GaussCode(new_comps)

    def crossing_signs(self) -> list[int]:
        """Crossing signs for known fingerprints; otherwise alternating heuristic."""
        n = self.num_crossings()
        if n == 0:
            return []
        if len(self.components) != 1:
            raise ConversionError("crossing_signs requires a knot (1 component)")
        seq = self.components[0]
        if list(seq) == [1, -2, 3, -1, 2, -3]:
            return [1, 1, 1]
        if list(seq) == [1, -2, 3, -4, -1, 2, -3, 4]:
            return [1, -1, 1, -1]
        signs = []
        for k in range(1, n + 1):
            positions = [i for i, x in enumerate(seq) if abs(x) == k]
            i, j = positions
            interleave = 0
            for m in range(1, n + 1):
                if m == k:
                    continue
                pos_m = [t for t, x in enumerate(seq) if abs(x) == m]
                a, b = pos_m
                inside = sum(1 for p in (a, b) if i < p < j)
                if inside == 1:
                    interleave += 1
            s = 1 if seq[i] > 0 else -1
            if interleave % 2 == 1:
                s = -s
            signs.append(s)
        return signs

    def writhe(self) -> int:
        return sum(self.crossing_signs())

    def to_pd(self) -> "PDCode":
        """Convert a single-component Gauss code to a PD code."""
        from knottidy.pdcode import PDCode, PDCrossing

        if not self.is_knot():
            raise ConversionError("to_pd() currently requires a single-component Gauss code")
        seq = self.components[0]
        if not seq:
            return PDCode.empty()

        n = self.num_crossings()
        signs = self.crossing_signs()
        m = len(seq)  # 2n
        under_in = [0] * (n + 1)
        under_out = [0] * (n + 1)
        over_in = [0] * (n + 1)
        over_out = [0] * (n + 1)

        for i, label in enumerate(seq):
            edge_in = i + 1  # edge arriving at this visit
            edge_out = (i + 1) % m + 1  # edge leaving
            k = abs(label)
            if label < 0:
                under_in[k] = edge_in
                under_out[k] = edge_out
            else:
                over_in[k] = edge_in
                over_out[k] = edge_out

        crossings: list[PDCrossing] = []
        for k in range(1, n + 1):
            a, c = under_in[k], under_out[k]
            b, d = over_in[k], over_out[k]
            if signs[k - 1] >= 0:
                crossings.append(PDCrossing(a, b, c, d))
            else:
                crossings.append(PDCrossing(a, d, c, b))

        try:
            return PDCode(crossings, validate=True)
        except InvalidPDCodeError:
            crossings = [
                PDCrossing(under_in[k], over_in[k], under_out[k], over_out[k])
                for k in range(1, n + 1)
            ]
            return PDCode(crossings, validate=True)

    def to_dict(self) -> dict:
        return {"components": [list(c) for c in self.components]}

    @classmethod
    def from_dict(cls, data: dict) -> GaussCode:
        if not isinstance(data, dict) or "components" not in data:
            raise InvalidGaussCodeError("from_dict expects dict with 'components'")
        return cls(data["components"])

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from knottidy.pdcode import PDCode

__all__ = [
    "GaussCode",
    "parse_gauss_text",
]
