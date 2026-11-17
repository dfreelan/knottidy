"""Dowker–Thistlethwaite (DT) codes."""

from __future__ import annotations

import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Optional

from knottidy.exceptions import ConversionError, InvalidDTCodeError, InvalidPDCodeError

def validate_dt_sequence(seq: Sequence[int]) -> tuple[int, ...]:
    """Validate a DT sequence and return it as a tuple."""
    if not seq:
        return ()
    frozen = tuple(int(x) for x in seq)
    n = len(frozen)
    for x in frozen:
        if x % 2 != 0:
            raise InvalidDTCodeError(f"DT entries must be even, got {x}")
        if x == 0:
            raise InvalidDTCodeError("DT entries must be nonzero")
    abs_vals = sorted(abs(x) for x in frozen)
    expected = list(range(2, 2 * n + 1, 2))
    if abs_vals != expected:
        raise InvalidDTCodeError(
            f"absolute values must be {expected}, got {abs_vals}"
        )
    return frozen

def parse_dt_text(text: str) -> tuple[int, ...]:
    """Parse a DT code from text such as ``4 6 2``, ``DT[4,6,2]``, or ``[4,-6,2]``."""
    if not isinstance(text, str):
        raise InvalidDTCodeError("DT text must be a string")
    cleaned = text.strip()
    if not cleaned or cleaned.lower() in ("e", "unknot", "dt[]", "[]"):
        return ()
    cleaned = re.sub(r"^[Dd][Tt]\s*", "", cleaned)
    cleaned = cleaned.strip().lstrip("[").rstrip("]")
    if not cleaned.strip():
        return ()
    tokens = re.split(r"[\s,;]+", cleaned.strip())
    try:
        vals = [int(t) for t in tokens if t]
    except ValueError as exc:
        raise InvalidDTCodeError(f"cannot parse DT text {text!r}") from exc
    return validate_dt_sequence(vals)

@dataclass(frozen=True)
class DTCode:
    """Dowker–Thistlethwaite code for a knot diagram."""

    code: tuple[int, ...]

    def __init__(self, code: Optional[Iterable[int]] = None, *, validate: bool = True) -> None:
        seq = tuple(code) if code is not None else ()
        if validate:
            seq = validate_dt_sequence(seq)
        else:
            seq = tuple(int(x) for x in seq)
        object.__setattr__(self, "code", seq)

    @classmethod
    def empty(cls) -> DTCode:
        return cls(())

    @classmethod
    def unknot(cls) -> DTCode:
        return cls.empty()

    @classmethod
    def from_string(cls, text: str) -> DTCode:
        return cls(parse_dt_text(text))

    @classmethod
    def trefoil(cls) -> DTCode:
        """DT code of the right-handed trefoil: ``DT[4, 6, 2]``."""
        return cls([4, 6, 2])

    @classmethod
    def figure_eight(cls) -> DTCode:
        """DT code of the figure-eight knot: ``DT[4, 6, 8, 2]``."""
        return cls([4, 6, 8, 2])

    @classmethod
    def from_pd(cls, pd: object) -> DTCode:
        """Compute a DT code from a single-component PD diagram."""
        from knottidy.pdcode import PDCode

        if not isinstance(pd, PDCode):
            raise InvalidDTCodeError("from_pd expects a PDCode")
        if not pd.crossings:
            return cls.empty()
        if pd.num_components() != 1:
            raise ConversionError("DT codes require a single-component diagram")

        return cls._from_pd_rotate(pd)

    @classmethod
    def _from_pd_rotate(cls, pd: object) -> DTCode:
        """Extract a DT code, trying each starting edge as basepoint."""
        from knottidy.pdcode import PDCode

        assert isinstance(pd, PDCode)
        n = len(pd.crossings)
        if n == 0:
            return cls.empty()
        signs = pd.signed_crossing_list()
        under_edges = {cr.a: i for i, cr in enumerate(pd.crossings)}
        over_edges = {
            (cr.b if signs[i] >= 0 else cr.d): i for i, cr in enumerate(pd.crossings)
        }
        component = pd.components()[0]
        if len(component) != 2 * n:
            raise InvalidDTCodeError("component length mismatch with crossing count")

        for shift in range(len(component)):
            rotated = component[shift:] + component[:shift]
            # label -> (crossing, is_over)
            visits: dict[int, tuple[int, bool]] = {}
            ok = True
            for step, edge in enumerate(rotated, start=1):
                if edge in under_edges:
                    visits[step] = (under_edges[edge], False)
                elif edge in over_edges:
                    visits[step] = (over_edges[edge], True)
                else:
                    ok = False
                    break
            if not ok:
                continue
            # Pair labels that share a crossing; require opposite parity.
            by_crossing: dict[int, list[tuple[int, bool]]] = {}
            for lab, (cidx, is_over) in visits.items():
                by_crossing.setdefault(cidx, []).append((lab, is_over))
            if any(len(v) != 2 for v in by_crossing.values()):
                continue
            odd_to_even: dict[int, int] = {}
            good = True
            for cidx, pair in by_crossing.items():
                (l1, o1), (l2, o2) = pair
                if (l1 % 2) == (l2 % 2):
                    good = False
                    break
                if l1 % 2 == 1:
                    odd, even, even_is_over = l1, l2, o2
                else:
                    odd, even, even_is_over = l2, l1, o1
                odd_to_even[odd] = even if even_is_over else -even
            if not good or len(odd_to_even) != n:
                continue
            dt_list = [odd_to_even[odd] for odd in range(1, 2 * n + 1, 2)]
            try:
                return cls(dt_list)
            except InvalidDTCodeError:
                continue
        raise ConversionError("could not extract a DT code from the given PD diagram")

    @classmethod
    def from_gauss(cls, gauss: object) -> DTCode:
        """Convert a knot Gauss code to DT by going through PD."""
        from knottidy.gauss import GaussCode

        if not isinstance(gauss, GaussCode):
            raise InvalidDTCodeError("from_gauss expects a GaussCode")
        return cls.from_pd(gauss.to_pd())

    def __len__(self) -> int:
        return len(self.code)

    def __iter__(self):
        return iter(self.code)

    def __repr__(self) -> str:
        return f"DTCode({list(self.code)!r})"

    def __str__(self) -> str:
        if not self.code:
            return "DT[]"
        return "DT[" + ", ".join(str(x) for x in self.code) + "]"

    def num_crossings(self) -> int:
        return len(self.code)

    def is_empty(self) -> bool:
        return not self.code

    def format(self) -> str:
        return str(self)

    def mirror(self) -> DTCode:
        """Mirror image: negate every entry."""
        return DTCode([-x for x in self.code])

    def reverse(self) -> "DTCode":
        """Orientation reverse via PD round-trip when possible."""
        if not self.code:
            return DTCode.empty()
        try:
            return DTCode.from_pd(self.to_pd().reverse())
        except Exception:
            return self.mirror()

    def normalised(self) -> "DTCode":
        """Lexicographically least among the code and its mirror."""
        if not self.code:
            return DTCode.empty()
        mirrored = tuple(-x for x in self.code)
        best = min(self.code, mirrored)
        return DTCode(best)

    def to_gauss(self) -> "GaussCode":
        """Convert DT code to a Gauss code."""
        from knottidy.gauss import GaussCode

        n = len(self.code)
        if n == 0:
            return GaussCode.empty()

        visits: dict[int, int] = {}  # step label 1..2n → signed crossing index 1..n
        for i, even in enumerate(self.code):
            odd = 2 * i + 1
            cidx = i + 1  # crossing number 1..n associated to this odd
            if even > 0:
                visits[odd] = -cidx
                visits[even] = cidx
            else:
                visits[odd] = cidx
                visits[abs(even)] = -cidx
        seq = [visits[s] for s in range(1, 2 * n + 1)]
        return GaussCode(seq)

    def to_pd(self) -> "PDCode":
        """Convert to PD via Gauss reconstruction."""
        return self.to_gauss().to_pd()

    def to_dict(self) -> dict:
        return {"code": list(self.code)}

    @classmethod
    def from_dict(cls, data: dict) -> DTCode:
        if not isinstance(data, dict) or "code" not in data:
            raise InvalidDTCodeError("from_dict expects dict with 'code'")
        return cls(data["code"])

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from knottidy.gauss import GaussCode
    from knottidy.pdcode import PDCode

__all__ = [
    "DTCode",
    "validate_dt_sequence",
    "parse_dt_text",
]
