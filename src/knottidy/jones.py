"""Jones polynomial via the Kauffman bracket."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Optional, Union

from knottidy.exceptions import InvalidPDCodeError, PolynomialError

@dataclass(frozen=True)
class LaurentPolynomial:
    """Laurent polynomial with integer coefficients.

    Internally stored as a frozen dict ``{exponent: coefficient}`` with all
    coefficients nonzero. The indeterminate is formal (usually ``A`` or ``t``).
    """

    terms: dict[int, int]
    variable: str = "A"

    def __init__(
        self,
        terms: Optional[Union[dict[int, int], Iterable[tuple[int, int]]]] = None,
        variable: str = "A",
    ) -> None:
        if terms is None:
            cleaned: dict[int, int] = {}
        elif isinstance(terms, dict):
            cleaned = {int(e): int(c) for e, c in terms.items() if c != 0}
        else:
            cleaned = {}
            for e, c in terms:
                e, c = int(e), int(c)
                if c == 0:
                    continue
                cleaned[e] = cleaned.get(e, 0) + c
                if cleaned[e] == 0:
                    del cleaned[e]
        object.__setattr__(self, "terms", dict(sorted(cleaned.items())))
        object.__setattr__(self, "variable", variable)

    @classmethod
    def zero(cls, variable: str = "A") -> LaurentPolynomial:
        return cls({}, variable=variable)

    @classmethod
    def one(cls, variable: str = "A") -> LaurentPolynomial:
        return cls({0: 1}, variable=variable)

    @classmethod
    def monomial(cls, coeff: int, exp: int, variable: str = "A") -> LaurentPolynomial:
        if coeff == 0:
            return cls.zero(variable)
        return cls({exp: coeff}, variable=variable)

    @classmethod
    def from_int(cls, n: int, variable: str = "A") -> LaurentPolynomial:
        return cls({0: n}, variable=variable) if n != 0 else cls.zero(variable)

    def is_zero(self) -> bool:
        return not self.terms

    def degree(self) -> Optional[int]:
        return max(self.terms) if self.terms else None

    def min_degree(self) -> Optional[int]:
        return min(self.terms) if self.terms else None

    def coefficient(self, exp: int) -> int:
        return self.terms.get(exp, 0)

    def support(self) -> list[int]:
        return list(self.terms.keys())

    def __bool__(self) -> bool:
        return bool(self.terms)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, int):
            if other == 0:
                return self.is_zero()
            return self.terms == {0: other}
        if not isinstance(other, LaurentPolynomial):
            return NotImplemented
        return self.terms == other.terms

    def __hash__(self) -> int:
        return hash((tuple(self.terms.items()), self.variable))

    def _alike(self, other: LaurentPolynomial) -> None:
        if self.variable != other.variable:
            raise PolynomialError(
                f"variable mismatch: {self.variable!r} vs {other.variable!r}"
            )

    def __add__(self, other: object) -> LaurentPolynomial:
        if isinstance(other, int):
            other = LaurentPolynomial.from_int(other, self.variable)
        if not isinstance(other, LaurentPolynomial):
            return NotImplemented
        self._alike(other)
        acc: dict[int, int] = dict(self.terms)
        for e, c in other.terms.items():
            acc[e] = acc.get(e, 0) + c
            if acc[e] == 0:
                del acc[e]
        return LaurentPolynomial(acc, self.variable)

    def __radd__(self, other: object) -> LaurentPolynomial:
        return self.__add__(other)

    def __neg__(self) -> LaurentPolynomial:
        return LaurentPolynomial({e: -c for e, c in self.terms.items()}, self.variable)

    def __sub__(self, other: object) -> LaurentPolynomial:
        return self + (-other)  # type: ignore[operator]

    def __rsub__(self, other: object) -> LaurentPolynomial:
        return (-self) + other

    def __mul__(self, other: object) -> LaurentPolynomial:
        if isinstance(other, int):
            if other == 0:
                return LaurentPolynomial.zero(self.variable)
            return LaurentPolynomial(
                {e: c * other for e, c in self.terms.items()}, self.variable
            )
        if not isinstance(other, LaurentPolynomial):
            return NotImplemented
        self._alike(other)
        acc: dict[int, int] = {}
        for e1, c1 in self.terms.items():
            for e2, c2 in other.terms.items():
                e = e1 + e2
                acc[e] = acc.get(e, 0) + c1 * c2
                if acc[e] == 0:
                    del acc[e]
        return LaurentPolynomial(acc, self.variable)

    def __rmul__(self, other: object) -> LaurentPolynomial:
        return self.__mul__(other)

    def __pow__(self, n: int) -> LaurentPolynomial:
        if not isinstance(n, int) or isinstance(n, bool):
            return NotImplemented
        if n < 0:
            raise PolynomialError("negative powers require units; use shift/monomial")
        result = LaurentPolynomial.one(self.variable)
        base = self
        while n:
            if n & 1:
                result = result * base
            base = base * base
            n >>= 1
        return result

    def shift(self, k: int) -> LaurentPolynomial:
        """Multiply by ``variable^k``."""
        return LaurentPolynomial({e + k: c for e, c in self.terms.items()}, self.variable)

    def rename(self, variable: str) -> LaurentPolynomial:
        return LaurentPolynomial(self.terms, variable=variable)

    def __repr__(self) -> str:
        return f"LaurentPolynomial({self.terms!r}, variable={self.variable!r})"

    def __str__(self) -> str:
        return self.format()

    def format(self, variable: Optional[str] = None) -> str:
        v = variable or self.variable
        if not self.terms:
            return "0"
        parts: list[str] = []
        for e, c in sorted(self.terms.items(), reverse=True):
            if e == 0:
                term = str(c)
            elif e == 1:
                if c == 1:
                    term = v
                elif c == -1:
                    term = f"-{v}"
                else:
                    term = f"{c}*{v}"
            else:
                if c == 1:
                    term = f"{v}^{e}"
                elif c == -1:
                    term = f"-{v}^{e}"
                else:
                    term = f"{c}*{v}^{e}"
            parts.append(term)
        out = parts[0]
        for p in parts[1:]:
            if p.startswith("-"):
                out += " - " + p[1:]
            else:
                out += " + " + p
        return out

    def to_dict(self) -> dict:
        return {"terms": dict(self.terms), "variable": self.variable}

    @classmethod
    def from_dict(cls, data: dict) -> LaurentPolynomial:
        return cls(data["terms"], variable=data.get("variable", "A"))

def _delta(A_var: str = "A") -> LaurentPolynomial:
    """δ = -A² - A^{-2}."""
    return LaurentPolynomial({2: -1, -2: -1}, variable=A_var)

def kauffman_bracket(pd: object, *, max_crossings: int = 16) -> LaurentPolynomial:
    """Compute the Kauffman bracket of a PD diagram by state sum."""
    from knottidy.pdcode import PDCode

    if not isinstance(pd, PDCode):
        raise PolynomialError("kauffman_bracket expects a PDCode")
    n = len(pd.crossings)
    if n > max_crossings:
        raise PolynomialError(
            f"refusing state-sum on {n} crossings (max_crossings={max_crossings})"
        )
    if n == 0:
        return LaurentPolynomial.one("A")

    crossings = pd.crossings
    # Handedness is already encoded in PD cyclic order for braid closures;
    # the state sum uses the stored (a,b,c,d) directly.
    labels = list(range(1, 2 * n + 1))

    def count_circles(state: int) -> int:
        parent = {e: e for e in labels}

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(x: int, y: int) -> None:
            rx, ry = find(x), find(y)
            if rx != ry:
                parent[ry] = rx

        for i, cr in enumerate(crossings):
            a, b, c, d = cr.a, cr.b, cr.c, cr.d
            if (state >> i) & 1 == 0:
                union(a, d)
                union(b, c)
            else:
                union(a, b)
                union(c, d)
        roots = {find(e) for e in labels}
        return len(roots)

    total = LaurentPolynomial.zero("A")
    for state in range(1 << n):
        a_count = sum(1 for i in range(n) if (state >> i) & 1 == 0)
        b_count = n - a_count
        circles = count_circles(state)
        exp_A = a_count - b_count
        contrib = LaurentPolynomial.monomial(1, exp_A, "A")
        if circles - 1 > 0:
            contrib = contrib * (_delta("A") ** (circles - 1))
        elif circles - 1 < 0:
            raise PolynomialError(f"non-positive circle count {circles}")
        total = total + contrib
    return total

def jones_from_bracket(
    bracket: LaurentPolynomial,
    writhe: int,
) -> LaurentPolynomial:
    """Normalise a Kauffman bracket to the Jones polynomial in ``t``."""
    sign = 1 if (writhe % 2 == 0) else -1
    normalised = bracket * sign
    normalised = normalised.shift(-3 * writhe)

    acc: dict[int, int] = {}
    for e, c in normalised.terms.items():
        if e % 4 == 0:
            t_exp = -e // 4
            acc[t_exp] = acc.get(t_exp, 0) + c
            if acc[t_exp] == 0:
                del acc[t_exp]
        elif e % 2 == 0:
            raise PolynomialError(
                f"bracket exponent {e} not divisible by 4; "
                "half-integer Jones (links) — use jones_polynomial_half"
            )
        else:
            raise PolynomialError(f"unexpected A-exponent {e} not divisible by 2")
    return LaurentPolynomial(acc, variable="t")

def jones_polynomial_half(
    bracket: LaurentPolynomial,
    writhe: int,
) -> LaurentPolynomial:
    """Jones polynomial as a Laurent polynomial in ``√t`` (variable ``q`` with ``q²=t``)."""
    sign = 1 if (writhe % 2 == 0) else -1
    normalised = (bracket * sign).shift(-3 * writhe)
    acc: dict[int, int] = {}
    for e, c in normalised.terms.items():
        if e % 2 != 0:
            raise PolynomialError(f"unexpected A-exponent {e}")
        if e % 2 == 0:
            q_exp = -e // 2
            acc[q_exp] = acc.get(q_exp, 0) + c
            if acc[q_exp] == 0:
                del acc[q_exp]
    return LaurentPolynomial(acc, variable="q")

def jones_polynomial(diagram: object, *, max_crossings: int = 16) -> LaurentPolynomial:
    """Compute the Jones polynomial of a knot/link diagram."""
    from knottidy.braid import Braid
    from knottidy.pdcode import PDCode

    if isinstance(diagram, Braid):
        pd = diagram.closure_pd()
    elif isinstance(diagram, PDCode):
        pd = diagram
    else:
        raise PolynomialError(f"unsupported diagram type {type(diagram)!r}")

    bracket = kauffman_bracket(pd, max_crossings=max_crossings)
    w = pd.writhe()
    if isinstance(diagram, Braid):
        w = diagram.writhe()
        pd = pd._with_writhe(w)

    try:
        return jones_from_bracket(bracket, w)
    except PolynomialError:
        return jones_polynomial_half(bracket, w)

def jones_unknot() -> LaurentPolynomial:
    """Known value: ``V(unknot) = 1``."""
    return LaurentPolynomial.one("t")

def jones_trefoil() -> LaurentPolynomial:
    """Known value: right-handed trefoil ``V = t + t^3 - t^4``."""
    return LaurentPolynomial({1: 1, 3: 1, 4: -1}, variable="t")

def jones_left_trefoil() -> LaurentPolynomial:
    """Left-handed trefoil ``V = t^{-1} + t^{-3} - t^{-4}``."""
    return LaurentPolynomial({-1: 1, -3: 1, -4: -1}, variable="t")

def jones_figure_eight() -> LaurentPolynomial:
    """Figure-eight ``V = t^{-2} - t^{-1} + 1 - t + t^2``."""
    return LaurentPolynomial({-2: 1, -1: -1, 0: 1, 1: -1, 2: 1}, variable="t")

def jones_hopf_link() -> LaurentPolynomial:
    """Positive Hopf link ``V = -t^{-1/2} - t^{1/2}`` as polynomial in ``q=√t``:
    ``-q^{-1} - q``.
    """
    return LaurentPolynomial({-1: -1, 1: -1}, variable="q")

def is_unknot_jones(poly: LaurentPolynomial) -> bool:
    """Return True if ``poly`` equals the unknot Jones polynomial."""
    return poly == jones_unknot() or (
        poly.variable in ("t", "A") and poly.terms == {0: 1}
    )

def verify_trefoil_jones(poly: LaurentPolynomial) -> bool:
    """Return True if ``poly`` matches right- or left-handed trefoil Jones."""
    return poly == jones_trefoil() or poly == jones_left_trefoil()

__all__ = [
    "LaurentPolynomial",
    "kauffman_bracket",
    "jones_from_bracket",
    "jones_polynomial",
    "jones_polynomial_half",
    "jones_unknot",
    "jones_trefoil",
    "jones_left_trefoil",
    "jones_figure_eight",
    "jones_hopf_link",
    "is_unknot_jones",
    "verify_trefoil_jones",
]
