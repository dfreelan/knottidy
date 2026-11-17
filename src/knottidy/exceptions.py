"""Exception hierarchy for knottidy.

All library-specific errors derive from :class:`KnottidyError` so callers can
catch a single base type when desired.
"""

from __future__ import annotations

class KnottidyError(Exception):
    """Base class for all knottidy errors."""

class InvalidBraidError(KnottidyError):
    """Raised when a braid word or strand count is inconsistent.

    Examples include generators outside the legal range for the stated number
    of strands, an empty strand count, or malformed word tokens.
    """

class InvalidPDCodeError(KnottidyError):
    """Raised when a planar diagram (PD) code fails structural checks.

    Typical causes: wrong arity of a crossing tuple, duplicated edge labels,
    gaps in the edge labelling, or an inconsistent over/under pairing.
    """

class InvalidDTCodeError(KnottidyError):
    """Raised when a Dowker–Thistlethwaite code is malformed.

    DT codes must consist of even integers with absolute values forming a
    permutation of ``{2, 4, ..., 2n}`` for an ``n``-crossing diagram.
    """

class InvalidGaussCodeError(KnottidyError):
    """Raised when a Gauss code fails combinatorial validity checks."""

class ConversionError(KnottidyError):
    """Raised when a representation cannot be converted to another format.

    Conversion may fail for multi-component links when a knot-only format is
    requested, or when an intermediate invariant is undefined.
    """

class PolynomialError(KnottidyError):
    """Raised for malformed Laurent polynomials or illegal substitutions."""

class MoveError(KnottidyError):
    """Raised when a requested Reidemeister / Markov move cannot be applied.

    The diagram may not admit the move at the given location, or the move
    parameters may be out of range.
    """

class IOFormatError(KnottidyError):
    """Raised when a serialised knot/link description cannot be parsed."""

class RandomGenerationError(KnottidyError):
    """Raised when random generation parameters are invalid."""

__all__ = [
    "KnottidyError",
    "InvalidBraidError",
    "InvalidPDCodeError",
    "InvalidDTCodeError",
    "InvalidGaussCodeError",
    "ConversionError",
    "PolynomialError",
    "MoveError",
    "IOFormatError",
    "RandomGenerationError",
]
