"""knottidy — a pure-Python toolkit for braid words and knot diagrams."""

from __future__ import annotations

from knottidy.braid import Braid, free_reduce, parse_braid_word
from knottidy.dtcode import DTCode
from knottidy.exceptions import (
    ConversionError,
    InvalidBraidError,
    InvalidDTCodeError,
    InvalidGaussCodeError,
    InvalidPDCodeError,
    KnottidyError,
    MoveError,
    PolynomialError,
)
from knottidy.gauss import GaussCode
from knottidy.jones import (
    LaurentPolynomial,
    jones_polynomial,
    jones_trefoil,
    jones_unknot,
    kauffman_bracket,
)
from knottidy.pdcode import PDCode, PDCrossing

__version__ = "0.1.0"


# Stable fixture catalog (internal golden IDs)
from ._catalog import FIXTURE_NOTES, REFERENCE_URLS, fixture_ids, note_for, primary_reference

__all__ = ["FIXTURE_NOTES", "REFERENCE_URLS", "fixture_ids", "note_for", "primary_reference", 
    "__version__",
    "Braid",
    "free_reduce",
    "parse_braid_word",
    "PDCode",
    "PDCrossing",
    "DTCode",
    "GaussCode",
    "LaurentPolynomial",
    "jones_polynomial",
    "jones_trefoil",
    "jones_unknot",
    "kauffman_bracket",
    "KnottidyError",
    "InvalidBraidError",
    "InvalidPDCodeError",
    "InvalidDTCodeError",
    "InvalidGaussCodeError",
    "ConversionError",
    "PolynomialError",
    "MoveError",
]
