"""Import / export helpers for knot and braid descriptions."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any, Optional, Union

from knottidy.braid import Braid, format_braid_word, parse_braid_word
from knottidy.dtcode import DTCode, parse_dt_text
from knottidy.exceptions import IOFormatError
from knottidy.gauss import GaussCode, parse_gauss_text
from knottidy.jones import LaurentPolynomial
from knottidy.pdcode import PDCode

PathLike = Union[str, Path]

def detect_format(text: str) -> str:
    """Heuristic format detection for a text payload."""
    s = text.strip()
    if not s:
        raise IOFormatError("empty input")
    if s[0] == "{":
        return "json"
    low = s.lower()
    if low.startswith("pd") or (s.startswith("[[") or s.startswith("[X")):
        return "pd"
    if low.startswith("dt"):
        return "dt"
    if low.startswith("gauss") or low.startswith("gc"):
        return "gauss"
    if low.startswith("b[") or low.startswith("braid") or low.startswith("b_"):
        return "braid"
    if re.fullmatch(r"[\d\s,\-\|\[\]OUou;]+", s):
        try:
            nums = [int(x) for x in re.findall(r"-?\d+", s)]
            if nums and all(abs(x) % 2 == 0 for x in nums) and len(nums) == max(abs(x) for x in nums) // 2:
                return "dt"
            if nums:
                return "gauss"
        except ValueError:
            pass
    if re.fullmatch(r"[\d\s,\-σsS_^\(\)]+", s):
        return "braid"
    raise IOFormatError(f"cannot detect format for text starting {s[:40]!r}")

def parse_pd_text(text: str) -> PDCode:
    """Parse ``PD[[a,b,c,d], ...]`` or ``[[a,b,c,d], ...]``."""
    s = text.strip()
    s = re.sub(r"^[Pp][Dd]\s*", "", s)
    s = s.strip()
    if not s or s in ("[]", "PD[]"):
        return PDCode.empty()
    s = re.sub(r"[Xx]\s*\(", "[", s)
    s = s.replace(")", "]")
    try:
        data = json.loads(s.replace("(", "[").replace(")", "]"))
    except json.JSONDecodeError as exc:
        raise IOFormatError(f"cannot parse PD text: {exc}") from exc
    if not isinstance(data, list):
        raise IOFormatError("PD text must be a list of crossings")
    return PDCode(data)

def parse_any(text: str, fmt: Optional[str] = None) -> Any:
    """Parse text into a Braid, PDCode, DTCode, or GaussCode."""
    if fmt is None:
        fmt = detect_format(text)
    fmt = fmt.lower()
    if fmt == "json":
        return from_json(text)
    if fmt == "braid":
        return parse_braid_text(text)
    if fmt == "pd":
        return parse_pd_text(text)
    if fmt == "dt":
        return DTCode.from_string(text)
    if fmt == "gauss":
        return GaussCode.from_string(_strip_prefix(text, ("gauss", "gc")))
    raise IOFormatError(f"unknown format {fmt!r}")

def _strip_prefix(text: str, prefixes: tuple[str, ...]) -> str:
    s = text.strip()
    low = s.lower()
    for p in prefixes:
        if low.startswith(p):
            s = s[len(p) :].strip()
            break
    return s.lstrip(":").strip()

def parse_braid_text(text: str) -> Braid:
    """Parse braid text including ``B[n](...)`` and ``braid n: ...`` forms."""
    s = text.strip()
    m = re.match(r"^[Bb]\[(\d+)\]\s*\((.*)\)$", s, re.DOTALL)
    if m:
        return Braid.from_string(m.group(2), strands=int(m.group(1)))
    m = re.match(r"^[Bb]_(\d+)\s*[<:]\s*(.*)$", s)
    if m:
        return Braid.from_string(m.group(2), strands=int(m.group(1)))
    m = re.match(r"^[Bb]raid\s+(\d+)\s*:\s*(.*)$", s, re.IGNORECASE)
    if m:
        return Braid.from_string(m.group(2), strands=int(m.group(1)))
    word, n = parse_braid_word(s)
    return Braid(word, n)

def dumps_braid(braid: Braid, style: str = "bracket") -> str:
    if style == "bracket":
        return f"B[{braid.strands}]({format_braid_word(braid.word)})"
    if style == "signed":
        return format_braid_word(braid.word, "signed")
    if style == "sigma":
        return format_braid_word(braid.word, "sigma")
    raise IOFormatError(f"unknown braid dump style {style!r}")

def dumps_pd(pd: PDCode) -> str:
    if not pd.crossings:
        return "PD[]"
    inner = ",".join("[" + ",".join(str(x) for x in c.as_tuple()) + "]" for c in pd.crossings)
    return f"PD[{inner}]"

def dumps_dt(dt: DTCode) -> str:
    return str(dt)

def dumps_gauss(gauss: GaussCode) -> str:
    return "Gauss[" + gauss.format() + "]"

def dumps(obj: Any, fmt: Optional[str] = None) -> str:
    """Serialise a knottidy object to text."""
    if fmt == "json" or (fmt is None and False):
        return to_json(obj)
    if isinstance(obj, Braid):
        return dumps_braid(obj)
    if isinstance(obj, PDCode):
        return dumps_pd(obj)
    if isinstance(obj, DTCode):
        return dumps_dt(obj)
    if isinstance(obj, GaussCode):
        return dumps_gauss(obj)
    if isinstance(obj, LaurentPolynomial):
        return obj.format()
    if fmt == "json":
        return to_json(obj)
    raise IOFormatError(f"cannot dump object of type {type(obj)!r}")

def to_json(obj: Any, *, indent: Optional[int] = 2) -> str:
    """Serialise to a JSON string with a ``type`` discriminator."""
    return json.dumps(to_jsonable(obj), indent=indent)

def to_jsonable(obj: Any) -> dict:
    if isinstance(obj, Braid):
        return {"type": "braid", **obj.to_dict()}
    if isinstance(obj, PDCode):
        return {"type": "pd", **obj.to_dict()}
    if isinstance(obj, DTCode):
        return {"type": "dt", **obj.to_dict()}
    if isinstance(obj, GaussCode):
        return {"type": "gauss", **obj.to_dict()}
    if isinstance(obj, LaurentPolynomial):
        return {"type": "laurent", **obj.to_dict()}
    raise IOFormatError(f"cannot JSON-encode {type(obj)!r}")

def from_json(text: str) -> Any:
    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        raise IOFormatError(f"invalid JSON: {exc}") from exc
    return from_jsonable(data)

def from_jsonable(data: dict) -> Any:
    if not isinstance(data, dict) or "type" not in data:
        raise IOFormatError("JSON object must have a 'type' field")
    t = data["type"]
    if t == "braid":
        return Braid.from_dict(data)
    if t == "pd":
        return PDCode.from_dict(data)
    if t == "dt":
        return DTCode.from_dict(data)
    if t == "gauss":
        return GaussCode.from_dict(data)
    if t == "laurent":
        return LaurentPolynomial.from_dict(data)
    raise IOFormatError(f"unknown JSON type {t!r}")

def load(path: PathLike, fmt: Optional[str] = None) -> Any:
    """Load a knot/braid description from a file."""
    p = Path(path)
    try:
        text = p.read_text(encoding="utf-8")
    except OSError as exc:
        raise IOFormatError(f"cannot read {p}: {exc}") from exc
    if fmt is None:
        suffix = p.suffix.lower()
        if suffix == ".json":
            fmt = "json"
        elif suffix in (".pd",):
            fmt = "pd"
        elif suffix in (".dt",):
            fmt = "dt"
        elif suffix in (".gauss", ".gc"):
            fmt = "gauss"
        elif suffix in (".braid", ".br"):
            fmt = "braid"
    return parse_any(text, fmt=fmt)

def dump(obj: Any, path: PathLike, fmt: Optional[str] = None) -> None:
    """Write a knot/braid description to a file."""
    p = Path(path)
    if fmt is None:
        suffix = p.suffix.lower()
        if suffix == ".json":
            fmt = "json"
    if fmt == "json":
        text = to_json(obj)
    else:
        text = dumps(obj, fmt=fmt)
    try:
        p.write_text(text + "\n", encoding="utf-8")
    except OSError as exc:
        raise IOFormatError(f"cannot write {p}: {exc}") from exc

__all__ = [
    "detect_format",
    "parse_any",
    "parse_pd_text",
    "parse_braid_text",
    "dumps",
    "dumps_braid",
    "dumps_pd",
    "dumps_dt",
    "dumps_gauss",
    "to_json",
    "to_jsonable",
    "from_json",
    "from_jsonable",
    "load",
    "dump",
]
