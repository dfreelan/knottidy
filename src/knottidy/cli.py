"""Command-line interface for knottidy."""

from __future__ import annotations

import argparse
import sys
from typing import Optional

from knottidy import __version__
from knottidy.braid import Braid
from knottidy.dtcode import DTCode
from knottidy.exceptions import KnottidyError
from knottidy.gauss import GaussCode
from knottidy.io import dumps, dumps_braid, dumps_pd, parse_any, parse_braid_text, parse_pd_text
from knottidy.jones import jones_polynomial, jones_trefoil, jones_unknot
from knottidy.moves import reduce_braid_word, simplify_braid, simplify_pd
from knottidy.pdcode import PDCode
from knottidy.random_knots import random_braid, random_knot_braid

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="knottidy",
        description="Knot theory toolkit: braids, PD/DT/Gauss codes, Jones polynomial.",
    )
    p.add_argument("--version", action="version", version=f"knottidy {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    j = sub.add_parser("jones", help="compute Jones polynomial")
    g = j.add_mutually_exclusive_group(required=True)
    g.add_argument("--braid", type=str, help="braid word, e.g. '1 1 1'")
    g.add_argument("--pd", type=str, help="PD code text")
    g.add_argument("--trefoil", action="store_true", help="use the right-handed trefoil")
    g.add_argument("--unknot", action="store_true", help="use the unknot")
    g.add_argument("--figure-eight", action="store_true", help="use the figure-eight knot")
    j.add_argument("--strands", type=int, default=None, help="strand count for --braid")

    b = sub.add_parser("braid", help="braid word utilities")
    b_sub = b.add_subparsers(dest="braid_cmd", required=True)
    br = b_sub.add_parser("reduce", help="free-reduce a braid word")
    br.add_argument("word", type=str)
    br.add_argument("--strands", type=int, default=None)
    bi = b_sub.add_parser("info", help="show braid invariants")
    bi.add_argument("word", type=str)
    bi.add_argument("--strands", type=int, default=None)
    bc = b_sub.add_parser("close", help="braid closure as PD")
    bc.add_argument("word", type=str)
    bc.add_argument("--strands", type=int, default=None)

    c = sub.add_parser("convert", help="convert between representations")
    c.add_argument("--braid", type=str, default=None)
    c.add_argument("--pd", type=str, default=None)
    c.add_argument("--dt", type=str, default=None)
    c.add_argument("--gauss", type=str, default=None)
    c.add_argument("--strands", type=int, default=None)
    c.add_argument("--to", type=str, required=True, choices=["pd", "dt", "gauss", "braid", "json"])

    r = sub.add_parser("random", help="sample a random braid")
    r.add_argument("--strands", type=int, default=3)
    r.add_argument("--length", type=int, default=6)
    r.add_argument("--seed", type=int, default=None)
    r.add_argument("--knot", action="store_true", help="require knot closure")

    info = sub.add_parser("info", help="summarise a diagram")
    info.add_argument("--braid", type=str, default=None)
    info.add_argument("--pd", type=str, default=None)
    info.add_argument("--strands", type=int, default=None)

    sm = sub.add_parser("simplify", help="greedy combinatorial simplification")
    sm.add_argument("--braid", type=str, default=None)
    sm.add_argument("--pd", type=str, default=None)
    sm.add_argument("--strands", type=int, default=None)

    cat = sub.add_parser("catalogue", help="print named example objects")
    cat.add_argument("--kind", choices=["braid", "pd", "jones"], default="braid")

    return p

def _braid_from_args(word: str, strands: Optional[int]) -> Braid:
    if strands is None:
        return Braid.from_string(word)
    return Braid.from_string(word, strands=strands)

def _cmd_jones(args: argparse.Namespace) -> int:
    if args.unknot:
        print(jones_unknot().format())
        return 0
    if args.trefoil:
        poly = jones_polynomial(Braid.trefoil())
        print(poly.format())
        return 0
    if args.figure_eight:
        poly = jones_polynomial(Braid.figure_eight())
        print(poly.format())
        return 0
    if args.braid is not None:
        b = _braid_from_args(args.braid, args.strands)
        print(jones_polynomial(b).format())
        return 0
    if args.pd is not None:
        pd = parse_pd_text(args.pd)
        print(jones_polynomial(pd).format())
        return 0
    return 1

def _cmd_braid(args: argparse.Namespace) -> int:
    b = _braid_from_args(args.word, args.strands)
    if args.braid_cmd == "reduce":
        print(dumps_braid(reduce_braid_word(b)))
        return 0
    if args.braid_cmd == "info":
        print(f"strands: {b.strands}")
        print(f"word: {b.format()}")
        print(f"length: {len(b)}")
        print(f"writhe: {b.writhe()}")
        print(f"permutation: {b.permutation()}")
        print(f"components: {b.num_components()}")
        print(f"pure: {b.is_pure()}")
        return 0
    if args.braid_cmd == "close":
        print(dumps_pd(b.closure_pd()))
        return 0
    return 1

def _cmd_convert(args: argparse.Namespace) -> int:
    obj: object
    if args.braid is not None:
        obj = _braid_from_args(args.braid, args.strands)
    elif args.pd is not None:
        obj = parse_pd_text(args.pd)
    elif args.dt is not None:
        obj = DTCode.from_string(args.dt)
    elif args.gauss is not None:
        obj = GaussCode.from_string(args.gauss)
    else:
        print("provide one of --braid/--pd/--dt/--gauss", file=sys.stderr)
        return 2

    target = args.to
    if target == "json":
        from knottidy.io import to_json

        print(to_json(obj))
        return 0
    if target == "braid":
        if not isinstance(obj, Braid):
            print("conversion to braid is only available from braid input", file=sys.stderr)
            return 2
        print(dumps_braid(obj))
        return 0
    if target == "pd":
        if isinstance(obj, Braid):
            print(dumps_pd(obj.closure_pd()))
        elif isinstance(obj, PDCode):
            print(dumps_pd(obj))
        elif isinstance(obj, DTCode):
            print(dumps_pd(obj.to_pd()))
        elif isinstance(obj, GaussCode):
            print(dumps_pd(obj.to_pd()))
        else:
            return 2
        return 0
    if target == "dt":
        if isinstance(obj, Braid):
            pd = obj.closure_pd()
        elif isinstance(obj, PDCode):
            pd = obj
        elif isinstance(obj, DTCode):
            print(str(obj))
            return 0
        elif isinstance(obj, GaussCode):
            pd = obj.to_pd()
        else:
            return 2
        print(str(DTCode.from_pd(pd)))
        return 0
    if target == "gauss":
        if isinstance(obj, Braid):
            pd = obj.closure_pd()
        elif isinstance(obj, PDCode):
            pd = obj
        elif isinstance(obj, DTCode):
            print(dumps(obj.to_gauss()))
            return 0
        elif isinstance(obj, GaussCode):
            print(dumps(obj))
            return 0
        else:
            return 2
        print(dumps(GaussCode.from_pd(pd)))
        return 0
    return 2

def _cmd_random(args: argparse.Namespace) -> int:
    if args.knot:
        b = random_knot_braid(args.strands, args.length, seed=args.seed)
    else:
        b = random_braid(args.strands, args.length, seed=args.seed)
    print(dumps_braid(b))
    return 0

def _cmd_info(args: argparse.Namespace) -> int:
    if args.braid is not None:
        b = _braid_from_args(args.braid, args.strands)
        pd = b.closure_pd()
        print(f"braid: {b}")
        print(f"writhe(braid): {b.writhe()}")
        print(f"components: {b.num_components()}")
        print(f"pd_crossings: {len(pd)}")
        return 0
    if args.pd is not None:
        pd = parse_pd_text(args.pd)
        print(f"crossings: {len(pd)}")
        print(f"components: {pd.num_components()}")
        print(f"writhe: {pd.writhe()}")
        return 0
    print("provide --braid or --pd", file=sys.stderr)
    return 2

def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "jones":
            return _cmd_jones(args)
        if args.command == "braid":
            return _cmd_braid(args)
        if args.command == "convert":
            return _cmd_convert(args)
        if args.command == "random":
            return _cmd_random(args)
        if args.command == "info":
            return _cmd_info(args)
    except KnottidyError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1
    parser.print_help()
    return 2

if __name__ == "__main__":
    sys.exit(main())
