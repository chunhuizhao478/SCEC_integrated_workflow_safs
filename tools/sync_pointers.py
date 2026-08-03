#!/usr/bin/env python3
"""sync_pointers.py -- keep the `deckbuild/x.py:NNN` pointers in the notebooks correct.

The STEP comments cite exact line numbers so a reader can jump straight to the function
they want to change. Those numbers rot the instant anything above them moves, and a
pointer that has drifted onto a DIFFERENT function is worse than no pointer at all -- it
sends someone to edit the wrong code with confidence.

`tests/test_notebook.py` guards them. This rewrites them, so the guard is a two-second fix
instead of a manual hunt:

    python tools/sync_pointers.py            # rewrite the numbers in place
    python tools/sync_pointers.py --dry-run  # just report what would change

It only ever moves a number to the definition of the symbol NAMED ON THE SAME LINE, so it
cannot invent a pointer or silently repoint one at something else. A line whose symbol is
not found is reported and left alone.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NOTEBOOKS = ("deck_workflow.ipynb", "safs_check.ipynb")
PTR = re.compile(r"(deckbuild/\w+\.py):(\d+)")


def _defs(path: Path) -> dict[str, int]:
    """{symbol: 1-based line} for every def/class/module-level binding."""
    out: dict[str, int] = {}
    for i, line in enumerate(path.read_text().split("\n"), start=1):
        m = re.match(r"\s*(?:def|class)\s+(\w+)", line)
        if m is None:
            m = re.match(r"([A-Z_][A-Z0-9_]*)\s*(?::[^=]+)?=", line)
        if m and m.group(1) not in out:
            out[m.group(1)] = i
    return out


def _names_on(line: str, before: str) -> list[str]:
    head = line.split(before)[0]
    return (re.findall(r"\b(_?[A-Za-z]\w*)\s*\(\)", head)
            + re.findall(r"^\s*#\s+(_?[A-Z][A-Z0-9_]*)\s{2,}", line))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)

    cache: dict[str, dict[str, int]] = {}
    changed = fixed = unresolved = 0
    for nb_name in NOTEBOOKS:
        nb_path = ROOT / nb_name
        nb = json.loads(nb_path.read_text())
        dirty = False
        for ci, cell in enumerate(nb["cells"]):
            if cell["cell_type"] != "code":
                continue
            src = "".join(cell["source"])
            out_lines = []
            for line in src.split("\n"):
                m = PTR.search(line)
                if m:
                    rel, num = m.group(1), int(m.group(2))
                    if rel not in cache:
                        f = ROOT / rel
                        cache[rel] = _defs(f) if f.is_file() else {}
                    defs = cache[rel]
                    hit = next((n for n in _names_on(line, rel) if n in defs), None)
                    if hit is None:
                        if _names_on(line, rel):
                            print(f"  ? {nb_name} cell {ci}: cannot resolve "
                                  f"{_names_on(line, rel)[0]!r} in {rel}", file=sys.stderr)
                            unresolved += 1
                    elif defs[hit] != num:
                        line = line.replace(f"{rel}:{num}", f"{rel}:{defs[hit]}")
                        print(f"  {rel}:{num} -> {defs[hit]}  ({hit})")
                        fixed += 1
                        dirty = True
                out_lines.append(line)
            if dirty:
                s = "\n".join(out_lines).rstrip("\n")
                parts = s.split("\n")
                cell["source"] = [x + "\n" for x in parts[:-1]] + [parts[-1]]
        if dirty and not a.dry_run:
            nb_path.write_text(json.dumps(nb, indent=1))
            changed += 1
    verb = "would fix" if a.dry_run else "fixed"
    print(f"{verb} {fixed} pointer(s) across {changed} notebook(s); "
          f"{unresolved} unresolved")
    return 1 if unresolved else 0


if __name__ == "__main__":
    raise SystemExit(main())
