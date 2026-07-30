#!/usr/bin/env python3
"""check_numbers.py - is every number printed in Part IV greppable in an archived output?

  Author: Carles Marin <karlesmarin@gmail.com>  (with Claude, Anthropic, as assistant)

The rule adopted for the companion paper and applied here: a number that appears in the text must
be locatable in anc_iv/outputs/, or be flagged.  Anything else is hearsay.  Numbers that are part
of a formula rather than a measurement (indices, years, section numbers) are listed separately
rather than silently dropped.
"""
import os
import re

TEX = "ghu_secondfixed.tex"
OUT = "anc_iv/outputs"

corpus = ""
for fn in sorted(os.listdir(OUT)):
    corpus += open(os.path.join(OUT, fn), encoding="utf-8", errors="ignore").read()

s = open(TEX, encoding="utf-8").read()
s = re.sub(r"%.*", "", s)                       # drop comments
s = re.sub(r"\\begin\{thebibliography\}.*", "", s, flags=re.S)   # drop the bibliography

# numbers written as $3060$ or 3060 or 10{,}920, at least three digits: those are measurements
raw = re.findall(r"(?<![\w\\^_{])(\d[\d]{2,}(?:\{,\}\d{3})*)(?![\w}])", s)
seen, rows = set(), []
for n in raw:
    plain = n.replace("{,}", "")
    if plain in seen:
        continue
    seen.add(plain)
    with_comma = format(int(plain), ",")
    found = plain in corpus or with_comma in corpus
    rows.append((plain, found))

print("numbers of three digits or more in %s, checked against %s/*\n" % (TEX, OUT))
missing = []
for plain, found in sorted(rows, key=lambda r: -int(r[0])):
    print("  %-8s %s" % (plain, "greppable" if found else "NOT FOUND  <-- check"))
    if not found:
        missing.append(plain)
print("\n%d distinct numbers, %d not backed by an archived run" % (len(rows), len(missing)))
