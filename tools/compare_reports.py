#!/usr/bin/env python3
"""Compare two diffReps reports (e.g. diffReps.pl vs diffreps-rs).

Usage: compare_reports.py REFERENCE OTHER [--ptol 1e-13]

Prints:
  * header values that differ (timings, versions, paths and the command line
    are ignored);
  * regions present in only one report (matched by chrom/start/end);
  * for shared regions, how many differ per column. Columns pval, padj, winP
    and winQ are compared with a relative tolerance (--ptol), all others as
    exact strings.

Exit status 0 if the reports agree (within --ptol), 1 otherwise.

Note: unless the Rust run used --bg-from-report REFERENCE, the background
(mean/deviation) comes from a different random sample than diffReps', so a
few borderline windows can legitimately differ.
"""
import argparse
import sys

P_COLS = {"pval", "padj", "winP", "winQ"}
IGNORED_HEADER = ("seconds", "version", "Full command", "Output file", "files", "processors")


def load(path):
    header, rows, cols = {}, {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                if ": " in line:
                    key, value = line.lstrip("# ").split(": ", 1)
                    header[key] = value
                continue
            fields = line.split("\t")
            if cols is None:
                cols = fields
                continue
            rows[(fields[0], int(fields[1]), int(fields[2]))] = fields
    return header, cols, rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("reference")
    ap.add_argument("other")
    ap.add_argument("--ptol", type=float, default=1e-13, help="relative tolerance for p-value columns")
    args = ap.parse_args()

    h1, cols, r1 = load(args.reference)
    h2, cols2, r2 = load(args.other)
    if cols != cols2:
        sys.exit("table columns differ:\n  %s\n  %s" % (cols, cols2))
    ok = True

    print("== Header")
    for key in h1:
        if any(word in key for word in IGNORED_HEADER):
            continue
        if h1.get(key) != h2.get(key):
            ok = False
            print("  %s:\n    reference: %s\n    other:     %s" % (key, h1.get(key), h2.get(key)))

    only1 = sorted(set(r1) - set(r2))
    only2 = sorted(set(r2) - set(r1))
    shared = sorted(set(r1) & set(r2))
    print("== Regions: reference %d, other %d, shared %d, only in reference %d, only in other %d"
          % (len(r1), len(r2), len(shared), len(only1), len(only2)))
    for label, keys in (("only in reference", only1), ("only in other", only2)):
        for k in keys[:5]:
            print("  %s: %s:%d-%d" % (label, *k))
        if len(keys) > 5:
            print("  ... and %d more" % (len(keys) - 5))
    ok &= not only1 and not only2

    diff_counts = {c: 0 for c in cols}
    worst = {c: 0.0 for c in P_COLS}
    examples = []
    for k in shared:
        a, b = r1[k], r2[k]
        for c, x, y in zip(cols, a, b):
            if x == y:
                continue
            if c in P_COLS:
                rel = abs(float(x) - float(y)) / max(abs(float(x)), 1e-300)
                worst[c] = max(worst[c], rel)
                if rel <= args.ptol:
                    continue
            diff_counts[c] += 1
            if len(examples) < 5:
                examples.append("%s:%d-%d %s: %s vs %s" % (*k, c, x, y))
    print("== Shared regions, columns that differ (p-values: relative difference > %g)" % args.ptol)
    any_diff = False
    for c in cols:
        if diff_counts[c]:
            any_diff = True
            print("  %-14s %d" % (c, diff_counts[c]))
    if not any_diff:
        print("  none")
    print("  worst relative p-value differences: " + ", ".join("%s %.1e" % (c, worst[c]) for c in sorted(worst)))
    for e in examples:
        print("  e.g. " + e)
    ok &= not any_diff

    print("== RESULT: " + ("reports agree" if ok else "reports differ"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
