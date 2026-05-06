#!/usr/bin/env python3
"""
expand_network.py — postprocessor for SupNet's HC pipeline.

Reads a network produced by SupNet on a *compressed* gene-tree set, plus
the JSON map written by compress_gtrees.py, and replaces each synthetic
leaf with the original subtree it represents (or a sub-HC mini-network
inferred from the original GT-fragments — see --sub-hc).

The output is a network newick over the original species. Reticulation
markers (#1, #A, etc.) in the input network are preserved literally —
this tool operates at the string level on synthetic-label tokens, so it
doesn't need to parse the full network grammar. Sub-HC mini-network
reticulations are renumbered to avoid colliding with main-network IDs.

Run with `--help` for options.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile


SYNTHETIC_LABEL_RE = re.compile(r"(?<![A-Za-z0-9_])(_cmp_[0-9a-f]+)(?![A-Za-z0-9_])")
RETIC_RE           = re.compile(r"#([A-Za-z]?\d+)")


def collect_main_retic_ids(network_newick: str) -> set:
    """Extract reticulation labels (the part after `#`) used by the main
    network. Sub-HC mini-networks must use disjoint IDs after renumeration."""
    return set(RETIC_RE.findall(network_newick))


def renumber_subnet_retics(subnet_newick: str, used_ids: set, prefix: str) -> tuple:
    """Renumber every reticulation label in a sub-network so it doesn't
    collide with `used_ids`. Returns (new_newick, new_ids_added). Strategy:
    each old retic label is replaced by `<prefix><n>` where n grows; if any
    such name still collides, n is incremented until it's free."""
    old_ids = set(RETIC_RE.findall(subnet_newick))
    mapping = {}
    next_n = 1
    for old in sorted(old_ids):
        while True:
            new = f"{prefix}{next_n}"
            next_n += 1
            if new not in used_ids and new not in mapping.values():
                mapping[old] = new
                break
    new_newick = subnet_newick
    # Build replace pairs sorted longest-first so e.g. `#10` doesn't get
    # touched by `#1` (descending length avoids prefix overlap).
    for old, new in sorted(mapping.items(), key=lambda x: -len(x[0])):
        new_newick = re.sub(r"#" + re.escape(old) + r"(?![A-Za-z0-9_])",
                            f"#{new}", new_newick)
    return new_newick, set(mapping.values())


def run_sub_hc(fragments: list, supnet_path: str, sub_R: int,
               sub_cost: str, sub_climb: int, randseed: int = 7,
               verbose: bool = False) -> str:
    """Run SupNet HC on a list of GT-fragment newick strings. Returns the
    inferred mini-network newick (without trailing `;`), or None on failure.

    The fragments are written to a temp file, supnet is invoked with
    `--HC --noodtfiles -O <prefix>`, and we read the first network from
    `<prefix>.log`."""
    with tempfile.TemporaryDirectory() as td:
        gts_path = os.path.join(td, "frag.txt")
        out_prefix = os.path.join(td, "sub")
        with open(gts_path, "w") as f:
            for frag in fragments:
                f.write(frag.strip() + "\n")
        cmd = [
            supnet_path, "-G", gts_path, "-q1", f"-R{sub_R}",
            f"-C{sub_cost}", "--HC", f"--hcstopclimb={sub_climb}",
            "-O", out_prefix, f"--randseed={randseed}", "-v0",
        ]
        if verbose:
            print(f"[sub-hc] {' '.join(cmd)}", file=sys.stderr)
        try:
            subprocess.run(cmd, check=True, capture_output=True, timeout=600)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
            if verbose:
                print(f"[sub-hc] FAILED: {e}", file=sys.stderr)
            return None
        log_path = out_prefix + ".log"
        if not os.path.exists(log_path):
            if verbose:
                print(f"[sub-hc] no .log produced", file=sys.stderr)
            return None
        with open(log_path) as f:
            line = f.readline().strip()
        return line.rstrip(";") or None


def expand_one_network(network_newick: str, label_to_subtree: dict) -> str:
    """Replace every synthetic label in `network_newick` with the
    corresponding subtree text. The replacement is purely textual; we use
    a word-boundary regex so partial matches inside other identifiers are
    not affected. Subtree newicks are stripped of trailing `;`.
    Repeated until no more synthetic labels remain (handles the case
    where representatives themselves contain synthetic labels — should not
    happen, but defensive)."""
    cur = network_newick
    for _ in range(8):  # at most 8 expansion passes for safety
        def repl(m):
            label = m.group(1)
            sub = label_to_subtree.get(label)
            if sub is None:
                return label  # unknown synthetic, leave it
            sub = sub.strip().rstrip(";")
            return sub
        nxt = SYNTHETIC_LABEL_RE.sub(repl, cur)
        if nxt == cur:
            break
        cur = nxt
    return cur


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Expand synthetic leaves in a network back to the "
                    "subtrees recorded by compress_gtrees.py.")
    p.add_argument("input_network",
                   help="Path to the network file (one newick per line; "
                        "use - for stdin). Reticulation markers are preserved.")
    p.add_argument("-m", "--map", required=True,
                   help="Path to the compress_map.json written by "
                        "compress_gtrees.py.")
    p.add_argument("-o", "--out", default="-",
                   help="Output path (default: stdout).")
    p.add_argument("--sub-hc", choices=["off", "all", "concentrated"],
                   default="off",
                   help="Sub-HC strategy for soft clusters. 'off' (default): "
                        "use the modal-topology representative. 'all': run a "
                        "small SupNet HC on every soft cluster's fragments to "
                        "infer a mini reticulation network. 'concentrated': "
                        "run sub-HC only on clusters tagged 'soft-concentrated' "
                        "(requires compress_gtrees.py --retic-detection).")
    p.add_argument("--supnet", default="./supnet",
                   help="Path to the supnet binary used by sub-HC (default: "
                        "./supnet).")
    p.add_argument("--sub-R", type=int, default=2,
                   help="Reticulations per sub-HC run (default 2).")
    p.add_argument("--sub-cost", default="RF",
                   help="Cost for sub-HC (default RF).")
    p.add_argument("--sub-climb", type=int, default=50,
                   help="--hcstopclimb for sub-HC (default 50).")
    p.add_argument("--sub-prefix", default="H_",
                   help="Prefix for sub-HC reticulation IDs after renumeration. "
                        "Each sub-HC also gets a per-cluster suffix to keep IDs "
                        "globally unique. Default 'H_'.")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Print per-line replacement counts to stderr.")
    args = p.parse_args(argv)

    # Load map
    with open(args.map) as f:
        cmap = json.load(f)
    label_to_subtree = {
        entry["label"]: entry["representative_newick"]
        for entry in cmap.get("compressed", [])
    }
    if args.verbose:
        print(f"[expand] {len(label_to_subtree)} synthetic label(s) loaded "
              f"from {args.map}", file=sys.stderr)
        n_soft = sum(1 for e in cmap.get("compressed", [])
                     if e["kind"].startswith("soft"))
        print(f"[expand] sub-hc={args.sub_hc}, soft cluster(s)={n_soft}",
              file=sys.stderr)

    # Read input lines
    if args.input_network == "-":
        lines = sys.stdin.read().splitlines()
    else:
        with open(args.input_network) as f:
            lines = f.read().splitlines()

    # Sub-HC: per soft cluster, run a small supnet HC to derive a mini-network,
    # then renumber its reticulations to avoid colliding with main-network IDs.
    # We resolve sub-HC LAZILY per main-line so each main line gets its own
    # collision-free numbering. (Per-line because the same synthetic label may
    # appear with different surrounding retics across lines.)
    sub_hc_targets = set()
    if args.sub_hc != "off":
        if not shutil.which(args.supnet) and not os.path.exists(args.supnet):
            print(f"[expand] ERROR: supnet binary not found at "
                  f"{args.supnet!r}", file=sys.stderr)
            return 2
        for entry in cmap.get("compressed", []):
            kind = entry["kind"]
            if args.sub_hc == "all" and kind.startswith("soft"):
                sub_hc_targets.add(entry["label"])
            elif args.sub_hc == "concentrated" and kind == "soft-concentrated":
                sub_hc_targets.add(entry["label"])
        if args.verbose:
            print(f"[expand] sub-hc targets: {len(sub_hc_targets)} cluster(s)",
                  file=sys.stderr)
    label_to_entry = {e["label"]: e for e in cmap.get("compressed", [])}

    out_lines = []
    for ln in lines:
        if not ln.strip() or ln.startswith("#"):
            out_lines.append(ln)
            continue

        # Build per-line label_to_subtree, optionally running sub-HC for
        # targeted clusters. We iterate over labels actually present in this
        # line (cheap) so we don't run sub-HC on labels that don't appear.
        line_label_to_sub = dict(label_to_subtree)
        present = set(SYNTHETIC_LABEL_RE.findall(ln))
        used_ids = collect_main_retic_ids(ln)
        n_subhc = 0
        for label in present:
            if label not in sub_hc_targets:
                continue
            entry = label_to_entry.get(label, {})
            fragments = entry.get("fragments", [])
            if not fragments:
                if args.verbose:
                    print(f"[expand]   {label}: no fragments stored, "
                          f"falling back to representative", file=sys.stderr)
                continue
            sub_net = run_sub_hc(
                fragments, args.supnet, args.sub_R, args.sub_cost,
                args.sub_climb, randseed=7, verbose=args.verbose)
            if not sub_net:
                if args.verbose:
                    print(f"[expand]   {label}: sub-HC returned no result, "
                          f"falling back to representative", file=sys.stderr)
                continue
            # Per-cluster prefix uses a short slice of the synthetic label
            # to keep retic ids globally unique without bloating the newick.
            short = label.split("_cmp_", 1)[-1][:6] if "_cmp_" in label else label[:6]
            new_net, added_ids = renumber_subnet_retics(
                sub_net, used_ids, prefix=f"{args.sub_prefix}{short}_")
            used_ids |= added_ids
            line_label_to_sub[label] = new_net + ";"
            n_subhc += 1
            if args.verbose:
                print(f"[expand]   {label}: sub-HC ok, "
                      f"{len(added_ids)} retic(s) added", file=sys.stderr)

        n_before = len(SYNTHETIC_LABEL_RE.findall(ln))
        expanded = expand_one_network(ln, line_label_to_sub)
        n_after = len(SYNTHETIC_LABEL_RE.findall(expanded))
        if args.verbose:
            replaced = n_before - n_after
            print(f"[expand] line: replaced {replaced}/{n_before} synthetic "
                  f"leaves (sub-HC: {n_subhc})", file=sys.stderr)
            if n_after:
                print(f"[expand]   warning: {n_after} synthetic labels not "
                      f"resolved in this line", file=sys.stderr)
        out_lines.append(expanded)

    out_text = "\n".join(out_lines) + ("\n" if out_lines else "")
    if args.out == "-":
        sys.stdout.write(out_text)
    else:
        with open(args.out, "w") as f:
            f.write(out_text)
        if args.verbose:
            print(f"[expand] wrote {args.out}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
