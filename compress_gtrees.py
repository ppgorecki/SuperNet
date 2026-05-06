#!/usr/bin/env python3
"""
compress_gtrees.py — external preprocessor for SupNet's HC pipeline.

Scans a collection of gene trees, finds clusters that are well-supported
across the input set and have consistent topology in their supporting GTs,
and compresses each such subtree into a single synthetic leaf. Outputs
reduced gene trees plus a JSON map for later expansion.

Pipeline:
  GT_set ──(this tool)──▶ compressed GT_set + map.json
                                 │
                                 ▼
                              supnet HC
                                 │
                                 ▼
                          expand_network.py
                                 │
                                 ▼
                          full network

Cluster classification per internal node of the guide tree (auto-built from
strict consensus, or user-supplied):
  - rigid: ≥ rigid-threshold of GTs containing the cluster have the modal
           subtree topology. Compressed; sub-HC NOT needed.
  - soft:  ≥ soft-threshold but < rigid-threshold. Compressed; flagged for
           later expand-time sub-HC.
  - open:  below soft-threshold. Not compressed.

Run with `--help` for options.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from collections import Counter, defaultdict
from typing import Optional


# ---------------------------------------------------------------------------
# Newick parser/serializer (minimal, no support for branch lengths/comments
# beyond stripping; supports unrooted trees by treating top trifurcation as
# rooted at the implicit root). Trees are stored as nested tuples:
#   leaf:     ("leaf", "label")
#   internal: ("node", [children])
# Reticulation markers (#1, #ABC) are NOT supported in inputs (this tool is
# for tree-shaped GTs); the network output is handled by expand_network.py.
# ---------------------------------------------------------------------------

class ParseError(Exception):
    pass


def parse_newick(s: str) -> tuple:
    s = s.strip()
    if s.endswith(";"):
        s = s[:-1]
    pos = [0]
    tree = _parse_node(s, pos)
    if pos[0] != len(s):
        raise ParseError(f"trailing garbage at {pos[0]}: {s[pos[0]:]!r}")
    return tree


def _parse_node(s: str, pos: list) -> tuple:
    if pos[0] >= len(s):
        raise ParseError("unexpected end")
    if s[pos[0]] == "(":
        pos[0] += 1
        children = [_parse_node(s, pos)]
        while pos[0] < len(s) and s[pos[0]] == ",":
            pos[0] += 1
            children.append(_parse_node(s, pos))
        if pos[0] >= len(s) or s[pos[0]] != ")":
            raise ParseError(f"missing ) at {pos[0]}")
        pos[0] += 1
        # optional internal label / branch length — strip silently
        _skip_label(s, pos)
        return ("node", children)
    else:
        # leaf
        start = pos[0]
        while pos[0] < len(s) and s[pos[0]] not in "(),:;":
            pos[0] += 1
        label = s[start:pos[0]]
        # strip branch length
        if pos[0] < len(s) and s[pos[0]] == ":":
            pos[0] += 1
            while pos[0] < len(s) and s[pos[0]] not in "(),;":
                pos[0] += 1
        if not label:
            raise ParseError(f"empty leaf label at {start}")
        return ("leaf", label)


def _skip_label(s: str, pos: list) -> None:
    while pos[0] < len(s) and s[pos[0]] not in "(),;":
        pos[0] += 1


def serialize_newick(tree: tuple) -> str:
    """Render to newick (no branch lengths)."""
    if tree[0] == "leaf":
        return tree[1]
    return "(" + ",".join(serialize_newick(c) for c in tree[1]) + ")"


def canonical_newick(tree: tuple) -> str:
    """Canonical form: sort children by their canonical strings. Used as a
    topology hash key — two subtrees with the same canonical newick are
    isomorphic on the labelled-leaf level."""
    if tree[0] == "leaf":
        return tree[1]
    parts = sorted(canonical_newick(c) for c in tree[1])
    return "(" + ",".join(parts) + ")"


def leaves_of(tree: tuple) -> frozenset:
    if tree[0] == "leaf":
        return frozenset((tree[1],))
    s = set()
    for c in tree[1]:
        s |= leaves_of(c)
    return frozenset(s)


def all_internal_subtrees(tree: tuple):
    """Yield (subtree, leafset) for every internal node, root included."""
    if tree[0] == "leaf":
        return
    yield tree, leaves_of(tree)
    for c in tree[1]:
        if c[0] == "node":
            yield from all_internal_subtrees(c)


# ---------------------------------------------------------------------------
# Cluster analysis
# ---------------------------------------------------------------------------

def collect_subtree_index(gtrees: list) -> dict:
    """Index: cluster (frozenset of leaves) -> list of (gt_idx, subtree).
    Only internal subtrees are recorded (skip singletons)."""
    index = defaultdict(list)
    for i, gt in enumerate(gtrees):
        for sub, leaves in all_internal_subtrees(gt):
            if len(leaves) >= 2:
                index[leaves].append((i, sub))
    return index


def majority_topology_percentage(entries: list) -> tuple:
    """Given list of (gt_idx, subtree) for one cluster, compute:
       - modal_canon: canonical newick of the most common topology
       - modal_count: how many GTs have that topology
       - total: |entries|
       - representative: an example subtree with the modal topology
       - canon_counts: full Counter of canonical topologies (for downstream
                       concentration analysis).
    Returns (modal_canon, modal_count, total, representative, canon_counts)."""
    canons = [canonical_newick(sub) for _, sub in entries]
    counter = Counter(canons)
    modal_canon, modal_count = counter.most_common(1)[0]
    representative = next(sub for (_, sub), c in zip(entries, canons)
                          if c == modal_canon)
    return modal_canon, modal_count, len(entries), representative, counter


def classify_concentration(canon_counts: Counter, total: int,
                           min_alt_coverage: float, max_alts: int) -> str:
    """Given the topology-count distribution for a soft cluster, decide
    whether the alternatives are 'concentrated' (few specific outlier
    topologies → plausibly a real reticulation) or 'dispersed' (many small
    topologies → more likely ILS or noise). Returns 'soft-concentrated' or
    'soft-dispersed'.
    The modal topology itself is excluded from the alternatives count."""
    sorted_topos = canon_counts.most_common()
    # Skip [0] (modal). For the rest, count how many have ≥ min_alt_coverage.
    n_significant = sum(
        1 for canon, c in sorted_topos[1:]
        if c / total >= min_alt_coverage
    )
    return "soft-concentrated" if n_significant <= max_alts else "soft-dispersed"


def is_laminar(a: frozenset, b: frozenset) -> bool:
    """Two clusters are compatible (laminar) iff one is a subset of the
    other or they are disjoint."""
    return a.isdisjoint(b) or a <= b or b <= a


def greedy_max_laminar(candidates: list) -> list:
    """candidates is a list of (cluster, ...metadata...) sorted by some
    priority (we use cluster size descending — larger first preserves
    biggest savings). Returns the maximal laminar subset."""
    accepted = []
    accepted_clusters = []
    for entry in candidates:
        c = entry[0]
        if all(is_laminar(c, a) for a in accepted_clusters):
            accepted.append(entry)
            accepted_clusters.append(c)
    return accepted


# ---------------------------------------------------------------------------
# Compression
# ---------------------------------------------------------------------------

def synthetic_label(canon: str, prefix: str = "_cmp_") -> str:
    """Deterministic synthetic species label. Short hash of the canonical
    topology so two compressions of the same cluster collapse to the same
    label across runs."""
    h = hashlib.sha1(canon.encode("utf-8")).hexdigest()[:8]
    return f"{prefix}{h}"


def replace_subtree(tree: tuple, target_leafset: frozenset, replacement: tuple) -> tuple:
    """Replace the (highest) subtree whose leaf set equals target_leafset
    with `replacement`. Returns a new tree (no in-place mutation).
    If no such subtree exists, returns tree unchanged."""
    if tree[0] == "leaf":
        return tree
    if leaves_of(tree) == target_leafset:
        return replacement
    new_children = [replace_subtree(c, target_leafset, replacement)
                    for c in tree[1]]
    return ("node", new_children)


def compress_one_gt(gt: tuple, accepted: list) -> tuple:
    """Apply all accepted compressions to a single GT in cluster-size
    descending order (largest first to avoid nested replacements)."""
    for cluster, label, *_ in accepted:
        gt = replace_subtree(gt, cluster, ("leaf", label))
    return gt


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main(argv=None):
    p = argparse.ArgumentParser(
        description="Compress similar subtrees in a gene-tree collection. "
                    "Outputs reduced GTs plus a JSON map for later expansion.")
    p.add_argument("input_gts",
                   help="Path to gene-tree file (one newick per line). Use - for stdin.")
    p.add_argument("-o", "--out-gts", default="compressed_gts.txt",
                   help="Output reduced GT file (default: compressed_gts.txt).")
    p.add_argument("-m", "--out-map", default="compress_map.json",
                   help="Output map file (default: compress_map.json).")
    p.add_argument("--rigid-threshold", type=float, default=0.95,
                   help="Min fraction of supporting GTs that must agree on "
                        "the modal topology to mark a cluster rigid (default 0.95).")
    p.add_argument("--soft-threshold", type=float, default=0.70,
                   help="Min fraction for 'soft' classification "
                        "(below rigid). Default 0.70.")
    p.add_argument("--cluster-freq", type=float, default=1.0,
                   help="Min fraction of GTs that must contain the cluster "
                        "as an internal cluster (default 1.0). "
                        "RF cost requires bijective leaf labelling which means "
                        "all GTs must end up with the same leaf set after "
                        "compression — keep this at 1.0 unless using a cost "
                        "(DC/D/DL) that tolerates missing leaves.")
    p.add_argument("--min-cluster", type=int, default=3,
                   help="Skip clusters smaller than this (default 3).")
    p.add_argument("--prefix", default="_cmp_",
                   help="Prefix for synthetic species labels.")
    p.add_argument("--retic-detection", action="store_true",
                   help="Mark each soft cluster as 'soft-concentrated' or "
                        "'soft-dispersed' based on the spread of alternative "
                        "topologies. Concentrated (≤2 alt-topologies with "
                        "support ≥ --concentration-min) is plausibly a "
                        "reticulation; dispersed is more likely ILS/noise. "
                        "expand_network.py --sub-hc=concentrated then targets "
                        "only the concentrated ones.")
    p.add_argument("--concentration-min", type=float, default=0.10,
                   help="Min coverage for an alternative topology to count "
                        "as 'significant' in the concentration test "
                        "(default 0.10, i.e. 10%% of supporting GTs).")
    p.add_argument("--concentration-max-alts", type=int, default=2,
                   help="Max number of significant alternatives for a soft "
                        "cluster to be classified concentrated (default 2).")
    p.add_argument("-v", "--verbose", action="store_true",
                   help="Print diagnostics to stderr.")
    args = p.parse_args(argv)

    # Read input GTs
    if args.input_gts == "-":
        lines = sys.stdin.read().splitlines()
    else:
        with open(args.input_gts) as f:
            lines = f.read().splitlines()
    gtrees = []
    for ln in lines:
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        try:
            gtrees.append(parse_newick(ln))
        except ParseError as e:
            print(f"parse error in line: {ln!r}: {e}", file=sys.stderr)
            return 1
    if not gtrees:
        print("no gene trees parsed", file=sys.stderr)
        return 1

    n_gts = len(gtrees)
    if args.verbose:
        all_leaves = set()
        for g in gtrees:
            all_leaves |= leaves_of(g)
        print(f"[compress] {n_gts} GTs, {len(all_leaves)} distinct species",
              file=sys.stderr)

    # 1. Index every internal subtree by leafset
    index = collect_subtree_index(gtrees)

    # Universe of all species across all GTs. We never compress a cluster
    # equal to U — that would replace whole gene trees with a single leaf
    # and erase the structure HC is supposed to learn.
    universe = set()
    for g in gtrees:
        universe |= leaves_of(g)

    # 2. For each cluster: classify
    #    cluster-freq = entries / n_gts
    #    agreement   = modal_count / entries
    candidates = []  # (cluster, label, kind, modal_canon, representative)
    for cluster, entries in index.items():
        if len(cluster) < args.min_cluster:
            continue
        if cluster == universe:
            continue  # don't collapse the entire tree
        freq = len(entries) / n_gts
        if freq < args.cluster_freq:
            continue
        modal_canon, modal_count, total, rep, canon_counts = \
            majority_topology_percentage(entries)
        agreement = modal_count / total
        if agreement >= args.rigid_threshold:
            kind = "rigid"
        elif agreement >= args.soft_threshold:
            kind = "soft"
            if args.retic_detection:
                kind = classify_concentration(
                    canon_counts, total,
                    args.concentration_min, args.concentration_max_alts)
        else:
            continue  # open — leave uncompressed
        label = synthetic_label(modal_canon, prefix=args.prefix)
        # Save fragments only for non-rigid kinds (rigid: all isomorphic by
        # definition, representative is enough). Fragments are needed for
        # sub-HC at expand time.
        fragments = []
        if kind != "rigid":
            fragments = [serialize_newick(sub) + ";" for _, sub in entries]
        candidates.append({
            "cluster": cluster,
            "label": label,
            "kind": kind,
            "modal_canon": modal_canon,
            "representative": rep,
            "fragments": fragments,
            "freq": freq,
            "agreement": agreement,
            "size": len(cluster),
            "topology_counts": dict(canon_counts.most_common()),
        })

    # 3. Greedy laminar — large clusters first
    candidates.sort(key=lambda e: -e["size"])
    laminar = greedy_max_laminar([
        (c["cluster"], c["label"], c) for c in candidates
    ])
    accepted = [meta for _, _, meta in laminar]

    if args.verbose:
        kind_counts = Counter(a["kind"] for a in accepted)
        kinds_summary = ", ".join(f"{k}={v}" for k, v in kind_counts.items())
        print(f"[compress] candidates={len(candidates)} "
              f"accepted={len(accepted)} ({kinds_summary})", file=sys.stderr)
        for a in accepted:
            print(f"  {a['kind']:18} size={a['size']:3} "
                  f"freq={a['freq']:.2f} agreement={a['agreement']:.2f} "
                  f"-> {a['label']}", file=sys.stderr)

    # 4. Compress each GT (largest accepted first to avoid nested overwrite)
    accepted_for_compress = [(a["cluster"], a["label"]) for a in accepted]
    accepted_for_compress.sort(key=lambda e: -len(e[0]))
    out_gts = [compress_one_gt(g, accepted_for_compress) for g in gtrees]

    # 5. Write outputs
    with open(args.out_gts, "w") as f:
        for g in out_gts:
            f.write(serialize_newick(g) + ";\n")
    if args.verbose:
        print(f"[compress] wrote {args.out_gts}", file=sys.stderr)

    cmap = {
        "version": 1,
        "rigid_threshold": args.rigid_threshold,
        "soft_threshold": args.soft_threshold,
        "min_cluster": args.min_cluster,
        "input_gt_count": n_gts,
        "retic_detection": args.retic_detection,
        "compressed": [
            {
                "label": a["label"],
                "kind": a["kind"],
                "size": a["size"],
                "freq": a["freq"],
                "agreement": a["agreement"],
                "leaves": sorted(a["cluster"]),
                "representative_newick": serialize_newick(a["representative"]) + ";",
                "modal_canon": a["modal_canon"],
                "fragments": a["fragments"],
                "topology_counts": a["topology_counts"],
            }
            for a in accepted
        ],
    }
    with open(args.out_map, "w") as f:
        json.dump(cmap, f, indent=2)
    if args.verbose:
        print(f"[compress] wrote {args.out_map}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    sys.exit(main())
