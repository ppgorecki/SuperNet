
// DP3 algorithm for the Robinson-Foulds score (Gorecki et al., RF Conflict
// Resolution). Computes max_xi rfsim(G, N, xi) via the paper's pair-based
// formulation: each entry stores <value, h, retusage>, where
//   value   = max_{h' ⪰ g} δ(g, s, h') (paper's ζ_1*),
//   h       = the lowest h' ⪰ g realising that max (paper's ζ_2*).
// Bonus +1 in the δ recursion is added when h = g, which by the
// cluster-to-lca-equality lemma corresponds exactly to C_xi(g) = C(g).
//
// Tie-break in pickmax prefers deeper h (chain-structure lemma guarantees
// candidates lie on a common root path so depth comparison is well-defined).
//
// LCA in G is O(1) via Euler tour + sparse table (Bender–Farach-Colton-style,
// here using sparse table for range-minimum, giving O(|G| log |G|) build and
// O(1) query). Total DP cost is the paper's O(|G||N|) per gene tree.

#ifndef _DP_RF_H
#define _DP_RF_H

#include "dag.h"
#include "rtree.h"
#include "network.h"
#include "bitcluster.h"

#define DPRF_ZETA   1
#define DPRF_ZETAU  2

class DP_RF
{
	RootedTree &genetree;
	Network &network;
	size_t sz;          // genetree.nn * network.nn

	// zeta(g, s) = <zeta_val, zeta_h, zeta_ret>
	COSTT    *zeta_val;
	NODEID   *zeta_h;
	RETUSAGE *zeta_ret;

	// zetaupup(g, s) = <zetau_val, zetau_h, zetau_ret>
	COSTT    *zetau_val;
	NODEID   *zetau_h;
	RETUSAGE *zetau_ret;

	int8_t *computed;          // bitmask of DPRF_ZETA / DPRF_ZETAU

	// gamma[s] = lca_G(C_t(s)) where C_t is the *tree* cluster (skipping
	// reticulation children, i.e. embretnet's Node.treecluster()). For a
	// reticulation s gamma[s] = MAXNODEID (the empty-cluster sentinel,
	// which lcag treats as identity). For nodes contracted away in
	// ContractedNetwork (mapdn[s] != s) this entry is never written.
	NODEID *gamma;
	int8_t *gamma_done;        // memoization flag (only "live" nodes set)

	// Euler tour + sparse table for O(1) LCA in genetree.
	NODEID *euler;        // length etlen, node id at each tour position
	int    *etdepth;      // length etlen, depth at each tour position
	int    *etfirst;      // length genetree.nn, first occurrence in euler
	int   **sparse;       // sparse[k][i] = position of min depth in
	                      // euler[i .. i + 2^k - 1]; k in [0..sparsek]
	int    *logtab;       // logtab[i] = floor(log2(i))
	int     etlen;
	int     sparsek;

	void buildrmq();
	void eulertour(NODEID v, int d, int &pos);
	NODEID rmqlca(NODEID a, NODEID b);

	inline size_t idx(NODEID g, NODEID s) const { return (size_t)g + (size_t)s * (size_t)genetree.nn; }

	NODEID computegamma(NODEID s);   // recursive, memoized
	void computezeta(NODEID g, NODEID s);
	void computezetau(NODEID g, NODEID s);

	NODEID lcag(NODEID a, NODEID b);

	void pickmax(COSTT &bestv, NODEID &besth, RETUSAGE &bestr,
	             COSTT cv, NODEID ch, const RETUSAGE &cr);

public:
	DP_RF(RootedTree &g, Network &n);
	~DP_RF();

	void preprocess();
	void clean();

	// Returns max_xi rfsim(G, N, xi); -INFTY if no scenario exists.
	// Fills retusage with the corresponding usage of reticulation edges.
	COSTT maxrfsim(RETUSAGE &retusage);

	// RF-score: 2|I(G)| - 2*max_rfsim + 2 (paper's lemma noname1 with
	// |I(G)| excluding the root). Returns INFTY if no scenario exists.
	COSTT minrfscore(RETUSAGE &retusage);
};

#endif
