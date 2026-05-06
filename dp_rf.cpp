
#include "dp_rf.h"
#include "tools.h"
#include "rtree.h"
#include "network.h"
#include "bitcluster.h"

#include <cstring>

DP_RF::DP_RF(RootedTree &g, Network &n) : genetree(g), network(n)
{
	sz = (size_t)genetree.nn * (size_t)network.nn;
	zeta_val  = new COSTT[sz];
	zeta_h    = new NODEID[sz];
	zeta_ret  = new RETUSAGE[sz];
	zetau_val = new COSTT[sz];
	zetau_h   = new NODEID[sz];
	zetau_ret = new RETUSAGE[sz];
	computed  = new int8_t[sz];
	gamma     = new NODEID[network.nn];
	gamma_done = new int8_t[network.nn];
	euler     = NULL;
	etdepth   = NULL;
	etfirst   = NULL;
	sparse    = NULL;
	logtab    = NULL;
	etlen     = 0;
	sparsek   = 0;
	clean();
}

DP_RF::~DP_RF()
{
	delete[] zeta_val;
	delete[] zeta_h;
	delete[] zeta_ret;
	delete[] zetau_val;
	delete[] zetau_h;
	delete[] zetau_ret;
	delete[] computed;
	delete[] gamma;
	delete[] gamma_done;
	delete[] euler;
	delete[] etdepth;
	delete[] etfirst;
	if (sparse)
	{
		for (int k = 0; k <= sparsek; k++) delete[] sparse[k];
		delete[] sparse;
	}
	delete[] logtab;
}

void DP_RF::clean()
{
	memset(computed,   0, sizeof(int8_t) * sz);
	memset(gamma_done, 0, sizeof(int8_t) * network.nn);
}

// Recursive memoized lca_G(C_t(s)). Uses Network's virtual getleftchild /
// getrightchild / getretchild — for ContractedNetwork these route through
// `mapdn` so we never visit nodes that were contracted away.
NODEID DP_RF::computegamma(NODEID s)
{
	if (gamma_done[s]) return gamma[s];
	gamma_done[s] = 1;
	if (s < network.lf)
	{
		gamma[s] = genetree.findlab(network.lab[s], 0);
	}
	else if (s >= network.rtstartid)
	{
		gamma[s] = MAXNODEID;
	}
	else
	{
		NODEID s0 = network.getleftchild(s);
		NODEID s1 = network.getrightchild(s);
		gamma[s] = lcag(computegamma(s0), computegamma(s1));
	}
	return gamma[s];
}

void DP_RF::preprocess()
{
	genetree.initdepth();
	buildrmq();
	clean();
	// gamma is computed lazily on first access (computegamma(s)) — this is
	// robust to ContractedNetwork's id layout where live nodes may not all
	// be reachable from network.root via a single descent.
}

// Euler-tour DFS: writes euler[pos] = current node id and etdepth[pos] = d at
// every entry/return; etfirst[v] is the first time v is touched.
void DP_RF::eulertour(NODEID v, int d, int &pos)
{
	etfirst[v] = pos;
	euler[pos] = v;
	etdepth[pos] = d;
	pos++;
	if (v >= genetree.lf)
	{
		eulertour(genetree.leftchild[v],  d + 1, pos);
		euler[pos] = v;
		etdepth[pos] = d;
		pos++;
		eulertour(genetree.rightchild[v], d + 1, pos);
		euler[pos] = v;
		etdepth[pos] = d;
		pos++;
	}
}

void DP_RF::buildrmq()
{
	int n = genetree.nn;
	etlen = 2 * n - 1;        // exact for a binary tree
	euler   = new NODEID[etlen];
	etdepth = new int[etlen];
	etfirst = new int[n];

	int pos = 0;
	eulertour(genetree.root, 0, pos);
	// pos == etlen; if not (degenerate tree), recompute etlen — keeps query
	// arithmetic correct on edge cases.
	etlen = pos;

	// floor-log2 table
	logtab = new int[etlen + 1];
	logtab[0] = 0;
	logtab[1] = 0;
	for (int i = 2; i <= etlen; i++) logtab[i] = logtab[i / 2] + 1;

	sparsek = logtab[etlen];
	sparse = new int*[sparsek + 1];
	for (int k = 0; k <= sparsek; k++) sparse[k] = new int[etlen];

	// Each cell stores an index into euler[] (not a depth) so we can decode
	// the LCA node directly.
	for (int i = 0; i < etlen; i++) sparse[0][i] = i;

	for (int k = 1; k <= sparsek; k++)
	{
		int half = 1 << (k - 1);
		int span = 1 << k;
		for (int i = 0; i + span <= etlen; i++)
		{
			int a = sparse[k - 1][i];
			int b = sparse[k - 1][i + half];
			sparse[k][i] = (etdepth[a] <= etdepth[b]) ? a : b;
		}
	}
}

NODEID DP_RF::rmqlca(NODEID a, NODEID b)
{
	if (a == b) return a;
	int l = etfirst[a];
	int r = etfirst[b];
	if (l > r) { int t = l; l = r; r = t; }
	int k = logtab[r - l + 1];
	int span = 1 << k;
	int x = sparse[k][l];
	int y = sparse[k][r - span + 1];
	int idxmin = (etdepth[x] <= etdepth[y]) ? x : y;
	return euler[idxmin];
}

NODEID DP_RF::lcag(NODEID a, NODEID b)
{
	if (a == MAXNODEID) return b;
	if (b == MAXNODEID) return a;
	return rmqlca(a, b);
}

void DP_RF::pickmax(COSTT &bestv, NODEID &besth, RETUSAGE &bestr,
                    COSTT cv, NODEID ch, const RETUSAGE &cr)
{
	if (cv > bestv)
	{
		bestv = cv; besth = ch; bestr = cr;
		return;
	}
	if (cv == bestv && cv > -INFTY)
	{
		// tie-break: prefer the deeper h. By the chain-structure lemma the
		// candidates lie on a common root-to-g path, so depth comparison is
		// well-defined.
		bool deeper_new;
		if (besth == MAXNODEID) deeper_new = (ch != MAXNODEID);
		else if (ch == MAXNODEID) deeper_new = false;
		else deeper_new = (genetree.depth[ch] > genetree.depth[besth]);
		if (deeper_new) { bestv = cv; besth = ch; bestr = cr; }
	}
}

void DP_RF::computezeta(NODEID g, NODEID s)
{
	size_t i = idx(g, s);
	if (computed[i] & DPRF_ZETA) return;
	computed[i] |= DPRF_ZETA;

	bool gleaf = g < genetree.lf;
	bool sleaf = s < network.lf;
	bool sret  = s >= network.rtstartid;

	if (gleaf && sleaf)
	{
		if (genetree.lab[g] == network.lab[s])
		{
			zeta_val[i] = 0;
			zeta_h[i]   = g;
			emptyretusage(zeta_ret[i]);
		}
		else
		{
			zeta_val[i] = -INFTY;
			zeta_h[i]   = MAXNODEID;
			emptyretusage(zeta_ret[i]);
		}
		return;
	}

	if (gleaf || sleaf || sret)
	{
		zeta_val[i] = -INFTY;
		zeta_h[i]   = MAXNODEID;
		emptyretusage(zeta_ret[i]);
		return;
	}

	// g internal, s tree-internal — paper's four cases
	NODEID gp  = genetree.leftchild[g];
	NODEID gpp = genetree.rightchild[g];
	NODEID sp  = network.getleftchild(s);
	NODEID spp = network.getrightchild(s);

	// theta(s): mark direct reticulation children of s
	bool spret  = (sp  >= network.rtstartid);
	bool sppret = (spp >= network.rtstartid);
	RETUSAGE theta;
	emptyretusage(theta);
	if (spret)
	{
		NODEID rtid = sp - network.rtstartid;
		if (network.getparent(sp) == s) addleftretusage (theta, rtid);
		else                            addrightretusage(theta, rtid);
	}
	if (sppret)
	{
		NODEID rtid = spp - network.rtstartid;
		if (network.getparent(spp) == s) addleftretusage (theta, rtid);
		else                             addrightretusage(theta, rtid);
	}

	COSTT    bestv = -INFTY;
	NODEID   besth = MAXNODEID;
	RETUSAGE bestr; emptyretusage(bestr);

	// Combine children's <value, h, retusage>; bonus is paper's [g = h].
	auto tryCase = [&](COSTT vL, NODEID hL, const RETUSAGE &rL,
	                   COSTT vR, NODEID hR, const RETUSAGE &rR,
	                   bool addtheta)
	{
		if (vL <= -INFTY || vR <= -INFTY) return;
		NODEID h = lcag(hL, hR);
		COSTT v = vL + vR + (h == g ? 1 : 0);
		RETUSAGE r;
		unionretusage(r, rL, rR);
		if (addtheta) unionretusage(r, r, theta);
		pickmax(bestv, besth, bestr, v, h, r);
	};

	computezetau(gp,  sp);
	computezetau(gpp, spp);
	computezetau(gp,  spp);
	computezetau(gpp, sp);
	computezeta (gp,  s);
	computezetau(gpp, s);
	computezeta (gpp, s);
	computezetau(gp,  s);

	size_t ia = idx(gp,  sp),  ib  = idx(gpp, spp);
	size_t ic = idx(gp,  spp), id  = idx(gpp, sp);
	size_t ie = idx(gp,  s),   ife = idx(gpp, s);
	size_t ig_= idx(gpp, s),   ih  = idx(gp,  s);

	tryCase(zetau_val[ia], zetau_h[ia], zetau_ret[ia],
	        zetau_val[ib], zetau_h[ib], zetau_ret[ib], true);   // C1
	tryCase(zetau_val[ic], zetau_h[ic], zetau_ret[ic],
	        zetau_val[id], zetau_h[id], zetau_ret[id], true);   // C2
	tryCase(zeta_val[ie],  zeta_h[ie],  zeta_ret[ie],
	        zetau_val[ife], zetau_h[ife], zetau_ret[ife], false); // C3
	tryCase(zeta_val[ig_], zeta_h[ig_], zeta_ret[ig_],
	        zetau_val[ih], zetau_h[ih],  zetau_ret[ih], false); // C4

	zeta_val[i] = bestv;
	zeta_h[i]   = besth;
	zeta_ret[i] = bestr;
}

void DP_RF::computezetau(NODEID g, NODEID s)
{
	size_t i = idx(g, s);
	if (computed[i] & DPRF_ZETAU) return;
	computed[i] |= DPRF_ZETAU;

	if (s >= network.rtstartid)
	{
		// reticulation: pass-through to the unique child; the ret edge is
		// marked by the caller (theta / thetabar at the tree-node above).
		NODEID sc = network.getretchild(s);
		computezetau(g, sc);
		size_t ic = idx(g, sc);
		zetau_val[i] = zetau_val[ic];
		zetau_h[i]   = zetau_h[ic];
		zetau_ret[i] = zetau_ret[ic];
		return;
	}

	if (s < network.lf)
	{
		computezeta(g, s);
		zetau_val[i] = zeta_val[i];
		zetau_h[i]   = zeta_h[i];
		zetau_ret[i] = zeta_ret[i];
		return;
	}

	// s tree-node
	NODEID sp  = network.getleftchild(s);
	NODEID spp = network.getrightchild(s);

	computezetau(g, sp);
	computezetau(g, spp);
	computezeta (g, s);

	size_t isp  = idx(g, sp);
	size_t ispp = idx(g, spp);
	size_t izet = i;

	COSTT    bestv = -INFTY;
	NODEID   besth = MAXNODEID;
	RETUSAGE bestr; emptyretusage(bestr);

	// thetabar(s, descended_child, sibling): mark direct ret edge on the
	// descent side, and the sibling reticulation's "other parent" edge
	// when the sibling is a reticulation.
	auto thetabar = [&](NODEID child, NODEID sibling) {
		RETUSAGE r;
		emptyretusage(r);
		if (child >= network.rtstartid)
		{
			NODEID rtid = child - network.rtstartid;
			if (network.getparent(child) == s) addleftretusage(r, rtid);
			else                                addrightretusage(r, rtid);
		}
		if (sibling >= network.rtstartid)
		{
			NODEID rtid = sibling - network.rtstartid;
			if (network.getparent(sibling) == s) addrightretusage(r, rtid);
			else                                  addleftretusage (r, rtid);
		}
		return r;
	};

	// Cand stop-at-s: zeta(g, s)
	if (zeta_val[izet] > -INFTY)
		pickmax(bestv, besth, bestr, zeta_val[izet], zeta_h[izet], zeta_ret[izet]);

	// Cand descend-through-s' — h = lca(zetau_h[g, s'], gamma(s''))
	if (zetau_val[isp] > -INFTY)
	{
		NODEID h = lcag(zetau_h[isp], computegamma(spp));
		RETUSAGE r, tb = thetabar(sp, spp);
		unionretusage(r, zetau_ret[isp], tb);
		pickmax(bestv, besth, bestr, zetau_val[isp], h, r);
	}

	// Cand descend-through-s'' (symmetric)
	if (zetau_val[ispp] > -INFTY)
	{
		NODEID h = lcag(zetau_h[ispp], computegamma(sp));
		RETUSAGE r, tb = thetabar(spp, sp);
		unionretusage(r, zetau_ret[ispp], tb);
		pickmax(bestv, besth, bestr, zetau_val[ispp], h, r);
	}

	zetau_val[i] = bestv;
	zetau_h[i]   = besth;
	zetau_ret[i] = bestr;
}

COSTT DP_RF::maxrfsim(RETUSAGE &retusage)
{
	NODEID rg = genetree.root;
	NODEID rn = network.root;
	computezetau(rg, rn);
	size_t i = idx(rg, rn);
	retusage = zetau_ret[i];
	return zetau_val[i];
}

COSTT DP_RF::minrfscore(RETUSAGE &retusage)
{
	COSTT rfsim = maxrfsim(retusage);
	if (rfsim <= -INFTY) return INFTY;
	// Paper's |I(G)| excludes the root; rfsim sums over T(G) (incl. root,
	// which always matches), so:  min RFs = 2|I(G)| - 2*rfsim + 2.
	long ig = (long)genetree.nn - (long)genetree.lf - 1;
	return 2 * (COSTT)ig - 2 * rfsim + 2;
}
