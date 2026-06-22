// vcr_anchor.cpp — Step 3 of the plan, on MEGAHIT's real SDBG.
// Find anchors = edges with indegree >= 2 whose 3 UPSTREAM nodes spell a STOP
// codon, using the exact label geometry read from sdbg.h.
//
// Geometry (from sdbg.h GetLabel): seq[0] is the 5'-most base; walking backward
// fills lower indices. For node v and an incoming edge p, GetLabel(p)[0] is the
// base immediately upstream of v. So:
//   1st IncomingEdges -> base -> codon position seq[2] (3', closest)
//   2nd IncomingEdges -> seq[1]
//   3rd IncomingEdges -> seq[0] (5', furthest)
// codon read 5'->3' is seq[0]seq[1]seq[2]; check vs TAA/TAG/TGA.
// Encoding (1-based): A=1 C=2 G=3 T=4 ; TAA=(4,1,1) TAG=(4,1,3) TGA=(4,3,1).
#include <sdbg/sdbg.h>
#include <cstdio>
#include <vector>
#include <string>

static int firstBase(const SDBG& g, uint64_t node) {   // GetLabel(node)[0], 1..4 or 0
    static thread_local std::vector<uint8_t> buf;
    buf.resize(g.k());
    g.GetLabel(node, buf.data());
    return buf[0];
}
static const char* stopName(int b0, int b1, int b2) {
    if (b0==4 && b1==1 && b2==1) return "TAA";
    if (b0==4 && b1==1 && b2==3) return "TAG";
    if (b0==4 && b1==3 && b2==1) return "TGA";
    return nullptr;
}
// any incoming 3-step path that spells a stop -> return its name, else nullptr.
static const char* upstreamStop(const SDBG& g, uint64_t v) {
    uint64_t in1[8], in2[8], in3[8];
    int n1 = g.IncomingEdges(v, in1);
    for (int a = 0; a < n1; ++a) {
        if (!g.IsValidEdge(in1[a])) continue;
        int b2 = firstBase(g, in1[a]);             // seq[2] (closest)
        int n2 = g.IncomingEdges(in1[a], in2);
        for (int b = 0; b < n2; ++b) {
            if (!g.IsValidEdge(in2[b])) continue;
            int b1 = firstBase(g, in2[b]);         // seq[1]
            int n3 = g.IncomingEdges(in2[b], in3);
            for (int c = 0; c < n3; ++c) {
                if (!g.IsValidEdge(in3[c])) continue;
                int b0 = firstBase(g, in3[c]);     // seq[0] (furthest)
                if (const char* s = stopName(b0, b1, b2)) return s;
            }
        }
    }
    return nullptr;
}

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s <sdbg_prefix>\n", argv[0]); return 1; }
    SDBG g; g.LoadFromFile(argv[1]);
    int K = g.k();
    std::vector<uint8_t> seq(K);
    unsigned anchors = 0;
    for (uint64_t i = 0; i < g.size(); ++i) {
        if (!g.IsValidEdge(i)) continue;
        if (g.EdgeIndegree(i) < 2) continue;
        const char* st = upstreamStop(g, i);
        if (!st) continue;
        ++anchors;
        g.GetLabel(i, seq.data());
        std::string lab(K, 'N');
        for (int j = 0; j < K; ++j) lab[j] = (seq[j]>=1 && seq[j]<=4) ? "ACGT"[seq[j]-1] : 'N';
        printf("ANCHOR edge=%llu indeg=%d upstream=%s kmer=%s\n",
               (unsigned long long)i, g.EdgeIndegree(i), st, lab.c_str());
    }
    fprintf(stderr, "[vcr_anchor] anchors (indeg>=2 + upstream STOP): %u\n", anchors);
    return 0;
}
