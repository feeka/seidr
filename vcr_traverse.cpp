// vcr_traverse.cpp — Step 5 on the real SDBG.
// anchors (indeg>=2 + upstream STOP) -> greedy profile-guided walk -> candidate VCR.
//   weight = linear CM profile (per-position log-odds, from cm_to_profile.py).
//   at each step pick the GRAPH successor whose added base maximises the profile
//   score at the current position (greedy argmax OVER GRAPH BRANCHES), advance
//   the profile position, stop when the profile is consumed.
// Candidates are written to FASTA for cmsearch validation; this binary does NOT
// validate or claim correctness — that is done by cmsearch + identity check.
#include <sdbg/sdbg.h>
#include <cstdio>
#include <vector>
#include <string>
#include <array>

// ---- profile: per-position log-odds (A C G T), from cm_to_profile.py ----
struct Profile {
    std::vector<std::array<double,4>> lo;
    bool load(const char* path) {
        FILE* f = fopen(path, "r"); if (!f) return false;
        int p; double a,c,g,t;
        while (fscanf(f, "%d %lf %lf %lf %lf", &p, &a, &c, &g, &t) == 5) lo.push_back({a,c,g,t});
        fclose(f); return !lo.empty();
    }
    int M() const { return (int)lo.size(); }
    double s(int pos, int base) const {
        return (pos<0||pos>=M()||base<0||base>3) ? -4.0 : lo[pos][base];
    }
};

// ---- anchor detection (same geometry as vcr_anchor.cpp) ----
static int firstBase(const SDBG& g, uint64_t node) {        // GetLabel(node)[0], 1..4
    static thread_local std::vector<uint8_t> b; b.resize(g.k());
    g.GetLabel(node, b.data()); return b[0];
}
static int lastBase(const SDBG& g, uint64_t node) {         // GetLabel(node)[k-1], 1..4
    static thread_local std::vector<uint8_t> b; b.resize(g.k());
    g.GetLabel(node, b.data()); return b[g.k()-1];
}
static const char* stopName(int b0,int b1,int b2){
    if(b0==4&&b1==1&&b2==1) return "TAA";
    if(b0==4&&b1==1&&b2==3) return "TAG";
    if(b0==4&&b1==3&&b2==1) return "TGA";
    return nullptr;
}
static const char* upstreamStop(const SDBG& g, uint64_t v) {
    uint64_t i1[8],i2[8],i3[8];
    int n1=g.IncomingEdges(v,i1);
    for(int a=0;a<n1;++a){ if(!g.IsValidEdge(i1[a]))continue; int b2=firstBase(g,i1[a]);
      int n2=g.IncomingEdges(i1[a],i2);
      for(int b=0;b<n2;++b){ if(!g.IsValidEdge(i2[b]))continue; int b1=firstBase(g,i2[b]);
        int n3=g.IncomingEdges(i2[b],i3);
        for(int c=0;c<n3;++c){ if(!g.IsValidEdge(i3[c]))continue; int b0=firstBase(g,i3[c]);
          if(const char* s=stopName(b0,b1,b2)) return s; } } }
    return nullptr;
}

// ---- greedy profile-guided walk from an anchor ----
static std::string traverse(const SDBG& g, uint64_t anchor, const Profile& P, double& score) {
    int K = g.k();
    std::vector<uint8_t> seq(K); g.GetLabel(anchor, seq.data());
    std::string path; path.reserve(P.M());
    double sc = 0;
    for (int j = 0; j < K; ++j) { int b = seq[j]-1; path.push_back((b>=0&&b<4)?"ACGT"[b]:'N'); sc += P.s(j,b); }
    int pos = K; uint64_t cur = anchor;
    while (pos < P.M()) {
        uint64_t outs[8]; int n = g.OutgoingEdges(cur, outs);
        int bestB = -1; uint64_t bestE = 0; double bestS = -1e18;
        for (int i = 0; i < n; ++i) {
            if (!g.IsValidEdge(outs[i])) continue;
            int b = lastBase(g, outs[i]) - 1;          // 0..3
            double s = P.s(pos, b);
            if (s > bestS) { bestS = s; bestB = b; bestE = outs[i]; }
        }
        if (bestB < 0) break;                          // dead end
        path.push_back("ACGT"[bestB]); sc += bestS; cur = bestE; ++pos;
    }
    score = sc; return path;
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s <sdbg_prefix> <profile.txt> [out.fa]\n", argv[0]); return 1; }
    SDBG g; g.LoadFromFile(argv[1]);
    Profile P; if (!P.load(argv[2])) { fprintf(stderr, "cannot load profile %s\n", argv[2]); return 1; }
    fprintf(stderr, "[traverse] k=%d edges=%llu profile_len=%d\n", g.k(), (unsigned long long)g.size(), P.M());
    FILE* out = (argc >= 4) ? fopen(argv[3], "w") : stdout;
    int n = 0;
    for (uint64_t i = 0; i < g.size(); ++i) {
        if (!g.IsValidEdge(i) || g.EdgeIndegree(i) < 2) continue;
        const char* st = upstreamStop(g, i);
        if (!st) continue;
        double sc; std::string vcr = traverse(g, i, P, sc);
        fprintf(out, ">cand%d anchor=%llu upstream=%s len=%zu score=%.1f\n%s\n",
                ++n, (unsigned long long)i, st, vcr.size(), sc, vcr.c_str());
    }
    if (out != stdout) fclose(out);
    fprintf(stderr, "[traverse] candidates: %d\n", n);
    return 0;
}
