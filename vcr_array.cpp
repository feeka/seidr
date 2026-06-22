// vcr_array.cpp — locked spec (PLAN.md §9): cross the ORF (gene) between VCRs.
//
// Pipeline per detected VCR:
//   anchor (indeg>=2 + upstream STOP)  --greedy CM walk-->  VCR, end k-mer
//   end k-mer  --find ATG within +-50 nodes (§9.2)-->       START (frame fixed)
//   START      --bounded BEAM traversal (§9.3-9.5)-->        cassette ORF(s)
//
// IMPORTANT (SDBG geometry): the succinct DBG stores edges ((k+1)-mers); a single
// edge-row's OutgoingEdges only follows THAT row's lineage, so it misses sibling
// branches at a fork. We therefore navigate by K-MER MEMBERSHIP instead: to extend
// a path we try each next base A/C/G/T and keep it iff that k-mer exists in the
// graph (IndexBinarySearch). This exposes every real branch at a node.
//
// Beam (§9.4): keep top 2000 partial paths ranked by moving-average multiplicity-
//   consistency (penalty = |mult - running_avg|, accumulated; lower is better),
//   with stop legality as a HARD GATE:
//     in-frame STOP before 210 bp -> discard; STOP in [210,1200] -> COMPLETE;
//     reach 1200 bp with no STOP   -> discard.
// Cycle rule (§9.5): a path revisiting k-mer N may continue only if outdeg(N) >= 2.
//
// This binary does NOT claim correctness — recovered ORFs are compared to the
// planted truth.txt by the test step (run.sh).
#include <sdbg/sdbg.h>
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
#include <array>
#include <set>
#include <unordered_set>
#include <algorithm>

// ---------- VCR CM profile (identical to vcr_traverse.cpp) ----------
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

// ---------- label helpers ----------
static int firstBase(const SDBG& g, uint64_t node) {        // GetLabel(node)[0], 1..4
    static thread_local std::vector<uint8_t> b; b.resize(g.k());
    g.GetLabel(node, b.data()); return b[0];
}
static int lastBase(const SDBG& g, uint64_t node) {         // GetLabel(node)[k-1], 1..4
    static thread_local std::vector<uint8_t> b; b.resize(g.k());
    g.GetLabel(node, b.data()); return b[g.k()-1];
}
static char baseChar(int code1) { return (code1>=1&&code1<=4) ? "ACGT"[code1-1] : 'N'; }
static int  baseCode(char c) {     // 1-based encoding used by GetLabel/IndexBinarySearch
    switch (c) { case 'A': return 1; case 'C': return 2; case 'G': return 3; case 'T': return 4; }
    return 0;
}
static std::string labelOf(const SDBG& g, uint64_t node) {     // full k-mer label, 5'->3'
    int K = g.k(); std::vector<uint8_t> b(K); g.GetLabel(node, b.data());
    std::string s; for (int i=0;i<K;++i) s.push_back(baseChar(b[i])); return s;
}
// look a k-mer up by sequence; returns edge id or kNullID
static uint64_t lookupKmer(const SDBG& g, const std::string& kmer) {
    int K = g.k(); if ((int)kmer.size() != K) return SDBG::kNullID;
    std::vector<uint8_t> seq(K);
    for (int i = 0; i < K; ++i) { int c = baseCode(kmer[i]); if (!c) return SDBG::kNullID; seq[i] = (uint8_t)c; }
    return g.IndexBinarySearch(seq.data());
}
// out-degree of a k-mer node = number of next bases whose k-mer exists (for §9.5)
static int kmerOutdeg(const SDBG& g, const std::string& kmer) {
    int d = 0; std::string nk = kmer.substr(1) + 'A';
    for (char b : std::string("ACGT")) { nk[nk.size()-1] = b;
        uint64_t e = lookupKmer(g, nk); if (e != SDBG::kNullID && g.IsValidEdge(e)) ++d; }
    return d;
}

// ---------- anchor detection (same geometry as vcr_traverse.cpp) ----------
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
static bool isStopCodon(const std::string& c){ return c=="TAA"||c=="TAG"||c=="TGA"; }

// ---------- greedy CM walk -> VCR; return the VCR's 3' end node ----------
static void traverseVCR(const SDBG& g, uint64_t anchor, const Profile& P, uint64_t& endNode) {
    int K = g.k();
    std::vector<uint8_t> seq(K); g.GetLabel(anchor, seq.data());
    int pos = K; uint64_t cur = anchor;
    while (pos < P.M()) {
        uint64_t outs[8]; int n = g.OutgoingEdges(cur, outs);
        int bestB = -1; uint64_t bestE = 0; double bestS = -1e18;
        for (int i = 0; i < n; ++i) {
            if (!g.IsValidEdge(outs[i])) continue;
            int b = lastBase(g, outs[i]) - 1;
            double s = P.s(pos, b);
            if (s > bestS) { bestS = s; bestB = b; bestE = outs[i]; }
        }
        if (bestB < 0) break;
        cur = bestE; ++pos;
    }
    endNode = cur;
}

// ---------- §9.2 find START: nearest ATG within <=50 nodes downstream ----------
// Navigates by k-mer membership (so it sees every branch). Returns the 23-mer that
// ends in ...ATG (the ORF start node), in `startKmer`.
static bool findStart(const SDBG& g, const std::string& vcrEndKmer, std::string& startKmer) {
    int K = g.k();
    struct St { std::string kmer; char pp; char p; int d; };
    std::vector<St> q;
    std::unordered_set<std::string> seen;
    q.push_back({vcrEndKmer, vcrEndKmer[K-2], vcrEndKmer[K-1], 0});
    seen.insert(vcrEndKmer);
    size_t head = 0;
    while (head < q.size()) {
        St s = q[head++];
        if (s.d > 50) continue;
        for (char b : std::string("ACGT")) {
            std::string nk = s.kmer.substr(1) + b;
            uint64_t e = lookupKmer(g, nk); if (e == SDBG::kNullID || !g.IsValidEdge(e)) continue;
            if (s.pp=='A' && s.p=='T' && b=='G') { startKmer = nk; return true; } // ATG completed
            if (seen.insert(nk).second) q.push_back({nk, s.p, b, s.d+1});
        }
    }
    return false;
}

// ---------- §9.3-9.5 bounded beam ORF traversal (k-mer navigation) ----------
struct Path {
    std::string kmer;                      // current k-mer (graph position)
    std::string seq;                       // ORF bases, includes the leading "ATG"
    double multsum; long multn; double devsum;        // moving-average bookkeeping
    std::unordered_set<std::string> visited;          // visited k-mers, for §9.5
};
static void crossORF(const SDBG& g, const std::string& startKmer, std::set<std::string>& recovered) {
    const size_t WIDTH = 2000;
    const int LMIN = 210, LMAX = 1200;
    Path init; init.kmer = startKmer; init.seq = "ATG";
    init.multsum = 0; init.multn = 0; init.devsum = 0; init.visited.insert(startKmer);
    std::vector<Path> beam; beam.push_back(std::move(init));
    while (!beam.empty()) {
        std::vector<Path> nxt;
        for (auto& path : beam) {
            if ((int)path.seq.size() >= LMAX) continue;          // §9.3 cap: no STOP in window
            for (char b : std::string("ACGT")) {
                std::string nk = path.kmer.substr(1) + b;
                uint64_t e = lookupKmer(g, nk);
                if (e == SDBG::kNullID || !g.IsValidEdge(e)) continue;       // not a real branch
                bool revisit = path.visited.count(nk) > 0;
                if (revisit && kmerOutdeg(g, nk) < 2) continue;             // §9.5 cycle rule
                double mult = (double)g.EdgeMultiplicity(e);
                double avg  = (path.multn > 0) ? path.multsum / path.multn : mult;
                double dev  = std::fabs(mult - avg);
                std::string nseq = path.seq; nseq.push_back(b);
                int len = (int)nseq.size();
                if (len % 3 == 0) {                                          // completed a codon
                    std::string codon = nseq.substr(len-3, 3);
                    if (isStopCodon(codon)) {                               // §9.4 stop legality gate
                        if (len < LMIN) continue;                          // premature -> discard
                        if (len <= LMAX) { recovered.insert(nseq); continue; }  // COMPLETE -> emit
                        continue;
                    }
                }
                Path np;
                np.kmer = std::move(nk); np.seq = std::move(nseq);
                np.multsum = path.multsum + mult; np.multn = path.multn + 1; np.devsum = path.devsum + dev;
                np.visited = path.visited; np.visited.insert(np.kmer);
                nxt.push_back(std::move(np));
            }
        }
        if (nxt.size() > WIDTH) {                                            // §9.4 keep top 2000 by lowest devsum
            std::nth_element(nxt.begin(), nxt.begin()+WIDTH, nxt.end(),
                             [](const Path& a, const Path& b){ return a.devsum < b.devsum; });
            nxt.resize(WIDTH);
        }
        beam = std::move(nxt);
    }
}

int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s <sdbg_prefix> <vcr_profile.txt> [out_orfs.fa]\n", argv[0]); return 1; }
    SDBG g; g.LoadFromFile(argv[1]);
    Profile P; if (!P.load(argv[2])) { fprintf(stderr, "cannot load profile %s\n", argv[2]); return 1; }
    fprintf(stderr, "[array] k=%d edges=%llu vcr_profile_len=%d\n", g.k(), (unsigned long long)g.size(), P.M());
    FILE* out = (argc >= 4) ? fopen(argv[3], "w") : stdout;

    std::set<std::string> recovered;
    int anchors = 0, vcr_walks = 0, starts = 0;
    for (uint64_t i = 0; i < g.size(); ++i) {
        if (!g.IsValidEdge(i) || g.EdgeIndegree(i) < 2) continue;
        if (!upstreamStop(g, i)) continue;
        ++anchors;
        uint64_t endNode; traverseVCR(g, i, P, endNode); ++vcr_walks;
        std::string startKmer;
        if (!findStart(g, labelOf(g, endNode), startKmer)) continue;
        ++starts;
        crossORF(g, startKmer, recovered);
    }
    int idx = 0;
    for (const auto& s : recovered) fprintf(out, ">orf%d len=%zu\n%s\n", idx++, s.size(), s.c_str());
    if (out != stdout) fclose(out);

    int nc = (int)recovered.size();
    fprintf(stderr, "[array] anchors=%d vcr_walks=%d starts_found=%d recovered_ORFs=%d\n",
            anchors, vcr_walks, starts, nc);
    if (nc < 2 || nc > 300) fprintf(stderr, "[array] NOTE: cassette count %d is outside [2,300]\n", nc);
    return 0;
}
