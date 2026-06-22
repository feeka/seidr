// vcr_fold.cpp — idea B (PLAN.md §10): CM PAIR-emission inverted-repeat search on the SDBG.
//
// Forward beam fold-walk: align graph paths to the CM consensus columns, scoring base
// PAIRS via a stack — push the 5'-pair base (P5), pop + score the pair at its 3' partner
// (P3). This rewards a base that still PAIRS with its partner even where the sequence
// diverged ("symmetry, not sequence"), which is exactly what the flattened-CM greedy walk
// (vcr_traverse) cannot do. Candidates are written to FASTA for cmsearch validation.
//
// v1 = MATCH-ONLY (substitutions; NO indels). Tests the core hypothesis before adding
// insert/delete. Pair model produced by cm_to_pairmodel.py.
#include <sdbg/sdbg.h>
#include <cstdio>
#include <vector>
#include <string>
#include <array>
#include <set>
#include <algorithm>

// ---------- pair model (from cm_to_pairmodel.py) ----------
struct PairModel {
    int N = 0;
    std::vector<int> role;                    // 1..N : 0=SL 1=SR 2=P5 3=P3
    std::vector<int> partner;
    std::vector<std::array<double,4>>  sing;  // SL/SR/P5 : 4 log-odds (A C G T)
    std::vector<std::array<double,16>> cond;  // P3       : 16 conditional log-odds [x*4+y]
    bool load(const char* path) {
        FILE* f = fopen(path, "r"); if (!f) return false;
        std::vector<int> cols, roles, parts; std::vector<std::vector<double>> vals;
        int col, part; char rolebuf[16];
        while (fscanf(f, "%d %15s %d", &col, rolebuf, &part) == 3) {
            std::string r = rolebuf;
            int nv = (r == "P3") ? 16 : 4;
            std::vector<double> v(nv);
            bool okrow = true;
            for (int i = 0; i < nv; ++i) if (fscanf(f, "%lf", &v[i]) != 1) { okrow = false; break; }
            if (!okrow) { fclose(f); return false; }
            int rc = (r=="SL")?0 : (r=="SR")?1 : (r=="P5")?2 : 3;
            cols.push_back(col); roles.push_back(rc); parts.push_back(part); vals.push_back(v);
            if (col > N) N = col;
        }
        fclose(f);
        if (N <= 0) return false;
        role.assign(N+1, -1); partner.assign(N+1, -1);
        sing.assign(N+1, {0,0,0,0}); cond.assign(N+1, {});
        for (size_t i = 0; i < cols.size(); ++i) {
            int c = cols[i]; role[c] = roles[i]; partner[c] = parts[i];
            if (roles[i] == 3) for (int k=0;k<16;++k) cond[c][k] = vals[i][k];
            else               for (int k=0;k<4; ++k) sing[c][k] = vals[i][k];
        }
        return true;
    }
};

// ---------- SDBG helpers ----------
static int firstBase(const SDBG& g, uint64_t n){ static thread_local std::vector<uint8_t> b; b.resize(g.k()); g.GetLabel(n,b.data()); return b[0]; }
static char baseChar(int c1){ return (c1>=1&&c1<=4)?"ACGT"[c1-1]:'N'; }
static int  baseCode(char c){ switch(c){case 'A':return 1;case 'C':return 2;case 'G':return 3;case 'T':return 4;} return 0; }
static std::string labelOf(const SDBG& g, uint64_t node){ int K=g.k(); std::vector<uint8_t> b(K); g.GetLabel(node,b.data()); std::string s; for(int i=0;i<K;++i) s.push_back(baseChar(b[i])); return s; }
static uint64_t lookupKmer(const SDBG& g, const std::string& kmer){
    int K=g.k(); if((int)kmer.size()!=K) return SDBG::kNullID;
    std::vector<uint8_t> seq(K);
    for(int i=0;i<K;++i){ int c=baseCode(kmer[i]); if(!c) return SDBG::kNullID; seq[i]=(uint8_t)c; }
    return g.IndexBinarySearch(seq.data());
}
static const char* stopName(int b0,int b1,int b2){
    if(b0==4&&b1==1&&b2==1) return "TAA"; if(b0==4&&b1==1&&b2==3) return "TAG"; if(b0==4&&b1==3&&b2==1) return "TGA"; return nullptr;
}
static const char* upstreamStop(const SDBG& g, uint64_t v){
    uint64_t i1[8],i2[8],i3[8]; int n1=g.IncomingEdges(v,i1);
    for(int a=0;a<n1;++a){ if(!g.IsValidEdge(i1[a]))continue; int b2=firstBase(g,i1[a]); int n2=g.IncomingEdges(i1[a],i2);
      for(int b=0;b<n2;++b){ if(!g.IsValidEdge(i2[b]))continue; int b1=firstBase(g,i2[b]); int n3=g.IncomingEdges(i2[b],i3);
        for(int c=0;c<n3;++c){ if(!g.IsValidEdge(i3[c]))continue; int b0=firstBase(g,i3[c]); if(const char* s=stopName(b0,b1,b2)) return s; } } }
    return nullptr;
}

// ---------- emit one base (0..3) at column col; update stack + score; false if illegal ----------
static bool emit(const PairModel& M, int col, int b, std::vector<uint8_t>& stk, double& d){
    int r = M.role[col];
    if (r==0 || r==1) { d = M.sing[col][b]; }                       // SL / SR
    else if (r==2)    { d = M.sing[col][b]; stk.push_back((uint8_t)b); }  // P5: marginal + push
    else {                                                          // P3: pop + pair
        if (stk.empty()) return false;
        int x = stk.back(); stk.pop_back();
        d = M.cond[col][x*4 + b];
    }
    return true;
}

// ---------- beam fold-walk from one anchor ----------
struct State { std::string kmer; int col; std::vector<uint8_t> stack; double score; std::string seq; };

static void fold(const SDBG& g, const PairModel& M, uint64_t anchor, std::set<std::string>& cands){
    const size_t WIDTH = 300;
    int K = g.k();
    std::string lab = labelOf(g, anchor);
    // seed: align the anchor's K bases to columns 1..min(K,N) (match-only)
    std::vector<uint8_t> stk; double sc = 0; int col = 1;
    for (int j = 0; j < K && col <= M.N; ++j, ++col){
        int b = baseCode(lab[j]) - 1; if (b < 0) return;
        double d; if (!emit(M, col, b, stk, d)) return; sc += d;
    }
    std::vector<State> beam{ State{lab, col, stk, sc, lab} };
    std::string bestSeq; double bestScore = -1e18;                  // best completion for THIS anchor
    while (!beam.empty()){
        std::vector<State> nxt;
        for (auto& s : beam){
            if (s.col > M.N) { if (s.score > bestScore){ bestScore = s.score; bestSeq = s.seq; } continue; }
            for (char ch : std::string("ACGT")){
                std::string nk = s.kmer.substr(1) + ch;
                uint64_t e = lookupKmer(g, nk); if (e == SDBG::kNullID || !g.IsValidEdge(e)) continue;
                std::vector<uint8_t> stk2 = s.stack; double d;
                if (!emit(M, s.col, baseCode(ch)-1, stk2, d)) continue;
                nxt.push_back(State{nk, s.col+1, std::move(stk2), s.score+d, s.seq+ch});
            }
        }
        if (nxt.size() > WIDTH){                                    // keep top-WIDTH by HIGHEST score
            std::nth_element(nxt.begin(), nxt.begin()+WIDTH, nxt.end(),
                             [](const State&a, const State&b){ return a.score > b.score; });
            nxt.resize(WIDTH);
        }
        beam = std::move(nxt);
    }
    if (!bestSeq.empty()) cands.insert(bestSeq);                    // emit ONE best fold per anchor
}

int main(int argc, char** argv){
    if (argc < 3){ fprintf(stderr, "usage: %s <sdbg_prefix> <pairmodel.txt> [out.fa]\n", argv[0]); return 1; }
    setvbuf(stderr, nullptr, _IONBF, 0);
    SDBG g; g.LoadFromFile(argv[1]);
    PairModel M; if (!M.load(argv[2])){ fprintf(stderr, "cannot load pairmodel %s\n", argv[2]); return 1; }
    fprintf(stderr, "[fold] k=%d edges=%llu model_cols=%d\n", g.k(), (unsigned long long)g.size(), M.N);
    FILE* out = (argc >= 4) ? fopen(argv[3], "w") : stdout;
    std::set<std::string> cands;
    int anchors = 0;
    for (uint64_t i = 0; i < g.size(); ++i){
        if (!g.IsValidEdge(i) || g.EdgeIndegree(i) < 2) continue;
        if (!upstreamStop(g, i)) continue;
        ++anchors;
        if (anchors % 1000 == 0)
            fprintf(stderr, "  [progress] anchors=%d candidates=%d  edge %llu/%llu\n",
                    anchors, (int)cands.size(), (unsigned long long)i, (unsigned long long)g.size());
        fold(g, M, i, cands);
    }
    int idx = 0;
    for (const auto& s : cands) fprintf(out, ">fold%d len=%zu\n%s\n", idx++, s.size(), s.c_str());
    if (out != stdout) fclose(out);
    fprintf(stderr, "[fold] anchors=%d candidates=%d\n", anchors, (int)cands.size());
    return 0;
}
