// vcr_detect.cpp — minimal VCR / super-integron multicycle detector
// -----------------------------------------------------------------------------
// Detects repeat-cassette arrays (the Vibrio cholerae super-integron signature)
// directly from raw reads, with no assembly step.
//
// Biology in one line: a super-integron is  VCR-cassette-VCR-cassette-...-VCR ,
// where the ~120 bp VCR recombination site is near-IDENTICAL between cassettes.
// In a de Bruijn graph that identical VCR collapses into one high-multiplicity
// "hub" path, and every gene cassette becomes a distinct cycle leaving and
// returning to that hub — a CRISPR-style multicycle, just at cassette scale.
// This tool finds exactly that motif.
//
// This is the MINIMAL v0:
//   * plain in-memory k-mer de Bruijn graph (no succinct/MEGAHIT graph),
//   * forward strand only,
//   * bounded DFS cycle enumeration (fine on small graphs; see README roadmap).
//
// Subcommands
//   vcr_detect detect <reads.fa|.fq> [-k N] [-c MINCOPIES] [-n MINCASSETTES] [-o OUT]
//   vcr_detect gen    [-o OUT] [-s SEED] [--cassettes N]     # synthetic test data
//
// (If the first argument is a file rather than "detect"/"gen", detect is assumed.)
// -----------------------------------------------------------------------------

#include <cstdint>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

using std::string;
using std::vector;
using std::unordered_map;
using std::unordered_set;

// ----------------------------- DNA helpers -----------------------------------
static inline bool isACGT(char c) { return c == 'A' || c == 'C' || c == 'G' || c == 'T'; }

static bool allACGT(const string& s) {
    for (char c : s) if (!isACGT(c)) return false;
    return true;
}

static char comp(char c) {
    switch (c) { case 'A': return 'T'; case 'C': return 'G';
                 case 'G': return 'C'; case 'T': return 'A'; default: return 'N'; }
}

[[maybe_unused]] static string revcomp(const string& s) {   // reserved for strand-aware mode
    string r(s.size(), 'N');
    for (size_t i = 0; i < s.size(); ++i) r[i] = comp(s[s.size() - 1 - i]);
    return r;
}

// Parse a numeric CLI argument, exiting cleanly (not aborting) on garbage input.
static long parse_num(const string& s, const char* flag) {
    try {
        size_t pos = 0;
        long v = std::stol(s, &pos);
        if (pos != s.size()) throw std::invalid_argument("trailing characters");
        return v;
    } catch (const std::exception&) {
        std::cerr << "error: invalid value '" << s << "' for " << flag << "\n";
        std::exit(1);
    }
}

// ----------------------------- read input ------------------------------------
// Minimal FASTA/FASTQ reader. Auto-detects by the first non-blank character.
static vector<string> read_sequences(const string& path) {
    std::ifstream in(path);
    if (!in.is_open()) { std::cerr << "error: cannot open " << path << "\n"; std::exit(1); }
    vector<string> seqs;
    string line, cur;
    char fmt = 0;                          // '>' FASTA, '@' FASTQ
    int fq_line = 0;                       // position within a 4-line FASTQ record
    while (std::getline(in, line)) {
        if (!line.empty() && (line.back() == '\r')) line.pop_back();
        if (line.empty()) continue;
        if (fmt == 0) fmt = line[0];
        if (fmt == '>') {
            if (line[0] == '>') { if (!cur.empty()) { seqs.push_back(cur); cur.clear(); } }
            else cur += line;
        } else {                           // FASTQ: sequence is line 2 of every 4
            int p = fq_line % 4;
            if (p == 1) seqs.push_back(line);
            ++fq_line;
        }
    }
    if (fmt == '>' && !cur.empty()) seqs.push_back(cur);
    for (auto& s : seqs) for (char& c : s) c = (char)std::toupper((unsigned char)c);
    return seqs;
}

// ----------------------------- de Bruijn graph -------------------------------
// Node  = k-mer.            node[kmer]  = number of times the k-mer was observed
//                                          (its multiplicity / coverage).
// Edge  = (k+1)-mer present in the data. Adjacency u->v exists iff the (k+1)-mer
//         u + last_base_of(v) was observed.
struct DBG {
    int k = 23;
    unordered_map<string, uint32_t> node;   // k-mer  -> multiplicity
    unordered_set<string>           edge;    // (k+1)-mers that exist

    void build(const vector<string>& reads, int kk) {
        k = kk;
        const size_t k1 = (size_t)k + 1;
        for (const string& r : reads) {
            if (r.size() >= (size_t)k)
                for (size_t i = 0; i + k <= r.size(); ++i) {
                    string km = r.substr(i, k);
                    if (allACGT(km)) ++node[km];
                }
            if (r.size() >= k1)
                for (size_t i = 0; i + k1 <= r.size(); ++i) {
                    string e = r.substr(i, k1);
                    if (allACGT(e)) edge.insert(e);
                }
        }
    }

    uint32_t mult(const string& u) const {
        auto it = node.find(u);
        return it == node.end() ? 0u : it->second;
    }

    vector<string> succ(const string& u) const {
        vector<string> out;
        string e = u; e.push_back('X');
        for (char b : {'A', 'C', 'G', 'T'}) {
            e.back() = b;
            if (edge.count(e)) out.push_back(e.substr(1));   // drop first base
        }
        return out;
    }

    int outdeg(const string& u) const { return (int)succ(u).size(); }
};

// ----------------------------- cycle utilities -------------------------------
// Decode a list of consecutive nodes (each overlaps the previous by k-1) into a
// single sequence: full first k-mer, then one trailing base per subsequent node.
static string decode(const vector<string>& nodes) {
    if (nodes.empty()) return "";
    string s = nodes[0];
    for (size_t i = 1; i < nodes.size(); ++i) s += nodes[i].back();
    return s;
}

// Every cycle through a repeat hub has the form  [hub, cassette_i ..., VCR body].
// The VCR body is identical for all cassettes (they all run into the same next
// copy of the repeat), so it is the LONGEST COMMON SUFFIX of the cycles — a
// threshold-free way to separate repeat from cassette, robust to coverage and
// to bases shared by only some cassettes. Returns 0 if there is no shared tail.
static int common_suffix_len(const vector<vector<string>>& cyc) {
    if (cyc.size() < 2) return 0;
    size_t mn = (size_t)-1;
    for (const auto& c : cyc) mn = std::min(mn, c.size());
    int maxL = (int)mn - 2;                   // leave the hub (front) + >=1 cassette node
    int L = 0;
    for (int t = 1; t <= maxL; ++t) {
        const string& ref = cyc[0][cyc[0].size() - t];
        bool ok = true;
        for (size_t i = 1; i < cyc.size(); ++i)
            if (cyc[i][cyc[i].size() - t] != ref) { ok = false; break; }
        if (!ok) break;
        L = t;
    }
    return L;
}

// ----------------------------- detection -------------------------------------
// Enumerate simple cycles from a hub, collecting each as a node list
// [hub, ... , last-before-closing]. The search is an EXPLICIT-STACK depth-first
// walk (not recursion): depth is bounded by the heap, so a long looping graph
// cannot overflow the call stack, and the work-stack is exactly the structure a
// future version could spill to disk. Safety caps are constants, not user knobs.
struct Detector {
    const DBG& g;
    static constexpr long MAX_STEPS  = 5'000'000;
    static constexpr int  MAX_CYCLES = 4096;
    static constexpr int  MAX_NODES  = 8000;
    static constexpr int  MIN_NODES  = 3;

    string anchor;
    vector<vector<string>> found;

    explicit Detector(const DBG& g_) : g(g_) {}

    struct Frame { string node; vector<string> succ; size_t idx; };

    void run_from(const string& a) {
        anchor = a;
        found.clear();
        long steps = 0;
        int  cycles = 0;

        vector<Frame> stk;                         // the explicit DFS stack
        vector<string> path;                       // node list = current stack path
        unordered_set<string> inpath;              // membership, to keep cycles simple

        auto push = [&](const string& u) {
            stk.push_back({u, g.succ(u), 0});
            path.push_back(u);
            inpath.insert(u);
        };
        auto pop = [&]() {
            inpath.erase(stk.back().node);
            path.pop_back();
            stk.pop_back();
        };

        push(a);
        while (!stk.empty()) {
            if (steps > MAX_STEPS || cycles >= MAX_CYCLES) break;
            Frame& f = stk.back();
            if (f.idx >= f.succ.size()) { pop(); continue; }
            string v = f.succ[f.idx++];            // copy: push() may reallocate stk
            if (++steps > MAX_STEPS) break;
            if (v == anchor) {
                if ((int)path.size() >= MIN_NODES) found.push_back(path);
                if (++cycles >= MAX_CYCLES) break;
                continue;                          // closing edge: do not descend
            }
            if (inpath.count(v)) continue;         // already on the path
            if ((int)path.size() >= MAX_NODES) continue;   // depth cap: do not descend
            push(v);
        }
    }
};

struct Array { string repeat; vector<string> cassettes; };

static vector<Array> detect(const DBG& g, uint32_t min_copies, int min_cassettes) {
    // Anchors: high-multiplicity branch points — the end of a collapsed repeat
    // where it fans out to its cassettes (mirrors MCAAT's ChunkStartNodes).
    // Process them highest-multiplicity first so the true peak-coverage repeat
    // hub is found and "claimed" before any lower-coverage frayed sub-hub can
    // re-report a fragment of the same array.
    vector<std::pair<uint32_t, string>> anchors;
    for (const auto& kv : g.node)
        if (kv.second >= min_copies && g.outdeg(kv.first) >= 2)
            anchors.emplace_back(kv.second, kv.first);
    std::sort(anchors.begin(), anchors.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });

    vector<Array> out;
    unordered_set<string> claimed;
    Detector det(g);

    for (const auto& pr : anchors) {
        const string& a = pr.second;
        if (claimed.count(a)) continue;
        det.run_from(a);
        if ((int)det.found.size() < min_cassettes) continue;

        int L = common_suffix_len(det.found);
        if (L < 1) continue;                         // cycles share no repeat body

        // repeat = shared VCR body (the common suffix) + the hub's own base
        const vector<string>& ref = det.found[0];
        vector<string> suffix(ref.end() - L, ref.end());
        string repeat = decode(suffix) + a.back();

        unordered_set<string> cas;
        for (const auto& cyc : det.found) {
            vector<string> mid(cyc.begin() + 1, cyc.end() - L);
            if (!mid.empty()) cas.insert(decode(mid));
        }
        if ((int)cas.size() < min_cassettes) continue;

        Array A; A.repeat = repeat;
        A.cassettes.assign(cas.begin(), cas.end());
        std::sort(A.cassettes.begin(), A.cassettes.end(),
                  [](const string& x, const string& y) {
                      return x.size() != y.size() ? x.size() < y.size() : x < y; });
        out.push_back(std::move(A));
        for (const auto& cyc : det.found)            // claim nodes against re-report
            for (const auto& n : cyc) claimed.insert(n);
    }

    std::sort(out.begin(), out.end(), [](const Array& a, const Array& b) {
        return a.cassettes.size() > b.cassettes.size();
    });
    return out;
}

// ----------------------------- output ----------------------------------------
static void write_arrays(std::ostream& os, const vector<Array>& arrays,
                         const DBG& g, uint32_t min_copies, int min_cassettes) {
    os << "# vcr_detect v0 — super-integron / VCR array output\n";
    os << "# graph_nodes=" << g.node.size() << " k=" << g.k
       << " min_copies=" << min_copies << " min_cassettes=" << min_cassettes << "\n";
    os << "# arrays=" << arrays.size() << "\n\n";
    int n = 0;
    for (const Array& A : arrays) {
        os << ">Array_" << (++n) << "  repeat_len=" << A.repeat.size()
           << "  cassettes=" << A.cassettes.size() << "\n";
        os << "REPEAT\t" << A.repeat << "\n";
        int c = 0;
        for (const string& cas : A.cassettes)
            os << "  cassette_" << (++c) << "\tlen=" << cas.size() << "\t" << cas << "\n";
        os << "\n";
    }
}

// ----------------------------- subcommands -----------------------------------
static void usage() {
    std::cerr <<
      "vcr_detect — minimal VCR / super-integron multicycle detector\n\n"
      "  vcr_detect detect <reads.fa|.fq> [-k N] [-c MINCOPIES] [-n MINCASSETTES] [-o OUT]\n"
      "  vcr_detect gen    [-o OUT] [-s SEED] [--cassettes N]\n\n"
      "detect parameters (only three real knobs):\n"
      "  -k  k-mer size                 (default 23)\n"
      "  -c  min-copies   hub multiplicity threshold for the repeat (default 20)\n"
      "  -n  min-cassettes              distinct cassettes to call an array (default 2)\n";
}

static int cmd_detect(const vector<string>& a) {
    string input, out;
    int k = 23, min_cassettes = 2;
    long min_copies = 20;
    for (size_t i = 0; i < a.size(); ++i) {
        const string& t = a[i];
        if      (t == "-k" && i + 1 < a.size()) k = (int)parse_num(a[++i], "-k");
        else if (t == "-c" && i + 1 < a.size()) min_copies = parse_num(a[++i], "-c");
        else if (t == "-n" && i + 1 < a.size()) min_cassettes = (int)parse_num(a[++i], "-n");
        else if (t == "-o" && i + 1 < a.size()) out = a[++i];
        else if (!t.empty() && t[0] != '-' && input.empty()) input = t;
    }
    if (input.empty()) { usage(); return 1; }
    if (k < 4) { std::cerr << "error: k too small\n"; return 1; }

    vector<string> reads = read_sequences(input);
    DBG g; g.build(reads, k);
    std::cerr << "[vcr_detect] reads=" << reads.size()
              << " nodes=" << g.node.size() << " edges=" << g.edge.size() << "\n";

    vector<Array> arrays = detect(g, (uint32_t)min_copies, min_cassettes);

    if (out.empty()) {
        write_arrays(std::cout, arrays, g, (uint32_t)min_copies, min_cassettes);
    } else {
        std::ofstream os(out);
        write_arrays(os, arrays, g, (uint32_t)min_copies, min_cassettes);
        std::cerr << "[vcr_detect] wrote " << arrays.size() << " array(s) to " << out << "\n";
    }
    return 0;
}

static string rand_seq(int n) {
    static const char B[4] = {'A', 'C', 'G', 'T'};
    string s(n, 'A');
    for (int i = 0; i < n; ++i) s[i] = B[std::rand() & 3];
    return s;
}

// Synthetic super-integron: flank + (VCR + cassette) x N + VCR + flank, then
// shred into overlapping reads. Ground truth is printed to stderr.
static int cmd_gen(const vector<string>& a) {
    string out;
    unsigned seed = 42;
    int nCass = 6;
    const int vcrLen = 80, flank = 100, readLen = 120, stride = 12;
    const int cassMin = 250, cassMax = 450;
    for (size_t i = 0; i < a.size(); ++i) {
        const string& t = a[i];
        if      (t == "-o" && i + 1 < a.size()) out = a[++i];
        else if (t == "-s" && i + 1 < a.size()) seed = (unsigned)parse_num(a[++i], "-s");
        else if (t == "--cassettes" && i + 1 < a.size()) nCass = (int)parse_num(a[++i], "--cassettes");
    }
    std::srand(seed);

    string vcr = rand_seq(vcrLen);
    vector<string> cass;
    for (int i = 0; i < nCass; ++i)
        cass.push_back(rand_seq(cassMin + (std::rand() % (cassMax - cassMin + 1))));

    string arr = rand_seq(flank);
    for (int i = 0; i < nCass; ++i) { arr += vcr; arr += cass[i]; }
    arr += vcr; arr += rand_seq(flank);

    std::ostream* osp = &std::cout;
    std::ofstream ofs;
    if (!out.empty()) { ofs.open(out); osp = &ofs; }
    std::ostream& os = *osp;

    int id = 0;
    size_t p = 0;
    for (; p + readLen <= arr.size(); p += stride)
        os << ">read_" << id++ << "\n" << arr.substr(p, readLen) << "\n";
    if (arr.size() >= (size_t)readLen)                 // make sure the tail is covered
        os << ">read_" << id++ << "\n" << arr.substr(arr.size() - readLen) << "\n";

    std::cerr << "[gen] GROUND TRUTH: vcr_len=" << vcrLen
              << " cassettes=" << nCass << " array_len=" << arr.size()
              << " reads=" << id << "\n";
    std::cerr << "[gen] vcr=" << vcr << "\n";
    return 0;
}

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }
    vector<string> rest;
    string cmd = argv[1];
    if (cmd == "gen") { for (int i = 2; i < argc; ++i) rest.push_back(argv[i]); return cmd_gen(rest); }
    if (cmd == "detect") { for (int i = 2; i < argc; ++i) rest.push_back(argv[i]); return cmd_detect(rest); }
    for (int i = 1; i < argc; ++i) rest.push_back(argv[i]);   // implicit detect
    return cmd_detect(rest);
}
