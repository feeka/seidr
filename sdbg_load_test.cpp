// Minimal proof that we can compile against MEGAHIT's SDBG and load a graph
// built by `megahit_core read2sdbg`. Nothing more — just open it and report.
#include <sdbg/sdbg.h>
#include <cstdio>

int main(int argc, char** argv) {
    if (argc < 2) { fprintf(stderr, "usage: %s <sdbg_prefix>\n", argv[0]); return 1; }
    SDBG dbg;
    dbg.LoadFromFile(argv[1]);
    printf("LOADED SDBG: k=%u  edges=%llu\n", dbg.k(), (unsigned long long)dbg.size());
    // sanity: count valid edges and how many have indegree >= 2
    unsigned long long valid = 0, indeg2 = 0;
    for (uint64_t i = 0; i < dbg.size(); ++i) {
        if (!dbg.IsValidEdge(i)) continue;
        ++valid;
        if (dbg.EdgeIndegree(i) >= 2) ++indeg2;
    }
    printf("valid_edges=%llu  indeg>=2=%llu\n", valid, indeg2);
    return 0;
}
