# vcr_detect — minimal super-integron (VCR) detector

Finds **VCR / gene-cassette arrays** — the *Vibrio cholerae* super-integron —
directly from raw reads, with no assembly step. It is a sibling of MCAAT: same
graph-theoretic idea (a repeat-spacer array is a *multicycle* in a de Bruijn
graph), retargeted from CRISPR repeats to integron recombination sites.

## The idea in one paragraph

A super-integron has the structure

```
... flank ... VCR  cassette₁  VCR  cassette₂  VCR  ...  cassetteₙ  VCR ... flank ...
```

The **VCR** (~120 bp recombination site) is near-**identical** between cassettes.
In a k-mer de Bruijn graph that identical VCR collapses into a single
**high-multiplicity hub** (it is seen once per cassette, so its coverage is
~n× the cassette coverage). Every gene cassette then becomes a distinct
**cycle** that leaves the hub and returns to it. Detecting the array reduces to:
find a high-multiplicity hub and enumerate the cycles through it. That is the
exact CRISPR "multicycle", just at gene-cassette scale.

Why VCRs and not generic integrons: mobile-integron *attC* sites are conserved
in *secondary structure*, not primary sequence, so they would **not** collapse
to one hub. *V. cholerae* VCRs are unusually homogeneous within a genome
(in strain N16961, 149/179 sites are 123 bp at >95% identity), which is exactly
what makes them collapse cleanly — the right first target.

## Build

```bash
cd vcr_detector
cmake -S . -B build && cmake --build build
# binary: build/vcr_detect   (or build/vcr_detect.exe on Windows)
```

Or directly: `g++ -O2 -std=c++17 vcr_detect.cpp -o vcr_detect`.
No external dependencies — just a C++17 compiler.

## Quick test (generate synthetic data, then detect)

```bash
./vcr_detect gen -o test.fa        # writes reads; prints ground truth to stderr
./vcr_detect detect test.fa        # should report 1 array, 6 cassettes
```

`gen` builds a synthetic super-integron (random VCR + 6 random cassettes),
shreds it into overlapping reads, and prints the ground-truth VCR to stderr so
you can check the recovered `REPEAT` against it. **Keep test inputs small** —
this v0 enumerates cycles by bounded DFS and is meant for small graphs only.

> Note: the synthetic VCR is 80 bp (so the demo reports `repeat_len=80`), kept
> short just to keep the test graph small. Real *V. cholerae* VCRs are ~120 bp;
> nothing in the algorithm depends on the exact length.

## Usage

```
vcr_detect detect <reads.fa|.fq> [-k N] [-c MINCOPIES] [-n MINCASSETTES] [-o OUT]
vcr_detect gen    [-o OUT] [-s SEED] [--cassettes N]
```

Input may be FASTA or FASTQ, plain text (not gzipped in v0).

## Parameters — only three

| Flag | Default | Meaning |
|------|---------|---------|
| `-k` | 23 | k-mer size. Must be < the VCR length so the VCR shares k-mers across copies. |
| `-c` | 20 | **min-copies**: multiplicity a node needs to count as repeat/hub. Set it **above the per-cassette read depth** and **below depth × cassette-count**. This is the one parameter that matters. |
| `-n` | 2 | **min-cassettes**: how many distinct cassettes a hub needs before it is reported as an array. |

That is the whole knob set. Everything else (cassette length, repeat length,
coverage ratios, step caps) is either derived automatically or a fixed internal
safety bound — deliberately *not* exposed.

Choosing `-c`: if your reads cover the region at depth *d* and the array has *n*
cassettes, cassette k-mers sit near *d* and VCR k-mers near *n·d*. Any `-c`
between those works; `d·2` is a safe default-ish choice.

## Output

```
# vcr_detect v0 — super-integron / VCR array output
# graph_nodes=... k=23 min_copies=20 min_cassettes=2
# arrays=1

>Array_1  repeat_len=80  cassettes=6
REPEAT	<consensus VCR sequence>
  cassette_1	len=263	<sequence>
  cassette_2	len=301	<sequence>
  ...
```

`REPEAT` is the recovered VCR (the high-multiplicity hub); each `cassette_i` is
one cycle through it. Cassette boundaries are approximate to within ~k bases.

The repeat is reported as the **maximal sequence shared by all cassettes** — by
sequence alone the tool cannot tell the VCR apart from any sequence that always
sits next to it. Two cases:

- **Coincidental** sharing of a random junction base or two — vanishingly
  unlikely as cassette count grows (~4^−(n−1) for n cassettes). Adds 1–2 bases.
- **Genuinely conserved** sequence adjacent to the VCR (e.g. a conserved
  attC-proximal motif present in every cassette) — here the sharing probability
  is 1 and the reported repeat absorbs the **entire** conserved stretch,
  independent of cassette count. The true VCR is still contained in the reported
  repeat; it is just prefixed by the conserved flank.

A future version can cap this using the multiplicity profile the graph already
has — the true VCR sits at the full hub-coverage plateau (seen once per VCR
copy), while a shared cassette flank sits one coverage-step lower — but v0
deliberately keeps the threshold-free common-suffix rule.

## How it works (pipeline)

1. **Build** a plain in-memory de Bruijn graph: nodes are k-mers with a
   multiplicity (occurrence count); edges are observed (k+1)-mers.
2. **Anchor**: pick high-multiplicity nodes (`mult ≥ min_copies`) with
   out-degree ≥ 2 — the points where a collapsed VCR fans out to its cassettes.
3. **Enumerate cycles** from each anchor by bounded depth-first search
   (simple cycles back to the anchor).
4. **Separate** repeat from cassette by the **longest common suffix** shared
   across all of a hub's cycles, compared node-by-node (not by coverage): that
   shared tail plus the hub base is the VCR; the variable middle of each cycle is
   its cassette. Multiplicity is used only for anchor selection in step 2, not
   here.
5. **Report** a hub with ≥ `min_cassettes` distinct cassettes as an array; hubs
   are processed highest-multiplicity first and their nodes are "claimed" so a
   frayed sub-hub cannot re-report a fragment of the same array.

## Scope & limitations of this v0 (by design)

- **Small graphs only.** Cycle enumeration is bounded DFS; on a real metagenome
  the cassette interiors (~kb) make cycles thousands of nodes long and shared
  subsequences cause branching, which is expensive. Do not run it on large
  inputs yet.
- **Forward strand only.** Real reads come from both strands; v0 does not yet
  build canonical k-mers.
- **Plain hash-map graph.** Memory scales with distinct k-mers — fine for tests,
  not for metagenomic scale.
- **Identical VCRs merge.** Two distinct loci that share an identical VCR are
  reported as one array (the union of their cassettes) — the de Bruijn graph has
  no information to separate them without paired or long reads. For real
  *V. cholerae* this is mostly moot (a genome carries one super-integron), but
  near-identical VCRs across co-occurring strains will pool their cassettes.

## Roadmap (the scaling fixes, in order)

1. **Unitig compaction first.** Compress every non-branching chain to one edge
   *before* enumeration. A whole cassette interior becomes a single edge and a
   cassette cycle becomes a 2–3-hop cycle, so depth-bounded enumeration stops
   exploding. This is the single most important change for real data. (MEGAHIT's
   `unitig_graph.h` already does this and could be reused.)
2. **BFS reachability instead of full enumeration** — test "does a hub→hub cycle
   exist?" and reconstruct the cassette by a unitig walk, rather than listing
   every interior path.
3. **Disk-backed search frontier.** If branching still blows up memory, spill the
   DFS/BFS stack (the explicit work-list) to disk instead of recursing in RAM —
   an external-memory traversal. (Planned; the current DFS is in-RAM recursion.)
4. **Strand-aware (canonical) k-mers** so both-strand reads merge correctly.
5. **Succinct graph for scale.** Swap the hash-map graph for MEGAHIT's succinct
   de Bruijn graph (as MCAAT uses) once the algorithm is settled.
6. **attC fold check** to separate true VCRs from tandem-gene look-alikes
   (the VCR is an imperfect inverted repeat — require palindromic self-similarity
   of the recovered repeat).

## Relationship to MCAAT

Same multicycle principle as the parent CRISPR tool, with the parameters and
length scales retuned for integron cassettes and the heavy graph machinery
stripped out for a minimal, self-contained starting point.
