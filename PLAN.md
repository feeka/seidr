# VCR detection on the succinct de Bruijn graph — implementation plan

This is the method **as specified by the project lead**, written for someone in
the field to implement directly. Each step states what to do and, where the spec
left a detail open, says so explicitly instead of inventing a choice.

Legend:
- **[SPEC]** — exactly what was specified. Implement as written.
- **[ASSUMPTION]** — my reading of a detail that was *not* explicitly specified.
  Treat as provisional; confirm before relying on it.
- **[OPEN]** — a decision still needed (listed again at the end).

---

## 0. What we are detecting (background, factual)

- A **VCR** (Vibrio Cholerae Repeat) is the ~120–126 bp **attC recombination
  site** of the *V. cholerae* chromosomal super-integron. (Verified earlier:
  N16961 majority class is 123 bp.)
- The super-integron is an array: `… VCR — cassette — VCR — cassette — … VCR …`,
  i.e. one VCR between **every** gene cassette. Each **cassette = one ORF**
  (start codon `ATG` … stop codon `TAA/TAG/TGA`), with the VCR sitting at the
  cassette's 3′ end. So along the sequence the local pattern at each junction is:
  `… ORF_i [STOP]  VCR  [START=ATG] ORF_{i+1} …`.
- **Goal:** detect the VCRs (and, by extension, the array) **directly from raw
  reads on the de Bruijn graph, with no assembly step.**

### 0.1 The VCR forms a cycle (the foundation)

The VCR is the **same sequence between every cassette**, so in the de Bruijn
graph all copies of it **collapse onto one path**. The array

```
… [STOP] VCR [START=ATG] ORF_1 [STOP] VCR [START] ORF_2 [STOP] VCR …
```

therefore forms a **cycle**: leave the VCR, go through one cassette
`[START]ORF[STOP]`, and arrive **back at the same collapsed VCR**. Different
cassettes are different loops off that one VCR → a **multicycle with the VCR as
the hub** (the MCAAT shape: repeat = collapsed hub, spacers/cassettes = loops).
This is why the VCR start is a convergence point: every cassette's STOP feeds
back into the shared VCR, giving the VCR-start node **in-degree ≥ 2**.

---

## 1. Input

1. Raw sequencing reads (FASTA or FASTQ).
2. A **VCR profile** — the scoring model for the traversal (see Step 4).

---

## 2. Build the graph: use MEGAHIT's succinct de Bruijn graph (SDBG)

**[SPEC]** Use **MEGAHIT's SDBG**, the same graph MCAAT uses. Do **not** build a
custom/ad-hoc graph.

1. Build the SDBG from the reads with MEGAHIT's core, the two-stage call MCAAT
   uses:
   - `megahit_core buildlib <data.lib> <prefix>` (writes the read library), then
   - `megahit_core read2sdbg --host_mem <bytes> --read_lib_file <prefix> -m 1 -k 23 --num_cpu_threads <n> --o <graph_prefix>`
   - **[fact]** k = 23 and `-m 1` are what MCAAT uses (`-m 1` = keep every edge,
     no coverage pruning).
2. Load the resulting SDBG and access it through its API:
   - `OutgoingEdges(node)` / `IncomingEdges(node)` — successors / predecessors;
   - `EdgeIndegree(node)` / `EdgeOutdegree(node)` — in/out degree;
   - `EdgeMultiplicity(node)` — coverage of that (k+1)-mer;
   - `GetLabel(node)` — recover the nucleotide(s).
3. A "node" here is an SDBG edge/(k+1)-mer index (`uint64_t`), as in MCAAT.

**[OPEN]** Confirm k (23?) and whether `-m 1` is what we want for VCR detection.

---

## 3. Find the anchor nodes (the VCR start)

**[SPEC]** Anchor = a node with **indegree ≥ 2** whose **3 PREVIOUS (upstream)
nodes spell a STOP codon**.

Rationale: at a VCR start, the cassette ORFs that precede different copies
converge into the one shared VCR, so that node has in-degree ≥ 2; and the ORF
immediately before ends in a stop codon, which sits in the 3 nodes upstream.

**Exact SDBG geometry (read from `sdbg.h`, not assumed):**
- `GetLabel` (`sdbg.h:214`) reconstructs a k-mer by walking **`Backward`** and
  reading `GetW`, writing into a **decreasing** index: `seq[k-1]` first (1 step),
  `seq[0]` last (k steps). So **`seq[0]` is the 5′-most base** and going backward
  fills lower indices.
- Bases are **1-based**: `GetW` returns 1..4 (strip the "last" flag, `-=`
  `kAlphabetSize`), decoded `"ACGT"[seq[i]-1]` → **A=1, C=2, G=3, T=4** (0 = `$`).

Procedure, per node `v`:
1. If `EdgeIndegree(v) < 2`: skip.
2. Read the **3 bases immediately upstream** of `v` by walking `IncomingEdges`
   3 times. Mirroring `GetLabel`'s backward-decreasing fill:
   - 1st `IncomingEdges` → base immediately upstream → `seq[2]` (3′ of the codon)
   - 2nd `IncomingEdges` → `seq[1]`
   - 3rd `IncomingEdges` → furthest base → `seq[0]` (5′ of the codon)

   So the codon **read 5′→3′ is `seq[0] seq[1] seq[2]`** (deepest step = `seq[0]`).
   **Orientation matters:** reading in walk-order without this reversal tests the
   reverse codon and misses every stop.
3. Accept `v` as an anchor if `seq[0..2] ∈ { TAA, TAG, TGA }`, which in the
   1-based encoding are **`(4,1,1)`, `(4,1,3)`, `(4,3,1)`**.
   - **[ASSUMPTION, to confirm]** because in-degree ≥ 2 there are several incoming
     paths — accept if **any** incoming 3-step path spells a stop (every cassette
     ORF ends in one). Confirm "any" vs. a specific path.

Output of this step: the list of anchor nodes (VCR starts).

---

## 4. The profile (the traversal's weight)

**[SPEC]** The traversal is guided by a **profile** of the VCR (the model). The
**step geometry depends on the profile alphabet:**
- **bp (nucleotide) profile** → advance **1 base** per profile position;
- **aa (amino-acid) profile** → advance **1 codon (3 bases), translated to one
  amino acid,** per profile position.

For the VCR itself (a DNA recombination site) the profile is **nucleotide (bp)**.
(An amino-acid profile is the path for detecting cassette **ORFs/genes**, not the
VCR — relevant when we extend to gene detection.)

**[OPEN — needs your decision] What exactly is the VCR profile, and how is a
single step scored?** Candidates we discussed:
- (a) a nucleotide position-specific profile (e.g. log-odds per base per position)
  built from real VCR instances — e.g. emitted from `Vibrionales.cm` via
  `cmemit`; or
- (b) the covariance model itself.
I will not pick one of these on my own. Tell me which, and the per-step score
definition (e.g. log-odds of the next base at the current profile position).

---

## 5. Traverse to find the VCR — simple greedy weighted walk

**[SPEC]** From each anchor, extend the path by a **greedy** weighted traversal:
**at each step take the single highest-scoring next step and travel along it.**
No beam search, no enumeration of many paths.

Procedure, from an anchor `v`:
1. Maintain: current node, current profile position, accumulated path of bases.
2. At each step, look at the candidate next steps:
   - **bp profile:** each outgoing edge adds one base `b`; score `b` against the
     current profile position.
   - **aa profile:** read the next **codon** (3 bases) along outgoing edges,
     translate to an amino acid, score it against the current profile position.
3. **Pick the single successor with the highest score.** Move to it. Advance the
   profile position by one (bp: +1 base; aa: +1 codon).
4. Continue extending until you reach the **next ORF's START codon (`ATG`)**,
   which marks the VCR's end.
5. The walked path between the STOP (upstream of the anchor, Step 3) and this
   START is the VCR. Pattern recovered: **STOP — VCR — START.**

**[OPEN]** Termination: stop at the next `ATG`, or when the profile reaches its
end, or both? (Step 3 says STOP→VCR→START, i.e. stop at `ATG`; confirm.)
**[OPEN]** Greedy as described advances one profile position per step (no
insert/delete). Real VCRs vary in length by a few bases. Do we keep it strict
(one-position-per-step), or allow indel moves? — your call; I will not add indel
handling unless you ask.

---

## 6. Output

For each anchor that produced a walk to a START: emit the recovered VCR sequence
(and, optionally, its anchor/coordinates and the bounding STOP/START).

---

## 7. Validate (separate from detection)

**[SPEC]** Use `cmsearch` with the real attC covariance model
(`models/Vibrionales.cm`) **only as a validator** of recovered VCRs — never as
the detector.

- `cmsearch --tblout hits.tbl models/Vibrionales.cm recovered.fa`
- **[honest caveat]** If the profile in Step 4 is itself derived from this same
  model, a `cmsearch` hit is partly circular and is weak evidence. On synthetic
  tests, the real check is identity of recovered VCRs to the **planted** VCRs
  (Step 8.5), not cmsearch alone.

---

## 8. Test plan (your design)

1. **Generate VCRs** from the model — best consensus / instances (e.g. `cmemit`
   from `Vibrionales.cm`; convert RNA `U`→`T` for DNA).
2. **Build ORFs** with `ATG` … `STOP` in correct reading-frame geometry (no
   premature in-frame stop).
3. **Build the array:** `random_flank + (ORF + VCR) × N + random_flank`.
4. **Make reads** from the array; **build the SDBG** (Step 2); run detection
   (Steps 3–6).
5. **Check against ground truth:** align each recovered VCR to the **planted**
   VCR copies and report % identity. Report the failing numbers first. Also check
   whether the cassette ORFs are recoverable.

---

## 9. Cassette/ORF crossing — bounded beam traversal between VCRs (LOCKED SPEC)

Confirmed by the project lead (2026-06-22). This is the **gene step**, run after a
VCR is detected. **Problem class: bounded START→STOP beam traversal — NOT cycle
enumeration.** (The VCR collapses to one high-multiplicity hub; enumerating cycles
through it explodes combinatorially and most cycles are spurious repeat-induced
pairings. A bounded *local* beam from each START is feasible: cost ≈ width × depth,
independent of the whole-graph size.)

**9.1 Array bound.** A cassette array has **2–300 genes**. Accept only arrays whose
recovered cassette count is in **[2, 300]**; reject otherwise. (Literature: mobile
integrons 1–8 cassettes; Vibrio super-integrons up to ~300; N16961 ≈ 166 ORFs /
179 cassettes.)

**9.2 START.** After a detected VCR, find the cassette START codon **`ATG` within
±50 nodes** downstream of the VCR end. That `ATG` fixes the reading frame. The
recovered ORF includes the `ATG`.

**9.3 ORF crossing.** From the START, walk forward:
- length window **[210, 1200] bp** (`ATG`→`STOP`; the ORF only, not ORF+attC);
- a path is **valid only if an in-frame STOP is reached inside the window**:
  STOP before 210 bp → **discard**; reach 1200 bp with no STOP → **discard**.
- (Literature bounds: 210 bp ≈ the 70-aa "real ORF" floor below which cassettes
  are gene-less; 1200 bp ≈ the largest N16961 cassette ORF, VCA0308 = 1209 bp.)

**9.4 Beam.** Keep the **top 2000** partial paths, ranked by
**moving-average multiplicity-consistency + stop legality**:
- *multiplicity consistency:* prefer the successor whose `EdgeMultiplicity` stays
  closest to the path's **running average** multiplicity (penalty = `|mult − avg|`,
  accumulated; lower is better). Rationale: uniform depth ⇒ one true contiguous
  cassette path holds ~constant coverage; a jump signals leaving the cassette.
- *stop legality:* a path placing an in-frame STOP **before** 210 bp is illegal
  (discard); a path is **complete** when it places one in **[210, 1200]**.
- **[implementation reading — flagged, per rule 2]** multiplicity is scored as the
  negative running deviation; stop legality is a **hard gate**, not a soft score.
  Correct me if you meant a soft stop term in the ranking.

**9.5 Cycle rule.** If a path revisits an already-visited node `N`
(e.g. `1→2→3→4→5→3`), it may leave `N` by a *different* edge only if
**outdegree(`N`) ≥ 2**; otherwise **discard the path**.

**9.6 Geometry.** Codon/aa: 3 bp per step, in the frame fixed by 9.2.

### 9.7 ORF test plan (against planted ground truth)

1. Plant **distinct** ORFs (each `ATG` … in-frame … `STOP`, length in [210,1200])
   between **identical** VCRs: `flank + VCR + (ORF_i + VCR)×N + flank`.
2. Write every planted ORF sequence to `truth.txt` (`orf<i>\t<seq>`).
3. Build the SDBG, detect VCRs, run 9.2–9.5 from each VCR.
4. **Check vs ground truth:** the recovered ORF set must equal the planted ORF set
   **byte-for-byte**. Report the actual count recovered / planted and any mismatch
   first, before any success claim. This is the divergent/distinct-cassette case —
   the real test of whether the beam stays on one cassette at the VCR fan-out.

---

## Decisions I need from you before coding (so I don't substitute my own choices)

1. **Profile (Step 4):** nucleotide profile built from VCR instances, or the
   covariance model itself? And the exact per-step score formula.
2. **Anchor STOP check (Step 3):** "3 previous nodes" = 3 bases immediately
   upstream — correct? Accept if *any* upstream path spells a stop, or a specific
   one?
3. **Termination (Step 5):** stop at the next `ATG`, at profile end, or both?
4. **Indels (Step 5):** strict one-profile-position-per-step, or allow
   insert/delete moves for VCR length variation?
5. **k-mer size / SDBG params (Step 2):** k = 23 and `-m 1` as MCAAT, or other?
6. **Scope now:** VCR only (bp profile), or also cassette ORFs (aa profile)?

---

## Honest status (facts, not claims)

**Implemented:** VCR detection (`vcr_traverse.cpp`) and the §9 ORF crossing
(`vcr_array.cpp`). The ORF traversal navigates by **k-mer membership**
(`IndexBinarySearch`), not `OutgoingEdges` on a single edge-row — the latter only
follows one lineage at a fork in the succinct DBG and misses sibling branches.

**Verified against planted ground truth (synthetic, identical VCRs):**
- VCR: recovered **byte-identical** to the planted VCR; `cmsearch` E = 7e-35.
- ORF crossing, in-range: 8 datasets (seeds 1–5 ×6 cassettes; ×12; ×20; ×50) —
  **every planted ORF recovered byte-identical, 0 false positives** (`test_orf.sh`).
- §9.3 length gate: with cassettes spanning 177–1485 bp, the 4 out-of-[210,1200]
  ones are **rejected** and all 20 in-range recovered, 0 leaks (`test_reject.sh`).

**NOT yet verified (do not claim these work):**
- **Divergent VCR copies** — all tests use identical VCRs (the clean-collapse case).
- **Cassettes that share sub-sequences** (paralogs) — would create forks where the
  beam could chimerize; the synthetic ORFs are random/distinct, so untested.
- **Array bound [2,300] boundary rejection** (counts tested were 6–50, in range).
- **Cycle rule (§9.5)** — no cycles arise in the distinct-ORF synthetic, so the
  rule is implemented but unexercised.
- **Real metagenomic conditions** — uneven coverage, sequencing errors, strain mix.
- **Scaling** — the beam copies a per-path visited-set on each branch (O(path len));
  fine for tests, needs a persistent/immutable set for billion-node graphs.

**Superseded:** the old `vcr_graph.cpp` prototype (custom hash-map graph, wrong
anchor check, beam-not-greedy → chimeras, 0/8 recovered) has been removed.
