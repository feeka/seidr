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

- Not implemented yet against this plan.
- The earlier prototype (`vcr_graph.cpp`) did **not** follow this plan: it used a
  custom hash-map graph (not the SDBG), added a wrong "anchor begins with a stop"
  check, and used a beam search (not greedy) — which produced chimeras. Verified
  failure: 0/8 planted copies recovered at ≥90 % identity.
- Whether the greedy SDBG method in this plan recovers divergent copies cleanly
  is **unverified** — it will be tested against planted ground truth (Step 8.5)
  before any success claim.
