# seidr — VCR / super-integron detection on a de Bruijn graph

Detects *Vibrio cholerae* **VCR (attC) repeat–cassette arrays** directly from raw
reads, on **MEGAHIT's succinct de Bruijn graph (SDBG)** — no assembly.

The VCR forms a cycle `… STOP — VCR — START — ORF — STOP — VCR …`. Detection:
find a node with **in-degree ≥ 2** whose **3 upstream nodes spell a STOP codon**
(the VCR start), then a **greedy traversal weighted by the attC covariance model**
(Infernal CM, linearized to a per-position profile) walks out the VCR. Candidates
are validated with `cmsearch`. Method, step by step: **[PLAN.md](PLAN.md)**.

## Requirements (Linux / WSL)

`git`, `cmake` (≥3.12), `g++` (C++17), `make`, `zlib`, `curl`, `python3`.

## Quick start

```bash
bash pipeline.sh      # builds deps, generates data, compiles, runs — all of the below
```

## Step by step (same order the pipeline runs)

```bash
bash build_megahit.sh    # 1. clone + build MEGAHIT      -> megahit/build/megahit_core
bash build_infernal.sh   # 2. build Infernal (cmsearch,cmemit) -> deps/
bash generate_data.sh    # 3. VCR (cmemit) + CM profile + reads + SDBG -> build/
bash compile.sh          # 4. compile detector vs MEGAHIT's SDBG -> build/vcr_traverse
bash run.sh              # 5. detect on the SDBG + validate with cmsearch
```

## Outputs (in `build/`)

| file | what |
|---|---|
| `reads.fa` | synthetic reads (the input) |
| `graph.sdbg.*` | the MEGAHIT SDBG |
| `profile.txt` | linear per-position CM profile |
| `vcr.fa` | the planted VCR (ground truth) |
| `cands.fa` | recovered VCR candidate(s) |
| `cands.tbl` | `cmsearch` validation |

## Key files

| file | what |
|---|---|
| `models/Vibrionales.cm` | the attC covariance model (Infernal) |
| `cm_to_profile.py` | CM → linear per-position profile |
| `sdbg_gen.py` | synthetic super-integron generator |
| `vcr_anchor.cpp` | anchors: in-degree ≥ 2 + upstream STOP |
| `vcr_traverse.cpp` | anchors → greedy CM-profile walk → VCR candidates |
| `vcr_array.cpp` | VCR → START (ATG) → bounded beam ORF crossing → cassettes (PLAN.md §9) |
| `PLAN.md` | the method, in full |

## Tests (recovery vs planted ground truth)

```bash
bash test_orf.sh       # ORF crossing across seeds & cassette counts (6..50); byte-identical check
bash test_reject.sh    # §9.3 length gate: out-of-[210,1200] ORFs are rejected
```

Both compare recovered ORFs to the planted `build/truth.txt` byte-for-byte
(`check_orf.py` / `check_reject.py`) and print `OK`/`FAIL`. Verified status and the
remaining untested cases are in **[PLAN.md](PLAN.md)** ("Honest status").
