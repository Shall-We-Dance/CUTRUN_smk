# Replicate consensus, optional IDR, and shared fragments

These rules distinguish three questions: which **pooled candidate** is
replicate-supported, which **bases** are jointly covered by broad peaks, and
whether **ranked narrow peaks** are reproducible by pairwise IDR. They are
different outputs, not interchangeable definitions of consensus.

## Configuration

All settings below are optional. These are the defaults; set `idr.enabled:
true` only for suitable narrow-peak samples with at least two biological
replicates.

```yaml
consensus:
  min_replicates: majority       # floor(n/2)+1; alternatively all or integer 2..n
  narrow:
    match: reciprocal_overlap   # alternatively summit_distance
    min_overlap_fraction: 0.5
    summit_distance: 100         # bp, used only with summit_distance
  broad:
    min_coverage_fraction: 0.5
    min_core_length: 1           # bp after splitting on support identities
idr:
  enabled: false
  threshold: 0.05
  relaxed_pvalue: 0.1
  rank: p.value                 # alternatively signal.value
  min_peaks: 20                  # safety floor, not a biological QC threshold
```

Each sample can override `consensus` and `idr`; nested `consensus.narrow` and
`consensus.broad` settings inherit unspecified global values. For example,
add these keys alongside a narrow sample's existing `replicates`:

```yaml
consensus:
  min_replicates: 2
  narrow:
    match: summit_distance
    summit_distance: 100
idr:
  enabled: true
```

Use biological replicates, not sequencing lanes, as support units. The default
is 2/2, 2/3, 3/4, 3/5, etc.; a shallow replicate is not automatically discarded.
An explicit threshold greater than the number of replicates fails validation.
There is no reproducibility claim for a single-replicate sample.

## Candidate boundaries and matching rules

1. Call peaks in each replicate and in pooled fragments, using the sample's
   narrow/broad MACS3 mode and matched controls, when configured.
2. Use the **unaltered pooled peak intervals** as candidates. A candidate that
   passes support retains exactly those boundaries. Accepted candidates are
   not merged, extended, or replaced by a union of replicate peaks.
3. Evaluate each replicate independently; one replicate contributes at most
   one support vote, even if several of its peaks match.

Narrow matching has two explicit alternatives:

- `reciprocal_overlap`: at least one replicate peak must overlap by at least
  `min_overlap_fraction` of **both** its own length and the candidate length.
  A 1-bp touch does not qualify unless both intervals are themselves tiny.
- `summit_distance`: at least one replicate summit must be within the
  configured distance of the candidate summit, inclusive. Both files must
  contain valid narrowPeak summit offsets. This rule compares positions and
  does **not** also require interval overlap.

Broad candidate support is the number of bases covered by the **union of that
replicate's peaks**, clipped to the candidate, divided by candidate length.
It must meet `min_coverage_fraction`. Overlapping peaks in the same replicate
never double-count bases. This is a candidate-coverage rule, not reciprocal
overlap and not a claim that every candidate base has k-fold support.

All intervals are 0-based, half-open. These matching definitions are explicit
pipeline choices, not an implementation of the full ENCODE workflow.

### Outputs

Under `macs3/consensus/`:

| Output | Meaning |
| --- | --- |
| `<sample>.consensus.bed` | Accepted full pooled candidates, valid BED6 |
| `<sample>.support.tsv` | All candidates, including rejected candidates |
| `<sample>.core.bed` | Broad-only jointly covered core segments, BED6 |
| `<sample>.core.support.tsv` | Broad core support identities and candidate IDs |

The support TSV includes candidate/source peak IDs, support count, JSON
replicate-name list, required support, pass/fail, mode, matching rule, per-
replicate coverage/matched peak names, and the original pooled peak fields.
Candidate IDs link the tables to BED6 names; BED score is deliberately zero,
not an invented pooled statistical score. Supporter lists are deterministic
and lexicographically ordered.

Pairwise any-overlap counts and base-pair Jaccard remain separate QC in
`qc/replicates/`. Empty-vs-empty Jaccard is `NA`. A status of `ok` means the
calculation completed, not that the experiment passed biological QC.

## Broad jointly covered cores

The core sweep operates at **basewise** resolution over the union of each
replicate's broad peaks. Keep only segments covered by at least k distinct
replicates and located inside a pooled candidate.

- Split whenever the supporter **identities** change, even if their number is
  unchanged. Split also when the enclosing candidate-ID set changes.
- Merge adjacent pieces only when both identity sets are identical.
- Never bridge an uncovered gap or use a long peak to connect separately
  supported ends.
- Apply `min_core_length` to the resulting constant-support segments.

Core eligibility is **independent of whole-candidate coverage support**.
Pooling can broaden a candidate so much that its overall coverage fraction
fails, while a smaller region is genuinely shared. Such a core is retained
and its candidate ID points to the rejected row in `support.tsv`. Conversely,
a supported full candidate need not contain a k-fold core: different
replicates may support different parts. Use the core BED when downstream
analysis requires bases with explicit joint coverage.

This is peak-interval support, not a signal-intensity concordance test; inspect
individual replicate tracks and appropriate QC as well.

## Optional narrow-peak IDR

Enabling IDR adds **separate relaxed MACS3 calls** (`-p relaxed_pvalue`) for
every replicate and the pool. It never feeds the stringent consensus BED into
IDR. Relaxed calls use the same fragment tables and matched control mappings
as regular calls.

For each unordered pair of biological replicates, IDR uses the relaxed pooled
peak list as its candidate list and the configured ranking field. The wrapper
uses a fixed random seed and emits:

- `macs3/idr/<sample>/pairs/pair_NNN.idr.tsv`: original IDR output;
- `pair_NNN.narrowPeak`: true ten-column narrowPeak records passing the
  **global** IDR threshold (raw IDR column 12 is -log10(global IDR));
- pair logs and `macs3/idr/<sample>/summary.tsv`, mapping pair IDs to replicate
  names, input counts, shared candidates, passing counts and status.

Fewer than `min_peaks` peaks in either replicate or too few shared pooled
candidates produces empty pair outputs and `insufficient_peaks`. A genuine
IDR fitting error fails the job rather than masquerading as reproducibility.
The floor of 20 only guards obvious underpowered inputs; it does not ensure
a stable fit, especially with tied scores or low-complexity peak lists.

With >2 replicates all pairs are reported separately; the workflow does not
pretend there is a single multivariate IDR score. It also does not implement
ENCODE's pseudoreplicate, self-consistency, rescue-ratio, or optimal/conservative
peak-selection workflow. Broad peaks cannot enable this module. See the
[IDR implementation](https://github.com/nboley/idr) for the statistical model.

## One shared fragment population

`fragments/<sample>.<rep>.fragments.bed` is a sorted, persistent **BED3**
table. Every row is one accepted molecule; identical coordinates from
independent molecules are retained. Flat samples omit the replicate suffix.
MACS3, FRiP, raw/normalized/spike-in BigWigs and usable-fragment QC all consume
this same table.

Preparation starts after the configured MAPQ, blacklist and duplicate
processing of the standard BAM:

- PE: require exactly one primary R1 and R2 in each query-name/read-group
  template, both mapped and marked proper, both passing exclusion flags, on
  the same retained chromosome. Secondary/supplementary records do not create
  extra fragments. Missing mates and ambiguous primary pairs are rejected.
- PE coordinates span both aligned mates, calculated from their CIGARs;
  stale TLEN and a remaining proper-pair flag on an orphan cannot inflate
  fragment counts. Reference-consuming CIGAR operations determine the span.
- SE: extend from the aligned 5' end by `fragments.se_extension`, falling back
  to the sample's `macs3.extsize` (200 if absent). Clip at chromosome edges.
  This is an estimated span, not a measured paired-end insert. MACS3
  `nomodel: false` is rejected for SE.
- When blacklist filtering is enabled, reject any fragment **span**
  overlapping a blacklist interval, including intervals between its mates.

Optional shared settings:

```yaml
fragments:
  exclude_flags: 1804
  remove_chrM_and_scaffolds: true
  se_extension: 200              # SE only; omit to use sample macs3.extsize
```

If absent, `exclude_flags` defaults to the union of the global
`samtools_exclude_flags` and legacy `bigwig.samtools_exclude_flags`;
`remove_chrM_and_scaffolds` inherits the legacy BigWig setting (default true).
These are additional eligibility checks: they cannot restore reads already
removed upstream. The contig mask removes names with underscores, chrM/MT/M
(case-insensitive), and scaffold/random/alt/fix/hap names. It now applies
consistently to **all shared quantitative outputs**, not just BigWigs.

MACS3 reads these coordinates with `-f BEDPE --keep-dup all`. MACS3's BEDPE
format means **three columns describing the full fragment**, not bedtools'
ten-column paired-alignment format; see the
[MACS3 input-format documentation](https://macs3-project.github.io/MACS/docs/callpeak.html).
Do duplicate removal upstream only. Pooled tables concatenate/sort accepted
replicate fragments, without deduplicating identical positions across
biological replicates. Retained BAMs remain alignment archives and can contain
records excluded from the shared table.

### Coverage and count units

Each fragment contributes **one count to every bin it overlaps**, regardless
of the number or overlap of its read ends. This is fragment-overlapping-bin
coverage, not fractionally weighted per-base averaging within a bin.
`bigwig.bin_size` defaults to 10. At bin size 1 it is exact per-base fragment
depth. The track writer streams sorted fragments using bounded local-depth
memory; it does not allocate a whole-genome array.

| Normalization | Multiplier on fragment-bin counts |
| --- | --- |
| `None` | 1 |
| `CPM` | 1,000,000 / accepted fragments |
| `RPKM` (legacy option name) | 1,000,000,000 / (accepted fragments × bin size) |
| `BPM` | 1,000,000 / total fragment-bin incidences |
| `RPGC` | retained chromosome-size sum / sum of accepted fragment lengths |

The legacy `RPKM` option is now **fragment-RPKM (FPKM)**, not read-end RPKM.
The final partial bin uses the configured bin width for this normalization.
RPGC is exact 1× scaling at bin size 1; wider bins introduce endpoint rounding.
Its genome size is the retained reference length, not a mappability-adjusted
effective size. Empty fragment sets produce explicit zero tracks.

FRiP is accepted fragments overlapping at least one corresponding MACS3 peak
divided by all accepted fragments; each fragment counts once. The workflow
uses the per-replicate or pooled peak set, **not** the consensus/core/IDR set.
Read-end-only FRiP fields are `NA` because BED3 spans do not reconstruct read
end positions. The standalone legacy `frip.py --bam` interface remains
available; workflow runs use `--fragments`.

Spike-in scaling counts validated templates using the same PE pair/flag
eligibility, instead of halving a read count. Target-only chromosome masks
and blacklist intervals are not applied to the separate spike-in genome.
The aligned spike-in BAM retains its existing processing policy (no target
MAPQ/dedup filter); this does not imply biologically identical filtering of
the two references. The scale is `scale_to / spikein_fragments`.
Zero spike-in retains the legacy scale of 1 with `NO_SPIKEIN` status, **not**
a calibrated abundance estimate; low counts are flagged `LOW_SPIKEIN`.

Pooled CPM/RPKM tracks are depth-weighted, not equal-replicate means.
Individual replicate tracks remain available. The all-alignment TE diagnostic
is intentionally unchanged and is not a quantitative fragment/family track.

## Migration and testing

- Use a **fresh output directory** when upgrading an analysis. Fragment
  eligibility, SE extension, normalization denominators and consensus
  boundaries change scientific semantics; do not mix old and new tracks.
- Consensus is now `.consensus.bed` (BED6), not a three-column file falsely
  named `.narrowPeak`/`.broadPeak`. Update downstream paths explicitly.
- The default support threshold is majority, not all replicates. Use
  `min_replicates: all` for strict support, but note that the new matching and
  boundary definitions still differ from the old 1-bp-overlap union.
- A large pooled broad domain can fail candidate support and still have a
  valid smaller core. Do not interpret the two BEDs as nested pass/fail lists.
- Keep `macs3.keep_dup: all`; other values now fail validation. Explicitly
  review any inherited chromosome/exclusion settings and SE extension length.
- Peak calling with empty treatment fragments emits empty outputs with a
  no-fit log. A configured but empty control with nonempty treatment fails
  rather than silently becoming an uncontrolled call.

The README lists unit and end-to-end test commands. The three-replicate broad
fixture verifies a retained joint core even when the enlarged pooled
candidate fails whole-domain coverage. Its narrow branch exercises optional
IDR's insufficient-data result; the separate IDR unit test runs a real fit on
400 synthetic ranked peaks. The spike-in fixture reuses the tiny target index
only to test plumbing and MultiQC, not biological calibration. No synthetic
test establishes performance or scientific validity on a real experiment.
