Hey Michael,

Picking up on your "I always forget to test on Bruker data" — rather
than wait for a fix, I've put together some simulated `.d` fixtures
you can drop straight into Sage's CI. They come from TimSim, so
ground truth is exact: every `(sequence_unimod, charge)` precursor
that was injected is recorded in a `synthetic_data.db` next to the
`.d` directory.

There are three pieces:

**1. CI-grade smoke fixture (~235 MB / replicate, 2 replicates).**
Tiny HeLa DDA run — 1 500 peptides, 5-min gradient. Two `.d`
directories so the multi-input code path gets exercised (which is
exactly what tripped the canonicalize bug in this issue —
`processing files 0 .. 1`). Fast enough to run end-to-end in a
GitHub Actions job under a couple of minutes.

**2. HeLa 150K G30M evaluation set (~6 GB).** Full-size HeLa DDA
from the same recipe used in our paper's MBR section. Useful as a
once-per-release sanity check or for reviewers wanting to compare
Sage to other DDA-PASEF tools.

**3. HLA peptide eval sets (10K and 100K, ~825 MB / 3 GB per rep).**
From the existing TimSim-HLA Thunder TopN compilation. Good for
confirming Sage handles short, low-charge immunopeptidome data
correctly.

A small Python harness (`sagebench`) joins Sage's parquet output
against the ground-truth DB and computes true FDR / TPR — same
definition as in the TimSim paper's methods section, so it's a fair
yardstick across search engines.

I confirmed end-to-end on my side: applied a one-liner fix to
`read_tdf` (mirroring the file-URL → local-path workaround that
`bench/matteo-pmsms` already applies to `read_pmsms`) and re-ran
sage 0.15.0-beta.2 against the three datasets. Numbers at 1 %
peptide-FDR (`peptide_q ≤ 0.01`):

```
HeLa CI smoke (2 × .d)   :   871 precursors, true FDR 0.23 %, TPR 60.0 %
HeLa 150K G30M           : 3 230 precursors, true FDR 3.59 %, TPR  3.0 % (frag)
HLA 10K G40 (unspecific) : 2 233 precursors, true FDR 2.42 %, TPR 25.1 % (frag)
```

(`TPR (frag)` is over the DDA-fragmentable subset — TimSim
precursors that got an MS2 spectrum at all. The denominator most
people would call honest. `TPR (full)` against every simulated
peptide is much lower for the larger runs.)

---

One thing that came out of running this that I wanted to flag (not
a sage bug, but worth knowing): on the simulated runs, every
semi-supervised rescorer I tested — sage's per-run LDA, mokapot's
linear SVM, and mokapot's XGBoost — converges on ~3.4–3.8 % *true*
FDR at the 1 % claimed peptide-FDR cutoff. Re-rescoring the same
parquet with raw hyperscore + classical TDC (`sagepy.qfdr.tdc`)
brings it back to ~1.6 %, well within target. So:

- The drift is rescorer-family-agnostic.
- Feature-ablation (drop top-weighted by |LinearSVC coef|, retrain)
  is diffuse — no single feature is the leaker; removing matched_peaks
  and hyperscore doesn't move FDR; the model just re-routes weight.
- Protein-level FDR is unaffected (sage native protein_q gave
  0.05 % true on HeLa G30M) because protein grouping aggregates
  away peptide-level miscalibrations.

Whether this is **specific to simulated data** (the simulator's
target/decoy feature gap is slightly wider than on real DDA) or a
**general property of these rescorers** (only visible on simulated
data because that's where you have GT) — I can't tell from
simulated data alone. Distinguishing them needs an entrapment-FDR
experiment on a real Bruker `.d`. Worth flagging so a future user
of TimSim-as-CI-fixture doesn't read the `peptide_q` column at face
value.

Repo (work-in-progress, will polish + open up before posting publicly):
`<insert SAGEBench URL once published>`. The smoke fixture and
generation script are reproducible from a TimSim build (the
`from_findings` code path keeps it deterministic, no FASTA download
needed). The eval sets will go to a dedicated Zenodo deposit so they
don't bloat the repo.

If you'd like, I can open a PR wiring the smoke fixture into Sage's
CI directly — just point me at where you'd want the workflow to live.

Best,
David

---

<!--
Editing notes for David before posting:
  - Replace `<insert SAGEBench URL once published>` once the repo is
    public.
  - Sizes are real now (235 MB / smoke rep; 6 GB / HeLa 150K G30M;
    825 MB / HLA 10K rep; 3 GB / HLA 100K rep) — confirmed on the
    actual generation runs and Zenodo bundles.
  - The "rescorer drift" paragraph is the bit that requires the
    most editorial judgement. Frame conservatively if you don't
    want it to read as a sage criticism — it's really a finding
    about ML rescoring on simulated DDA in general, not about the
    sage implementation.
  - The full HTML report lives at runs/REPORT.html. Include or link
    it depending on whether you want Michael to see Venns + TDC
    plots inline.
-->
