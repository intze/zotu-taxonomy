# Project context for Claude

Species-level taxonomy for zOTUs (UNOISE3, usually V3–V4 16S amplicons, 350–500 bp) as a
free, batch replacement for manual EzBioCloud lookups. The user (Greek-speaking; reply in Greek)
runs it on Windows via the GUI. See README.md for usage and method.

## Files
- `zotu_taxonomy.py`: core pipeline (`run_pipeline`) + CLI. vsearch `--usearch_global` against
  NCBI RefSeq 16S Targeted Loci (bacteria + archaea, ~28k type-strain seqs), best hit per species,
  top-N species per zOTU, Excel output (Summary / Top hits / Settings).
- `zotu_gui.py`: Tkinter GUI; `--selftest IN.fasta OUT.xlsx` runs headless (used by CI).
- `.github/workflows/build-windows.yml`: builds `ZotuTaxonomy.exe` (PyInstaller, onefile, windowed)
  on windows-latest, bundles vsearch in `tools/`, self-tests the exe on `tests/example_zotus.fasta`,
  and replaces the GitHub release `latest` on every push to `main`.

## Verified so far
- Primer remnants: zOTUs keep the 341F tail `GCWGCAG` (7 bp) and 805R tail (9 bp). The degenerate W
  gave one artificial mismatch per zOTU (100% hits went 3 → 134 per 300 after trimming).
  Auto-trimming is on by default.
- Heuristic search (`maxaccepts 200 / maxrejects 2000`) matched exhaustive search
  (`0/0`) for all identities ≥ 94.5% on 150 random zOTUs; small differences only below 92%.
- User dataset (3,293 zOTUs): Species 478, Ambiguous 120, Genus only 436, below genus 1,382,
  no type strain ≥ 80% 877. The last group is real: 60–80% to nearest type strain
  (e.g. Patescibacteria/Minisyncoccus, Mycoplasma); not a bug.
- Runtime ≈ 10 min for 3,293 zOTUs on 4 cores.
- GUI tested on Linux under Xvfb (run, error dialog, double-click guard). Greek text could not be
  rendered there (Tk without Xft); not yet confirmed on real Windows.

## Open decisions / next steps
1. Tie tolerance default (`--tie-delta`, GUI "Ανοχή ισοπαλίας") is 0. Recommended 0.3
   (≈1 nt at 430 bp). On the user's data: 0 → 478 Species, 0.3 → 418, 0.5 → 369, 1.0 → 295.
   Waiting for the user's choice.
2. Compare with EzBioCloud on 10–20 zOTUs the user already identified there.
3. Confirm the Windows exe opens and shows Greek correctly on the user's PC.

## Environment notes (Claude Code cloud sessions)
- `ftp.ncbi.nlm.nih.gov` must be in the environment's allowed domains to download the database.
- Pushing git tags was blocked from the session; releases are therefore published by the workflow
  from `main`. GitHub Actions logs/artifacts are not downloadable from the session (proxy 403);
  release assets are, via `gh api -H "Accept: application/octet-stream" repos/Intze/zotu-taxonomy/releases/assets/<id>`.
- System Python has no tkinter; for GUI tests use `uv python install 3.12` and `xvfb-run`.
- Never commit user data (FASTA, Excel); `.gitignore` excludes outputs, `db/`, `tools/`.
