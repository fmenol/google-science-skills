# Worked Example: Ubiquitin Tracks (1UBQ)

Extract ESM3's five conditioning tracks from a real experimental structure,
inspect them, and mask them into a scaffolding prompt. The subject is **PDB
1UBQ** — Vijay-Kumar's 1.8 Å human ubiquitin, chain A, 76 residues — because it
is small, fully resolved, and a textbook beta-grasp fold, so every track has an
unambiguous right answer.

## What this example shows

-   `extract` turns 1UBQ into the sequence, coordinate, SS8 and SASA tracks with
    no model, no API key and no credits — track preparation is pure local
    computation.
-   `inspect` summarises the tracks and renders the SS8-ribbon-over-SASA figure.
-   `build-prompt --keep 10-20` masks the tracks into an ESM3 prompt, keeping a
    short motif and destroying everything else — the exact input a design run
    consumes.

## Reproduce it

All three commands are local and free (1UBQ is fetched once from RCSB and then
cached, so re-runs are offline):

    uv run --no-project scripts/tracks.py extract \
      --pdb-id 1UBQ --chain A --output ./out/1ubq_tracks.json

    uv run --no-project scripts/tracks.py inspect \
      --tracks ./out/1ubq_tracks.json \
      --output ./out/1ubq_inspect.json --plot ./out/ubiquitin_tracks.png

    uv run --no-project scripts/tracks.py build-prompt \
      --tracks ./out/1ubq_tracks.json --keep 10-20 \
      --output ./out/1ubq_prompt.json

## Files

-   [report.md](report.md) — the tracks, the prompt, and the conventions that
    matter, quoting the real numbers.
-   [ubiquitin_tracks.png](ubiquitin_tracks.png) — the SS8 ribbon over the SASA
    profile, straight from `inspect`.

## The one thing to take away

`--keep 10-20` keeps **11** residues (the 10th through the 20th), not 10 and not
residues 11–21. Ranges are 1-indexed and inclusive. Trust the script's
`kept_positions`; never hand-count.
