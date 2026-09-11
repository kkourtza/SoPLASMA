# History rewrites, and how to resolve a hash that no longer exists

## 2026-09-12 — run-output purge, to unblock the push

The branch could not be pushed: **seven blobs exceeded GitHub's 100 MB hard limit**
(largest 317.5 MB), all `validation/*/logs/log.soPlasmaFoam`, baked into 151 of the 208
unpushed commits. The `.gitignore` work of 2026-09-11 stopped them accumulating, but git
pushes HISTORY, so the whole push was rejected.

All run output was purged from the unpushed range with `git filter-repo`:
**2094 MB → 39.7 MB, 0 blobs over 50 MB.** The rewritten chain still descends from origin's
tip, so it went up as a normal fast-forward with no `--force`.

**The working tree was never changed.** The purge set was DERIVED, not hand-written:
`purge = output-shaped AND not present in the current tree`, so tree identity holds by
construction. Verified: tree `ce3db1a6f43dcaa42996bfb85bcf5b45bc0437d7` before and after.

*(A first attempt hand-wrote the path list and silently dropped 414 tutorial `plasmaTables`
files — which would have broken `positiveStreamer_LMEA_fast`, the ~2 s debugging bed — plus
14 `testSnesJFNK` mesh files that are the only copy. It was caught only because the tree
hash was compared before and after. The lesson is the one G1 already states: do not write a
second definition of something the repo already defines.)*

**152 commits changed SHA.** The 102 hash citations in docs and memories were rewritten in
the same change, so nothing in the live documentation is stale. But **commit MESSAGES are
immutable** — where one cites a pre-purge hash, it names a commit that no longer exists.

## Resolving an old hash

`commit-map-2026-09-12-log-purge.txt` maps `old new` (12 chars each), changed commits only:

```bash
grep ^d78a8d2 docs/history/commit-map-2026-09-12-log-purge.txt
```
