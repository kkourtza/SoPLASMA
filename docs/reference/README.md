# Option reference

Organised by **the file you edit**, because that is how you arrive here: you
open a dictionary, you want to know what a key does and what it defaults to.

| you are editing | reference |
|---|---|
| `configuration/config` | [configuration-config.md](configuration-config.md) |
| `system/controlDict` | [controlDict.md](controlDict.md) |
| `system/plasmaSimulationControls` | [plasmaSimulationControls.md](plasmaSimulationControls.md) |
| `system/fvSchemes`, `system/fvSolution` | [fvSchemes-fvSolution.md](fvSchemes-fvSolution.md) |
| `constant/regionProperties`, `constant/<region>/electricalProperties` | [regions-and-materials.md](regions-and-materials.md) |
| `constant/plasmaSpeciesProperties` | *pending* |
| `constant/plasmaTransportProperties` | *pending* |
| `etc/changeDictionary.<region>` | [changeDictionary.md](changeDictionary.md) |

## What goes in `configuration/config`, and what does not

`config` is **not** a place to put every number. A value used in ONE dictionary
does not belong there: `$var` in the dictionary plus a declaration in `config` is
**two files for one value**, which is worse than writing it where it is read.

Measured on the needle-DBD case: 22 of 37 variables were used exactly once. They
were moved into the dictionaries that read them, and `config` went from 37 to 15.

**Two things earn a place in `config`:**

1. **Values used in more than one place.** The five scheme variables are used
   **five times each** — once per mobile species — so one decision lives in one
   place instead of five. That is the case for indirection made concretely.
2. **The headline physics inputs**, even when used once: the drive waveform, the
   gas state. These are what a user changes between runs, so they belong
   together at the front rather than scattered across `system/` and `constant/`.

Everything else goes in its own dictionary.

## Why this reference exists separately from the case files

Two kinds of comment, with different lifecycles:

- **Reference** — what an option is, its type, its default, the menu of
  alternatives. Identical in every case, so it lives here, once.
- **Rationale** — why *this* case chose *this* value, and what was measured
  here. Unique per case, so it stays in that case's file.

Embedding the reference in each case means N copies drifting independently. That
failure has been removed three times in this repository:
`electromagneticsProperties` copied per region so every copy held every other
region's data; `fvSchemes` copied per region, replaced by `#include` stubs; and a
`tools/` directory copied into a case that was 117 lines stale within the hour.

## Keeping it honest

```bash
tools/checkConfigReference.py                 # docs vs src: missing readers, stale defaults
tools/checkConfigReference.py --case <dir>    # a case's config: dangling, wired-to-nothing
tools/checkConfigReference.py --strict        # undocumented options become fatal
```

It catches the class of error no careful reading finds: an option documented in
prose that **nothing reads**. `ePotentialControls { nonCoupledResidualControl }`
was documented in two files here and has never had a reader.
