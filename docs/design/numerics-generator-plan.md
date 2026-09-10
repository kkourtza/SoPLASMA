# A semantic numerics layer, and a generated `fvSchemes`

Status: DESIGN, 2026-09-06. Nothing implemented. Queued behind the Grubert
loss-imbalance thread at the user's direction.

## The question that produced this

"Are the current cases running with SG flux or ROUNDF?" — and the answer took
reading three files, which is the defect. `fluxScheme standard` does not name a
scheme. It means *assemble as `fvm::div(phi,n)` + `fvm::laplacian(D,n)`*
(`driftDiffusion.C:314-318`), so the actual convective scheme comes from
`system/fvSchemes`, which in the Grubert cases says `Gauss ROUNDF` for every
species. **The cases are running ROUNDF, via `standard`.** `standard` and
`ROUNDF` are not alternatives; they sit on different axes.

`ScharfetterGummel` is a different *kind* of choice: it replaces the
div+laplacian pair with one exponentially-fitted two-point operator
(`driftDiffusion.C:304`), so it **bypasses `divSchemes` and `laplacianSchemes`
entirely**. A case can carry a carefully chosen `Gauss ROUNDF` line, switch to
SG, and lose it with no message. That silent trap is undocumented today.

## What is wrong with the present state

1. **`fluxScheme` names nothing.** It is an ASSEMBLY choice that decides
   whether the per-term schemes are used at all.
2. **`fvSchemes` is hand-authored, with regex catch-alls.** Nothing generates
   it — `fvSolution-foam -> fvSolution` is a `cp` in `Allrun`, and there is no
   equivalent for schemes. The file's own comment records the cost: a case
   whose ions were mobile **died at step 1** on `laplacian(D_Ar2p,n_Ar2p)`
   because only the electron was enumerated. A generator that knows the species
   list cannot make that mistake.
3. **The vocabulary is OpenFOAM's, not ours.** `divSchemes`,
   `laplacianSchemes`, `interpolationSchemes` mean nothing to a user who does
   not already know OpenFOAM, and G3 says a concept gets one name across the
   framework. `fvSchemes` is no better as a user-facing name, which is why it
   was rejected as the replacement for `standard`.
4. **A real footgun in the scheme names.** The built `libROUNDSchemes.so`
   registers EIGHT scalar names, not the four the source's
   `makeLimitedSurfaceInterpolationScheme` lines suggest: ROUNDA, ROUNDAplus,
   ROUNDF, ROUNDL **and ROUNDA01, ROUNDAplus01, ROUNDF01, ROUNDL01**
   (registered through `Limited01Limiter`). The `01` variants clamp the field to
   **[0, 1]**. For a mass fraction that is the point; for a number density of
   ~1e16 m^-3 it is catastrophic — and "bounded" is exactly what a user reaches
   for when told a density must stay positive. Verified 2026-09-06 by reading
   the symbols in the built `.so`, not the source.

## Proposal

### 1. Terms named by physics

A `numerics` block in the semantic layer. Names describe the TERM, not the
OpenFOAM operator:

    numerics
    {
        driftScheme                 ROUNDF;        // convection of a density
        diffusionScheme             harmonic;      // face value of D
        transportCoefficientScheme  harmonic;      // mu, D onto faces
        timeScheme                  backwardEuler;

        electron { driftScheme ROUNDAplus; }       // per-species override
    }

Every name is reused wherever the same concept appears — dictionary key,
start-up report, doc — per G3.

### 2. A curated list, PLUS a verbatim escape

The curated names carry what they are for and what is actually claimed about
them, with a citation (G2). Anything NOT in the curated list is passed through
to `fvSchemes` verbatim, so **every OpenFOAM scheme stays reachable** without
this framework enumerating OpenFOAM. That is what makes a short curated list
safe: it is a recommendation, not a restriction.

Claims must stay inside what the sources support. For the ROUND family that is
Deng 2023a (*J. Comput. Phys.* 112052, the theory) and Deng 2023b (*J. Comput.
Sci.* 102150, the library): "high-resolution structure-preserving",
"essentially oscillation-free", at "a minor increased CPU cost". **No accuracy
ORDERING between ROUNDF / ROUNDA / ROUNDAplus / ROUNDL may be shipped until it
is measured on a plasma bed** — the library's own benchmarks are wave
convection and a Mach-3 step, neither of which is a sheath.

The `01` variants must be **refused** for number densities with a message
saying why, not merely omitted: omitting them leaves them reachable through the
verbatim escape, which is where a user would land after reading "strictly
bounded" in the library README.

### 3. `system/fvSchemes` becomes generated

Best home is `foamPlasmaCreateSpeciesFields`: it already knows every species —
it created the 13 in the Grubert case — and already emits per-species files.
Generating one explicit entry per actual species removes the regex catch-alls,
the `#include` depth problem, and the die-at-step-1 class of failure.

Follow the boundaries-generator method that worked: **hand-author the target
first**, then verify the generator reproduces it byte-for-byte before switching
any case over.

### 4. Rename the assembly key

`fluxScheme` -> `fluxAssembly`, values `splitOperators` and
`ScharfetterGummel`, with `standard` accepted as a documented alias so no
existing case breaks. It is not a scheme; it decides whether the per-term
schemes apply.

### 5. Report it at start-up, and warn when a choice is dead

    drift-diffusion assembly: splitOperators
      n_e     drift ROUNDF   diffusion harmonic(0.33)
      n_Arp   drift ROUNDF   diffusion harmonic(0.33)
      nEps_e  drift ROUNDF   diffusion harmonic(0.33)

and under `ScharfetterGummel`, an explicit warning that the drift and diffusion
choices are UNUSED. A user should never have to read source to answer "what
scheme am I running".

### 6. Defaults, so most users write nothing

drift `ROUNDF`; diffusion `harmonic corrected 0.33`; coefficients `harmonic`;
densities `linear`. These are what the cases already use, so this moves
existing choices into a defaulted library rather than inventing any. Per G1 the
user should be able to write NOTHING and get these.

## What this must not become

More knobs is the opposite of G1. The resolution is that every one of these is
an OPTIONAL OVERRIDE over a documented default — a case that says nothing gets
good numerics, and an expert who needs `Gauss limitedLinear 1` is not blocked.
If the block ends up mandatory, the design has failed.

## Open, and needing measurement rather than assertion

* Which ROUND variant suits a sheath. Not answerable from the library's
  benchmarks; needs a plasma bed with ground truth.
* Whether `ScharfetterGummel` or `splitOperators + ROUNDF` is better at a
  sheath edge. At the Grubert cases' refined electrode cell Pe ~ 0.1-0.2, where
  the two boundary `f` forms nearly coincide, so this is a question about the
  INTERIOR, and it has not been measured.
