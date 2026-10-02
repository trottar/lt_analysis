# E.8.4 Fix.5.4 — current-lineage numerical identity and SIMC normalization audit

**Status:** `DEVELOPMENT COMPLETE, FARM VALIDATION PENDING`.
Numerical source is `SOURCE REVIEWED` at committed/pushed `test` HEAD
`761fbb6c03d2d7a10bb911cf84e9ba898496fab6`, parent
`71a912bf14c59eee6d52956462aa0aa39c7f7c84`. Independent ChatGPT actual-diff
and pushed-state review passed, as supplied by the Fix.5.5 task contract.
The earlier implementation began at that parent, whose parent scientific
source is `f9d70732290ea461096374ca1270b47452644991`. Review/push establishes
source acceptance only; Fix.5.4 has no new farm numerical closure.

## Ownership and implementation

The [implementation contract](e8-4-fix5-4-current-lineage-identity-audit-and-fix-task-contract.md)
owns this numerical/provenance step. `calculate_yield.py` attaches one private
`f6_3_current_lineage_identity_audit/v1` record after public baseline yield
extraction. It reads the exact measured final histogram and existing detached
F.6.3 children. It validates geometry, content/error fingerprints, binwise
pion/MM algebra and histogram/scalar closure with fixed scaled `1e-12` identity
tolerance. Normal MM bins are summed in order without width weighting;
underflow/overflow are excluded from arithmetic but included in ownership
fingerprints. Signed positive, negative, signed and absolute support is retained
for baseline, Method A and their difference, including sparse and zero children.
No scalar or histogram is replaced to force closure.

Current aggregates sum only current canonical children without normalization,
retain complete/populated inventories and exact candidate authority. Candidate
F.4 parent sums include broader application support; cut-window pion-template
integrals are not asserted equal to those sums. Historical E.8.3/F.6.1
aggregates remain separate and cannot satisfy this record.

`full_background_subtraction_plots.py` checks record fingerprints, current
authority, inventories, object fingerprints and recorded arithmetic; clones
objects and detached records; rejects missing, malformed, stale, historical or
inconsistent identities. Explanatory wording explicitly separates current
aggregates from historical E.8.3. Existing page IDs/styles are retained for
consistent, comparable synthetic inputs.

## Source findings and scientific blocker

The active baseline snapshot is made immediately after pion subtraction.
Both subsequent empirical fit/subtraction/pruning blocks are conditional and
disabled in the active profile; final `H_MM_DATA` is cloned by `bin_data` and
measured by `integral_with_stat_error`. Equality is enforced against that actual
measured clone, rather than inferred solely from zero background scales.
No source-proven wrong-object linkage repair was warranted.

The existing SIMC producer fills same-cell `h10.missmass` with `iter_weight`
after cuts, scales once by stored `normfac_simc = simc_normfactor / Ncontribute`,
then clones to `_xsect_support_simc.mm`. Model iteration changes previous event
weight by the ratio of model cross sections. Its luminosity/charge units are
not declared by the traced source, while data yields use effective-charge
normalization. The audit copies source/root/tree identity, setting/geometry,
actual normalization inputs/factor, model semantics, histogram meaning,
fingerprints and normal-bin window integrals. It records absolute comparison
unavailable with `SIMC_normfac_luminosity_and_charge_units_not_source_proven`.

The [Fix.1 repair](e8-4-fix5-4-fix1-simc-availability-and-manifest-scope-task-contract.md)
preserves overall availability of the valid current-F.6.3 E.8.4 payload and
detached numerical/aggregate/SIMC audits. Only the two absolute-SIMC page
families are provenance-blocked when those units are unproven, retaining page
IDs, canonical inventories and literal reasons. Data-only/identity/yield and
parent-closure evidence remains renderable. No renderer scale, unit conversion, setting-wide fallback or
physics conclusion is introduced. This explicit unavailable scientific result
is within the contract; it does not establish incorrect SIMC normalization.
See the [source trace](../investigations/e8-4-fix5-left-lowe-post-farm-identity-lineage-visualization-2026-10-02.md).

## Local verification and boundaries

Python 3.12.10 was discovered locally. Deterministic focused numerical,
provenance, E.8.4, plotting, F.6.3, E.8.2 yield, ownership, parent-correction
and propagation tests pass; retained skip reasons distinguish unavailable
PyROOT from superseded historical procedure-tail assertions. Exact final
counts and memory health are recorded in the investigation. Syntax checks and
`git diff --check` pass. These are local checks only, not ROOT/PyROOT, rendered
PDF, full `main.py`, farm closure or actual observed signed-cancellation proof.

F.4.Refresh.2 and accepted historical/narrow runtime closures remain unchanged.
Method A remains detached/non-production; Method B remains numerically absent.
Frozen factors, weights, cuts, windows, fits, normalization, yield arithmetic,
production objects and uncertainties are unchanged. During Fix.5.4 implementation no SIMC helper, owner,
profile, collector, visualization style or farm execution changed. The user
subsequently committed/pushed the reviewed numerical source.

The separately contracted [Fix.5.5 presentation](e8-4-fix5-5-current-lineage-visualization-clarity.md)
now consumes these stored records from the reviewed numerical source, retaining
the absolute-SIMC blocker. Numerical ownership, validators and producer bytes
remain unchanged. CURRENT owns the sole NEXT; no farm execution is authorized.
