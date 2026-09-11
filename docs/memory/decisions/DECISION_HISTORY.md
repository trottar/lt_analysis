# KaonLT durable decision history

The source provenance for D001-D020 is the imported DECISION_HISTORY.md, with
farm-specific support linked below. [SUPERSEDED] decisions are retained because
they explain constraints and rejected paths; only their stated later decision
is current.

## Production background ownership

### D001 — slow protons and pion production are distinct

[SCIENTIFIC_REFERENCE] [CHAT_HISTORY] Do not collapse PID contamination in
kaon-selected data and pion-production background mapped from pion-control data
into one empirical missing-mass correction.

### D002 — proton cleaning precedes pion subtraction

[SCIENTIFIC_REFERENCE] [HANDOFF] Preserve random/dummy -> frozen binning
spectra -> slow proton -> pion subtraction.

### D003 — K-Lambda preservation is a setting-wide safety gate

[HANDOFF] [FARM_EVIDENCE] Compare proposed removal in protected K Lambda to the
default 10% limit. On failure, bypass cleaning for the entire setting
(applied=0, cleaned factor=1) while retaining proposed diagnostics. Do not
partially apply per-t results.

### D004 — proposed and applied proton semantics differ

[HANDOFF] Applied means the post-gate, committed production quantity; proposed
model output remains diagnostic when the setting gate rejects it.

### D005 — pion-alignment comparison uses a fixed envelope

[HANDOFF] Candidate fit windows may vary for amplitudes, but every candidate,
including the explicit baseline, is ranked on common evaluation bins with
support/pathology diagnostics.

### D006 — HGCer consumes, not redefines, the pion baseline

[CHAT_HISTORY] [HANDOFF] Do not change pion components merely to make a later
HGCer diagnostic agree.

## HGCer diagnostic independence and persistence

### D007 — Method A and Method B are independent

[HANDOFF] Method A is detector-response driven; Method B is frozen local
missing-mass closure. No Method-A number belongs in Method B.

### D008 — Method B is same-canonical-t relative closure

[HANDOFF] [FARM_EVIDENCE] Normalize only within the same canonical t parent.
Do not interpolate or pool across t, and do not use the neutron peak as an
absolute pion-normalization anchor.

### D009 — adaptive Method B is DO NOT PROMOTE

[FARM_EVIDENCE] It was stable where supported but too sparse and
setting-dependent. Legacy Method B remains the Phase-D reference; adaptive B
is diagnostic context only. No C.Fix.2.4, D.Fix.1, or production change
follows. See evidence/phase-c-five-setting-closure.md.

### D010 — persisted diagnostics require round-trip review

[FARM_EVIDENCE] Trace producer -> serializer/checkpoint -> checkpoint-first
payload -> consumer -> renderer; a live producer alone is insufficient.

### D011 — Phase E owns presentation only

[HANDOFF] Phase E consumes frozen values; it cannot recompute Method A, inject
Method-B numerical input, or construct a correction.

### D012 — generic collector plus JSON profile is validation architecture

[FARM_EVIDENCE] Add collector capabilities only when necessary; use
gate-specific profiles for different declared review gates.

### D013 — validation bundles include review artifacts

[FARM_EVIDENCE] Package declared pages, checkpoints, and provenance so farm
review is reproducible rather than manually reconstructed.

### D014 — source identity is a tooling responsibility

[FARM_EVIDENCE] Validation tooling must machine-check evaluated source
provenance and reject stale or unexpected analysis changes.

## Method-A future role and Phase-F boundaries

### D015 — an A+B correction plan is superseded

[HANDOFF] [SUPERSEDED] Earlier procedure language suggested a future A+B
correction. The later authoritative decision is Method B as
diagnostic/cross-check/historical comparison only, and Method A as the sole
candidate future numerical HGCer leakage input.

### D016 — current Method A is not absolute leakage probability

[HANDOFF] Because it excludes NPE=0, Method A is a positive-response relative
leakage diagnostic. Stronger absolute language requires a dedicated
zero/nonpositive-response checkpoint.

### D017 — Method-A training and application are distinct populations

[FARM_EVIDENCE] [HANDOFF] Train from prompt/noRF/nommcuts positive-NPE response
records; retain the authoritative physical NPE>2 pion-control cache for
application. The F.1 v1 farm result proved the latter cannot own response
training.

### D018 — response and yield coordinates intentionally differ

[HANDOFF] Learn response in canonical t plus acceptance coordinates, then apply
event by event before canonical (t,phi) child filling.

### D019 — normalize at the parent, never each child

[HANDOFF] Preserve each canonical-t parent baseline normalization; do not
independently renormalize (t,phi) children.

### D020 — no production promotion before F.6

[HANDOFF] F.1-F.5 are detached diagnostic/development stages. Production
promotion requires a separately validated and explicitly approved F.6 decision.

