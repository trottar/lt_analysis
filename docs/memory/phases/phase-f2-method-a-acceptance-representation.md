# Phase F.2 Method-A acceptance representation audit

## Source implementation

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — F.2 started from the exact C2
baseline `02f21886cff8a721df46ebf673594856583fc3c5` after F.1 closure. Its
implementation commit and exact F.2.Fix.1 source starting HEAD is
`5549098d2552b9c092b65e31aeb77edc0807ddda`.

The implementation changes only:

- `src/cuts/pion_hgcer_method_a_acceptance_representation.py`;
- `testing/analyze_pion_hgcer_method_a_acceptance_representation.py`; and
- their two focused tests.

It is a JSON/numeric postprocessor with numpy, scipy, and matplotlib only. It
does not import ROOT or modify `rand_sub.py`, F.1 producers, normal analysis,
subtraction, correction, yield, or production paths.

F.2.Fix.1 is a strict F.1 authority/provenance boundary repair. Before any
response or support metric it reconstructs the producer's canonical-JSON
SHA-256 fingerprints for both serialized populations, exact feature metadata,
the application child-assignment projection, and the complete v2 fingerprint
inputs. It requires exact detached flags, `reason = null`, nonempty frozen
provenance, `source_target_state = post_proton_noRF`, exact frozen metadata,
and finite strictly increasing t/delta/phi edges with record-to-t-edge
agreement. The deterministic CLI validates the declared setting, kinematic,
and particle before associating each source hash, and rejects colliding output
paths. This repair changes neither F.2 science nor F.1.

## Frozen audit contract

- Input is exactly the five canonical Q4p4W2p74 F.1 v2 artifacts: Left-lowe,
  Left-highe, Center-lowe, Center-highe, and Right-highe.
- F.1 schema, identity, detached flags, feature metadata, positive-response
  training classification, physical application mask, and prompt identity
  closure are re-audited before metrics.
- The only candidates are `delta_only`, `track3`, `hgcer3`, and non-promotable
  `full5_reference`. Canonical t is a parent partition, never a predictor.
- Every setting x canonical-t group requires at least 25 low and 100 control
  training records. Five stable-identity folds use fold-local median/IQR
  scaling with standard-deviation fallback and class-balanced analytic-gradient
  L-BFGS-B logistic probes (`lambda=1e-3`, `maxiter=1000`, `gtol=1e-8`).
- Reduced bases must meet the fixed full5-relative information gates and the
  same-parent cKDTree non-prompt application OOD gate. Sparse non-prompt groups
  are explicitly flagged for manual review.
- Output persists aggregate metrics, fixed configuration, input fingerprints,
  F.1 source hashes, and deterministic fingerprints only. No score,
  probability, correction, C_A, adjusted weight, map, or child normalization is
  persisted or constructed.

## Local validation

- F.2.Fix.1 module/CLI/tests `py_compile` passed.
- F.2.Fix.1 focused suite: 19 tests passed, including stale F.1 fingerprint,
  authority metadata, geometry, deterministic-file identity, sparse zero-row,
  optimizer failure, and output-safety coverage.
- Frozen F.1 contract/runtime suite: 10 tests passed.
- `git diff --check` passed before memory updates.

## Farm gate

NEXT — Run only the detached analyzer from a clean checkout containing the F.2
implementation:

```tcsh
cd /group/c-kaonlt/USERS/trottar/lt_analysis
python3 testing/analyze_pion_hgcer_method_a_acceptance_representation.py \
  --outdir /lustre24/expphy/volatile/hallc/c-kaonlt/trottar/OUTPUT/Analysis/KaonLT \
  --kinematic Q4p4W2p74 \
  --output-json Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-representation.json \
  --output-pdf Q4p4W2p74_kaon_pion-background_hgcer_method-a-acceptance-representation.pdf
```

Collect JSON, PDF, stdout/stderr, `git rev-parse HEAD`, and `git status
--short`. Inspect all per-group gates, source hashes, and the recommendation.
F.2 cannot become `CLOSED / RUNTIME VALIDATED` and F.3 cannot begin until one
reduced basis is explicitly accepted.
