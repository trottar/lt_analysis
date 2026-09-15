# Phase F.4 Method-A parent-preserving correction

DEVELOPMENT COMPLETE, FARM VALIDATION PENDING — F.4 begins from
`5382cfc1994b078c620b32c043938134c33ffa39` after accepted F.3 closure. It is
a detached, aggregate-only Method-A correction candidate over frozen `hgcer3`
F.3 models and the authoritative F.1 physical application population.

For each setting × canonical-t parent it applies F.3 response only inside the
reconstructed p99 support, uses raw `A = 1` outside support, then forms one
signed parent normalization `N = sum(b A) / sum(b)` and correction `C = A/N`.
This preserves each signed parent sum while redistributing it across the
existing application rows. Source labels and canonical-phi children are
diagnostic aggregates only and are never independently normalized.

F.4 persists no event identities/factors, templates, probabilities, weights,
or production objects. Local focused F.4 tests, frozen F.3/F.2/F.1 and
collector tests, syntax compilation, and F.4 profile JSON validation pass.
Farm validation must run only its detached analyzer and unchanged generic
collector against accepted F.1/F.3 evidence; no farm result is claimed here.
F.5 remains BLOCKED pending F.4 farm evidence.
