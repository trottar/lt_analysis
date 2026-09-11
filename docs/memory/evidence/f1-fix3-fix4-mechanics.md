# F.1.Fix.3 and Fix.4 farm mechanics

## Status

CLOSED / RUNTIME VALIDATED for the two named mechanical regressions only.
This is not scientific acceptance of F.1 v1 and not a F.1.Fix.5 farm result.

## Evidence

[FARM_EVIDENCE] [HANDOFF]

F.1 lineage through the repaired mechanics includes:

    d656e15761970d7d612bb028d2746d077795e9ad  F.1.Fix.2
    8b1ad5b735e1b5cc95a5d876526d321296cec16b  F.1.Fix.3
    eb253046ff6e23ec94c8638315c4bc712fa4f992  F.1.Fix.4

Fix.3 followed a farm ValueError stating that truth value of an array with more
than one element is ambiguous. NumPy-like edge arrays had been used in a
truth-value fallback expression; safe edge materialization repaired it.

Fix.4 followed parent_child_parity_mismatch:allcuts. Parent values were native
Python bool while child cache values could be NumPy boolean scalars used with
identity-style comparison. The F.1 read boundary now detaches NumPy-like
scalars to native Python values.

## Later five-setting mechanical gate

All five settings produced F.1 artifacts with status=available,
diagnostic_stage=complete, nonempty records, populated fingerprints, and five
F.1 pages. Focused source checks recorded:

- full-background tests: 90 PASS;
- F.1 contract tests: 10 PASS;
- Phase-F runtime tests: 3 PASS;
- collector tests: 14 PASS;
- py_compile PASS;
- git diff --check PASS.

The gate exposed the population-ownership defect below; it must not be read as
F.1 v1 scientific acceptance or as evidence for F.1.Fix.5.

