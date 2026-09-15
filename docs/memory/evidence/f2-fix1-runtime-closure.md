# F.2.Fix.1 accepted farm bundle and basis decision

## Evidence identity

- Bundle: `KaonLT_PhaseF2v1_validation_Q4p4W2p74.zip`.
- ZIP SHA-256: `9081b2741677e1ca4976eb0cabf3dab38c43ec59a9ec479bd1897707832177c6`.
- Manifest: `pion_hgcer_validation_bundle/v4`, generated
  `2026-09-15T03:16:02.494168Z`, `complete = true` with no errors or
  unexpected committed files.
- Farm analyzer/collector HEAD: `8e919fc618cea900227db5090d65da728c3aa555`;
  required source `170e6fae3d2fed1949fc6932b8eac9ad83e3e01c` is an ancestor.
- The farm `git status --short` recorded unrelated normal-analysis edits and
  untracked Kaon analysis outputs. They are evidence provenance, not F.2
  inputs or an F.2 source change; F.2 reads only the packaged JSON artifacts.

## Artifact continuity

- F.2 JSON SHA-256:
  `87e0a16870ea34a6611702e1bd4fc02c4dc6a8daa132e75b725a1b782a9a31da`;
  artifact fingerprint:
  `19bb6b55f622d11c63dc4a69932d8dc2403f5a1a9c3546b599a0661c25f317ab`.
- F.2 representation fingerprint:
  `8ffd44fca4ae5c0a53a057a44e1a1ad5a62363ce56c39bc393a4a2aecfa233d3`.
- F.2 PDF SHA-256:
  `5641515bd0eefc45bf3d12c4dff2620277323c6bf7b43f83d2df6065af821b55`.
- The five source F.1 JSON SHA-256 values are preserved by the manifest:
  Left-lowe `f16c5e89f62e848ec221335bd4916265757040af9970b005d8f67cc800e0077d`;
  Left-highe `203f4c76f1a251e3e8f231fa3a5e50c9a4fa337420efffffd803205c6d7ea218`;
  Center-lowe `1593e22b55382b4a9e831d3a1114584e2e3057fcbc4aeeea4aa74c948edf39f8`;
  Center-highe `5de64b850735ebe70040a703bac999bdd1ac84dc821aba8e1a21ec86ae4b9db3`;
  Right-highe `77d98006ff81e772c466bf8d10ec83088509a0b8a4d430bded1172bc118bc15f`.

## Accepted F.2 result

The owner accepted F.2 as `CLOSED / RUNTIME VALIDATED`. All 15 canonical
setting × t parents had valid, non-sparse `hgcer3` response and support gates.
`delta_only` and `track3` failed fixed information retention; `full5_reference`
remains diagnostic-only. The accepted reduced basis is frozen for F.3 only:

```text
hgcer3 = (SHMS_delta, P_hgcer_xAtCer, P_hgcer_yAtCer)
```

For `hgcer3`, median/max AUC loss were `-0.0007784076762388237` and
`0.020042411293277107`; median/max balanced-log-loss penalty were
`0.0010416346895916195` and `0.011262387818522734`; maximum application OOD
fraction was `0.0588235294117647`. The frozen F.2 thresholds were not retuned.

This acceptance changes neither the historical F.2 artifact field
`basis_frozen = false` nor the F.2 runtime identity. F.3 records its own
`basis_frozen = true` only because of this explicit human decision.
