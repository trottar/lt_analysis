#!/bin/tcsh -f
#
# KaonLT F.6.2.Fix.4 final validation-bundle creation
# Date: 2026-09-17
#
# This script ONLY packages existing accepted/rerendered artifacts.
# It does not rerun F.6.2 science and does not regenerate the F.6.2 JSON.
#

umask 002
setenv PYTHONDONTWRITEBYTECODE 1
setenv PYTHONUNBUFFERED 1

# ----------------------------------------------------------------------
# Frozen paths / identities
# ----------------------------------------------------------------------

set REPO_DIR = "/group/c-kaonlt/USERS/trottar/lt_analysis"
set OUTDIR = "${REPO_DIR}/OUTPUT/Analysis/KaonLT"
set GLOBUS_DIR = "/volatile/hallc/c-kaonlt/trottar/globus"

set KIN = "Q4p4W2p74"
set COLLECTOR = "${REPO_DIR}/testing/collect_pion_hgcer_validation_bundle.py"
set PROFILE = "${REPO_DIR}/testing/pion_hgcer_validation_bundle_profile_f6_2.json"

set EXPECTED_HEAD = "c3252f1ec43dcca0357aebe8c6ca75d7d198f99a"
set EXPECTED_SOURCE = "b929815517643edaa75949f7e61ce46c6e8f0d63"

set EXPECTED_JSON_SHA = "5fb52310b44c4fbba66bbbf868c0c7ee8894992a8f06f2d0bd209d1608310bb1"
set EXPECTED_VALIDATION_FP = "7edc73fce20ad7dc7622b8c23a7ba7e8986e367ace605884e010945595370b3b"
set EXPECTED_ARTIFACT_FP = "ee713b70de898bad8fa61164cbc8a54886af1ea1eb1e712ed95de17df8890cd0"
set OLD_PRE_FIX4_PDF_SHA = "4d07c2989c217d18d088b181390e9478fae51e4168dba2d64faea690b40f0f85"

set F62_JSON = "${OUTDIR}/${KIN}_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.json"
set F62_PDF = "${OUTDIR}/${KIN}_kaon_pion-background_hgcer_method-a-acceptance-refinement-validation.pdf"

set ERRMSG = ""

# ----------------------------------------------------------------------
# Stage 1: repository / executable inputs
# ----------------------------------------------------------------------

echo ""
echo "======================================================================"
echo "Stage 1/6: source identity"
echo "======================================================================"

if ( ! -d "${REPO_DIR}" ) then
    set ERRMSG = "repository directory missing: ${REPO_DIR}"
    goto FAIL
endif

cd "${REPO_DIR}"
if ( $status != 0 ) then
    set ERRMSG = "could not cd to repository"
    goto FAIL
endif

if ( ! -f "${COLLECTOR}" ) then
    set ERRMSG = "collector missing: ${COLLECTOR}"
    goto FAIL
endif

if ( ! -f "${PROFILE}" ) then
    set ERRMSG = "F.6.2 bundle profile missing: ${PROFILE}"
    goto FAIL
endif

echo "Repository status (record only; unrelated farm worktree files are allowed):"
git status --short --branch
if ( $status != 0 ) then
    set ERRMSG = "git status failed"
    goto FAIL
endif

set HEAD_NOW = `git rev-parse HEAD`
if ( $status != 0 ) then
    set ERRMSG = "git rev-parse HEAD failed"
    goto FAIL
endif

echo "Local HEAD: ${HEAD_NOW}"
if ( "${HEAD_NOW}" != "${EXPECTED_HEAD}" ) then
    set ERRMSG = "wrong local HEAD; expected ${EXPECTED_HEAD}"
    goto FAIL
endif

# Do not reject the already-known unrelated model/model-file farm changes.
# Do reject uncommitted edits to code/configuration that will execute NOW.
set CRITICAL_DIRTY = ( `git status --porcelain -- src/cuts/pion_hgcer_method_a_acceptance_refinement_validation.py testing/analyze_pion_hgcer_method_a_acceptance_refinement_validation.py testing/collect_pion_hgcer_validation_bundle.py testing/pion_hgcer_validation_bundle_profile_f6_2.json` )
if ( $#CRITICAL_DIRTY != 0 ) then
    echo ""
    echo "Critical validation files have uncommitted changes:"
    git status --short -- src/cuts/pion_hgcer_method_a_acceptance_refinement_validation.py testing/analyze_pion_hgcer_method_a_acceptance_refinement_validation.py testing/collect_pion_hgcer_validation_bundle.py testing/pion_hgcer_validation_bundle_profile_f6_2.json
    set ERRMSG = "refusing to bundle with uncommitted validation/collector changes"
    goto FAIL
endif

set PROFILE_SOURCE = `python -c 'import json,sys; print(json.load(open(sys.argv[1], encoding="utf-8"))["source_identity"]["required_analysis_commit"])' "${PROFILE}"`
if ( $status != 0 ) then
    set ERRMSG = "could not read required_analysis_commit from profile"
    goto FAIL
endif

echo "Profile required analysis source: ${PROFILE_SOURCE}"
if ( "${PROFILE_SOURCE}" != "${EXPECTED_SOURCE}" ) then
    set ERRMSG = "F.6.2 profile source pin mismatch"
    goto FAIL
endif

# ----------------------------------------------------------------------
# Stage 2: accepted artifact identity
# ----------------------------------------------------------------------

echo ""
echo "======================================================================"
echo "Stage 2/6: accepted F.6.2 scientific artifact"
echo "======================================================================"

if ( ! -d "${OUTDIR}" ) then
    set ERRMSG = "artifact directory missing: ${OUTDIR}"
    goto FAIL
endif

if ( ! -d "${GLOBUS_DIR}" ) then
    set ERRMSG = "Globus transfer directory missing: ${GLOBUS_DIR}"
    goto FAIL
endif

if ( ! -f "${F62_JSON}" ) then
    set ERRMSG = "F.6.2 JSON missing: ${F62_JSON}"
    goto FAIL
endif

if ( ! -s "${F62_PDF}" ) then
    set ERRMSG = "F.6.2 PDF missing or empty: ${F62_PDF}"
    goto FAIL
endif

set JSON_SHA = `sha256sum "${F62_JSON}" | awk '{print $1}'`
if ( $status != 0 ) then
    set ERRMSG = "sha256sum failed for F.6.2 JSON"
    goto FAIL
endif

echo "F.6.2 JSON SHA-256: ${JSON_SHA}"
if ( "${JSON_SHA}" != "${EXPECTED_JSON_SHA}" ) then
    set ERRMSG = "F.6.2 scientific JSON hash changed; refusing to bundle"
    goto FAIL
endif

python -c 'import json,sys; a=json.load(open(sys.argv[1], encoding="utf-8")); v=a["validation"]; assert v["fingerprint"] == sys.argv[2], v["fingerprint"]; assert a["artifact_fingerprint"] == sys.argv[3], a["artifact_fingerprint"]; print("F.6.2 validation fingerprint:", v["fingerprint"]); print("F.6.2 artifact fingerprint:", a["artifact_fingerprint"]); print("Scientific identity PASS")' "${F62_JSON}" "${EXPECTED_VALIDATION_FP}" "${EXPECTED_ARTIFACT_FP}"
if ( $status != 0 ) then
    set ERRMSG = "F.6.2 scientific fingerprint check failed"
    goto FAIL
endif

# ----------------------------------------------------------------------
# Stage 3: rerendered PDF sanity
# ----------------------------------------------------------------------

echo ""
echo "======================================================================"
echo "Stage 3/6: Fix.4 PDF sanity"
echo "======================================================================"

set PDF_SHA = `sha256sum "${F62_PDF}" | awk '{print $1}'`
if ( $status != 0 ) then
    set ERRMSG = "sha256sum failed for F.6.2 PDF"
    goto FAIL
endif

echo "Current F.6.2 PDF SHA-256: ${PDF_SHA}"

if ( "${PDF_SHA}" == "${OLD_PRE_FIX4_PDF_SHA}" ) then
    set ERRMSG = "PDF is still the pre-Fix.4 artifact; rerender did not replace it"
    goto FAIL
endif

set PAGE_COUNT = `python -c 'import sys; d=open(sys.argv[1],"rb").read(); print(d.count(b"/Type /Page") - d.count(b"/Type /Pages"))' "${F62_PDF}"`
if ( $status != 0 ) then
    set ERRMSG = "could not determine PDF page count"
    goto FAIL
endif

echo "F.6.2 PDF page count: ${PAGE_COUNT}"
if ( "${PAGE_COUNT}" != "132" ) then
    set ERRMSG = "unexpected F.6.2 PDF page count; expected 132"
    goto FAIL
endif

# ----------------------------------------------------------------------
# Stage 4: create one unique fresh bundle
# ----------------------------------------------------------------------

echo ""
echo "======================================================================"
echo "Stage 4/6: create validation bundle"
echo "======================================================================"

set STAMP = `date +%Y%m%d-%H%M%S`
set BUNDLE = "KaonLT_PhaseF6_2_Fix4_validation_${KIN}_${STAMP}.zip"
set BUNDLE_PATH = "${GLOBUS_DIR}/${BUNDLE}"

if ( -e "${BUNDLE_PATH}" ) then
    set ERRMSG = "bundle path already exists: ${BUNDLE_PATH}"
    goto FAIL
endif

echo "Bundle basename: ${BUNDLE}"
echo "Collector destination: ${BUNDLE_PATH}"

python "${COLLECTOR}" \
    --outdir "${OUTDIR}" \
    --kinematic "${KIN}" \
    --profile "${PROFILE}" \
    --output "${BUNDLE}"

if ( $status != 0 ) then
    set ERRMSG = "collector returned nonzero; do NOT transfer the generated bundle"
    goto FAIL
endif

if ( ! -s "${BUNDLE_PATH}" ) then
    set ERRMSG = "collector reported success but bundle is missing/empty"
    goto FAIL
endif

# ----------------------------------------------------------------------
# Stage 5: ZIP + manifest gates
# ----------------------------------------------------------------------

echo ""
echo "======================================================================"
echo "Stage 5/6: ZIP and manifest validation"
echo "======================================================================"

unzip -t "${BUNDLE_PATH}"
if ( $status != 0 ) then
    set ERRMSG = "ZIP integrity test failed"
    goto FAIL
endif

python -c 'import json,sys,zipfile; z=zipfile.ZipFile(sys.argv[1]); m=json.loads(z.read("manifest.json")); expected_settings=[("Left","lowe"),("Left","highe"),("Center","lowe"),("Center","highe"),("Right","highe")]; assert m.get("complete") is True, m.get("errors"); assert m.get("errors") == [], m.get("errors"); assert m.get("git_head") == sys.argv[2], m.get("git_head"); assert m.get("required_analysis_commit") == sys.argv[3], m.get("required_analysis_commit"); assert m.get("required_analysis_commit_is_ancestor") is True; assert m.get("unexpected_committed_files_after_required_analysis_commit") == [], m.get("unexpected_committed_files_after_required_analysis_commit"); assert [(x["phi"],x["epsilon"]) for x in m.get("requested_settings",[])] == expected_settings, m.get("requested_settings"); assert len(m.get("global_artifacts",{})) == 6; assert len(m.get("settings",[])) == 5; j=m["global_artifacts"]["f6_2_acceptance_refinement_validation_json"]; p=m["global_artifacts"]["f6_2_acceptance_refinement_validation_pdf"]; assert j.get("status") == "exists" and j.get("json_status") == "valid"; assert j.get("sha256") == sys.argv[4], j.get("sha256"); assert p.get("status") == "exists"; assert p.get("sha256") == sys.argv[5], p.get("sha256"); assert all(s["artifacts"]["method_a_acceptance_contract"].get("status") == "exists" and s["artifacts"]["method_a_acceptance_contract"].get("json_status") == "valid" for s in m["settings"]); print("Manifest complete:", m["complete"]); print("Manifest errors:", m["errors"]); print("Bundle HEAD:", m["git_head"]); print("Reviewed source:", m["required_analysis_commit"]); print("Manifest gate PASS")' "${BUNDLE_PATH}" "${EXPECTED_HEAD}" "${EXPECTED_SOURCE}" "${EXPECTED_JSON_SHA}" "${PDF_SHA}"

if ( $status != 0 ) then
    set ERRMSG = "bundle manifest gate failed"
    goto FAIL
endif

# ----------------------------------------------------------------------
# Stage 6: final transfer identity
# ----------------------------------------------------------------------

echo ""
echo "======================================================================"
echo "Stage 6/6: final bundle identity"
echo "======================================================================"

sha256sum "${BUNDLE_PATH}"
if ( $status != 0 ) then
    set ERRMSG = "final bundle sha256sum failed"
    goto FAIL
endif

ls -lh "${BUNDLE_PATH}"

echo ""
echo "======================================================================"
echo "PASS"
echo "======================================================================"
echo "Transfer/upload exactly this fresh ZIP:"
echo "${BUNDLE_PATH}"
echo ""
echo "Do not rename it to an older F.6.2 bundle filename."
exit 0

FAIL:
echo ""
echo "======================================================================"
echo "FAILED — bundle not accepted for transfer"
echo "======================================================================"
echo "ERROR: ${ERRMSG}"
echo ""
if ( $?BUNDLE_PATH ) then
    if ( -e "${BUNDLE_PATH}" ) then
        echo "A file exists at:"
        echo "${BUNDLE_PATH}"
        echo "Do NOT transfer it unless the failure is reviewed."
    endif
endif
exit 1
