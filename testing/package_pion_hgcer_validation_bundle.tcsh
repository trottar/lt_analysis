#!/bin/tcsh -f
# Package only existing profile-declared validation artifacts.  This script does
# not run an analyzer, renderer, scheduler, or arbitrary caller-supplied command.

set repo = "/group/c-kaonlt/USERS/trottar/lt_analysis"
set python = ""
set bundle_commit = ""
set profile = ""
set artifact_dir = ""
set kinematic = ""
set phi = ""
set epsilon = ""
set output = ""
set immutable_paths = ()
set immutable_hashes = ()
set worktree_added = 0
set exit_code = 1

while ( $#argv > 0 )
    switch ( "$1" )
    case --python:
        if ( $#argv < 2 ) goto usage
        set python = "$2"
        shift argv
        shift argv
        breaksw
    case --bundle-commit:
        if ( $#argv < 2 ) goto usage
        set bundle_commit = "$2"
        shift argv
        shift argv
        breaksw
    case --profile:
        if ( $#argv < 2 ) goto usage
        set profile = "$2"
        shift argv
        shift argv
        breaksw
    case --artifact-dir:
        if ( $#argv < 2 ) goto usage
        set artifact_dir = "$2"
        shift argv
        shift argv
        breaksw
    case --kinematic:
        if ( $#argv < 2 ) goto usage
        set kinematic = "$2"
        shift argv
        shift argv
        breaksw
    case --phi:
        if ( $#argv < 2 ) goto usage
        set phi = "$2"
        shift argv
        shift argv
        breaksw
    case --epsilon:
        if ( $#argv < 2 ) goto usage
        set epsilon = "$2"
        shift argv
        shift argv
        breaksw
    case --output:
        if ( $#argv < 2 ) goto usage
        set output = "$2"
        shift argv
        shift argv
        breaksw
    case --immutable:
        if ( $#argv < 3 ) goto usage
        set immutable_paths = ( $immutable_paths "$2" )
        set immutable_hashes = ( $immutable_hashes "$3" )
        shift argv
        shift argv
        shift argv
        breaksw
    default:
        echo "unknown argument: $1"
        goto usage
    endsw
end

if ( "$python" == "" || "$bundle_commit" == "" || "$profile" == "" || "$artifact_dir" == "" || "$kinematic" == "" || "$output" == "" ) goto usage
if ( "$phi" == "" && "$epsilon" != "" ) goto usage
if ( "$phi" != "" && "$epsilon" == "" ) goto usage
if ( ! -x "$python" ) then
    echo "Python executable is unavailable: $python"
    exit 2
endif
if ( ! -d "$repo/.git" ) then
    echo "repository is unavailable: $repo"
    exit 2
endif
if ( ! -d "$artifact_dir" ) then
    echo "artifact directory is unavailable: $artifact_dir"
    exit 2
endif
if ( "$artifact_dir" != "/group/c-kaonlt/USERS/trottar/lt_analysis/OUTPUT/Analysis/KaonLT" ) then
    echo "artifact directory must be the canonical KaonLT artifact root"
    exit 2
endif

switch ( "$profile" )
case "":
case /*:
case ../*:
case */../*:
    echo "profile must be a repository-relative path without parent traversal"
    exit 2
endsw

if ( "$output" =~ /* ) then
    set zip = "$output"
else
    set zip = "/volatile/hallc/c-kaonlt/trottar/globus/$output"
endif
if ( -e "$zip" ) then
    echo "refusing to overwrite existing ZIP: $zip"
    exit 2
endif
switch ( "$zip" )
case /volatile/hallc/c-kaonlt/trottar/globus/*:
    breaksw
default:
    echo "ZIP output must be in the canonical Globus transfer directory"
    exit 2
endsw

set index = 1
while ( $index <= $#immutable_paths )
    set immutable_path = "$immutable_paths[$index]"
    set immutable_hash = "$immutable_hashes[$index]"
    set immutable_parent = `dirname "$immutable_path"`
    if ( "$immutable_parent" != "$artifact_dir" ) then
        echo "immutable artifact must be under the canonical artifact root: $immutable_path"
        exit 2
    endif
    if ( ! -f "$immutable_path" ) then
        echo "immutable artifact is unavailable: $immutable_path"
        exit 2
    endif
    set actual_hash = `sha256sum "$immutable_path" | awk '{print $1}'`
    if ( $status != 0 || "$actual_hash" != "$immutable_hash" ) then
        echo "immutable artifact SHA-256 mismatch: $immutable_path"
        exit 2
    endif
    @ index = $index + 1
end

cd "$repo"
set resolved_commit = `git rev-parse --verify "$bundle_commit^{commit}"`
if ( $status != 0 || "$resolved_commit" != "$bundle_commit" ) then
    echo "bundle commit must resolve to its supplied full 40-character identity: $bundle_commit"
    exit 2
endif
set timestamp = `date -u +%Y%m%d-%H%M%S`
set worktree = "/tmp/kaonlt-validation-bundle-$timestamp-$$"
if ( -e "$worktree" ) then
    echo "refusing existing temporary-worktree path: $worktree"
    exit 2
endif
onintr cleanup_failure
git worktree add --detach "$worktree" "$resolved_commit"
if ( $status != 0 ) exit 1
set worktree_added = 1
if ( ! -f "$worktree/$profile" ) then
    echo "profile is absent from detached bundle commit: $profile"
    goto cleanup_failure
endif

cd "$worktree"
set collector_args = ( --outdir "$artifact_dir" --kinematic "$kinematic" --profile "$profile" --output "$output" )
if ( "$phi" != "" ) then
    set collector_args = ( $collector_args --phi "$phi" --epsilon "$epsilon" )
endif
"$python" testing/collect_pion_hgcer_validation_bundle.py $collector_args
if ( $status != 0 ) goto cleanup_failure

set index = 1
while ( $index <= $#immutable_paths )
    set immutable_path = "$immutable_paths[$index]"
    set immutable_hash = "$immutable_hashes[$index]"
    set actual_hash = `sha256sum "$immutable_path" | awk '{print $1}'`
    if ( $status != 0 || "$actual_hash" != "$immutable_hash" ) then
        echo "immutable artifact changed or is unreadable: $immutable_path"
        goto cleanup_failure
    endif
    @ index = $index + 1
end

unzip -t "$zip"
if ( $status != 0 ) goto cleanup_failure
unzip -p "$zip" manifest.json | "$python" -m json.tool
if ( $status != 0 ) goto cleanup_failure
sha256sum "$zip"
echo "validation ZIP: $zip"
set exit_code = 0
goto cleanup

cleanup_failure:
set exit_code = 1

cleanup:
if ( $worktree_added == 1 ) then
    cd "$repo"
    git worktree remove --force "$worktree"
    set cleanup_status = $status
    git worktree prune
    if ( $exit_code == 0 && $cleanup_status != 0 ) then
        echo "temporary detached worktree could not be removed: $worktree"
        set exit_code = 1
    endif
endif
exit $exit_code

usage:
echo "Usage: tcsh package_pion_hgcer_validation_bundle.tcsh --python PATH --bundle-commit SHA1 --profile REPO_RELATIVE_PATH --artifact-dir PATH --kinematic TOKEN --output ZIP_OR_GLOBUS_PATH [--phi VALUE --epsilon VALUE] [--immutable PATH SHA256] ..."
exit 2
