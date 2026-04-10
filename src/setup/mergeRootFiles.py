
#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-11-27 11:40:10 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import sys, os, subprocess
from pathlib import Path
import ROOT

input_root_path = Path(sys.argv[1]).expanduser()
output_root_path = Path(sys.argv[2]).expanduser()
inp_file_name = sys.argv[3]
inp_tree_names = sys.argv[4]
output_file_name = sys.argv[5]
string_run_nums = sys.argv[6]
particle = sys.argv[7]
err_fout = sys.argv[8]

###############################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import progres bar
from ltsep import Misc

##################################################################################################################################################
# Importing utility functions

sys.path.append("../utility")
from utility import open_root_file, process_lines

###############################################################################################################################################

# Overwrite error file if already exists
with open(err_fout, 'w') as f:
    f.write("Bad runs for {}...\n\n".format(output_file_name))

def log_bad_runs(inp_root_file, err_fout, warning):
    print(warning)
    with open(err_fout, 'a') as f:
        f.write(warning+'\n')

def resolve_absolute_path(path_obj):
    return path_obj.resolve(strict=False)


def absolute_path_has_cache(path_obj):
    return "cache" in str(path_obj).lower()


def cache_path_to_mss(path_obj):
    path_text = str(path_obj)
    kaonlt_cache_index = path_text.find("/cache/hallc/kaonlt")
    if kaonlt_cache_index != -1:
        suffix = path_text[kaonlt_cache_index + len("/cache/hallc/kaonlt"):]
        return Path("/mss/hallc/kaonlt" + suffix)

    generic_cache_index = path_text.find("/cache/")
    if generic_cache_index != -1:
        suffix = path_text[generic_cache_index + len("/cache"):]
        return Path("/mss" + suffix)

    return None


def run_jcache_in_batches(mss_files, batch_size=50):
    unique_files = []
    seen = set()
    for mss_file in mss_files:
        mss_text = str(mss_file)
        if mss_text not in seen:
            unique_files.append(mss_text)
            seen.add(mss_text)

    if not unique_files:
        return

    total_batches = (len(unique_files) + batch_size - 1) // batch_size
    print("\nCache-backed skim source detected; requesting jcache staging for combine inputs...")
    for batch_index in range(total_batches):
        start = batch_index * batch_size
        batch = unique_files[start:start + batch_size]
        command = ["jcache", "get"] + batch
        print("  Batch {}/{}: {}".format(batch_index + 1, total_batches, " ".join(command)))
        try:
            subprocess.run(command, check=True)
        except FileNotFoundError:
            print("ERROR: jcache command not found. Cannot stage cache-backed skim inputs.")
            sys.exit(1)
        except subprocess.CalledProcessError as exc:
            print("ERROR: jcache staging failed for combine inputs (exit code {}).".format(exc.returncode))
            sys.exit(exc.returncode if exc.returncode != 0 else 1)


input_root_abs = resolve_absolute_path(input_root_path)
output_root_abs = resolve_absolute_path(output_root_path)
print("\nCombining ROOT files from: {}".format(input_root_path))
print("Absolute skim source path: {}".format(input_root_abs))
print("Writing merged ROOT file to: {}".format(output_root_abs))

arr_run_nums = [int(x) for x in string_run_nums.split()]
input_root_files = [
    input_root_path / (particle + "_" + str(run_num) + inp_file_name + ".root")
    for run_num in arr_run_nums
]
absolute_input_root_files = [
    input_root_abs / (particle + "_" + str(run_num) + inp_file_name + ".root")
    for run_num in arr_run_nums
]

if absolute_path_has_cache(input_root_abs):
    mss_stage_files = []
    for path_obj in absolute_input_root_files:
        if path_obj.exists():
            print("Cache file already present, skipping jcache: {}".format(path_obj))
            continue
        mss_path = cache_path_to_mss(path_obj)
        if mss_path is None:
            print("WARNING: Could not map cache-backed input to MSS for jcache staging: {}".format(path_obj))
            continue
        mss_stage_files.append(mss_path)
    run_jcache_in_batches(mss_stage_files)

os.makedirs(output_root_path, exist_ok=True)
out_root_file = os.path.join(output_root_path, output_file_name + ".root")

outfile = ROOT.TFile(out_root_file, "RECREATE")
if not outfile.IsOpen():
    print("ERROR: Output file {} cannot be opened. Exiting the function.".format(outfile.GetName()))
    sys.exit(1)

for tree in inp_tree_names.split():

    chain = ROOT.TChain(tree)

    for i,n in enumerate(arr_run_nums):
        # Progress bar
        if len(arr_run_nums) > 1:
            Misc.progressBar(i, len(arr_run_nums)-1,bar_length=25)
        else:
            Misc.progressBar(len(arr_run_nums), len(arr_run_nums),bar_length=25)
        inp_root_file = os.path.join(input_root_path, particle + "_" + str(n) + inp_file_name + ".root")
        if not os.path.isfile(inp_root_file):
            warning = "WARNING: File {} not found. Removing...".format(inp_root_file)
            log_bad_runs(inp_root_file, err_fout, warning)
            continue
        tempfile = open_root_file(inp_root_file)
        if tempfile == None or not tempfile.IsOpen() or tempfile.TestBit(ROOT.TFile.kRecovered):
            warning = "WARNING: File {} not found or not opened or corrupted. Removing...".format(inp_root_file)
            log_bad_runs(inp_root_file, err_fout, warning)
            if os.path.exists(inp_root_file):        
                os.remove(inp_root_file)            
            continue
        # Get the tree from the temporary file using the tree_name
        tree_temp = tempfile.Get(tree)
        # Check if the tree exists
        if tree_temp:
            # Get the number of entries in the tree
            num_entries = tree_temp.GetEntries()
            if num_entries == 0:
                warning = "WARNING: Tree {} in file {} is empty. Removing...".format(tree, inp_root_file)
                log_bad_runs(inp_root_file, err_fout, warning)
                continue
        #print("Adding {}...".format(inp_root_file))
        chain.Add(inp_root_file)

    if os.path.exists(inp_root_file):
        
        if chain.GetEntries() == 0:
            warning = "WARNING: No entries found for tree {}.  Removing...".format(tree)        
            log_bad_runs(inp_root_file, err_fout, warning)
            continue

        outfile.cd()
        chain.Write(tree, ROOT.TObject.kWriteDelete)
        print("\n\tTree {} added to {}.root".format(tree,output_file_name))
    else:
        continue
    
outfile.Close()

# Remove superfluous statements in error file
for n in arr_run_nums:
    process_lines(str(n), err_fout)
