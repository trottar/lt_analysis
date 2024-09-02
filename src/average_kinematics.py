#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-09-01 20:12:20 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import subprocess
import sys

# Read command-line arguments
inp_pid = sys.argv[1]
inp_pol = sys.argv[2]
inp_Q2 = sys.argv[3]
inp_W = sys.argv[4]
inp_loeps = sys.argv[5]
inp_hieps = sys.argv[6]

# Command to execute the Fortran program
cmd = ['./average_kinematics']

# Run the Fortran program and interact with it
process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

# Wait for the prompt and send input
stdout, stderr = process.communicate(input=f"{inp_pid}\n{inp_pol}\n{inp_Q2}\n{inp_W}\n{inp_loeps}\n{inp_hieps}\n")

# Print the output for verification
print(stdout)
if stderr:
    print("Errors:", stderr)
