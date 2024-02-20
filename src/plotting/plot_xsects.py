#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-02-20 00:29:26 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import pandas as pd
import root_pandas as rpd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import re
import sys, os

##################################################################################################################################################
# Check the number of arguments provided to the script

if len(sys.argv)-1!=10:
    print("!!!!! ERROR !!!!!\n Expected 10 arguments\n Usage is with - ParticleType POL Q2 W LOEPS HIEPS NumtBins NumPhiBins KIN OutUnsepxsectsFilename\n!!!!! ERROR !!!!!")
    sys.exit(1)

###############################################################################################################################################

ParticleType = sys.argv[1]
POL = sys.argv[2]

Q2 = sys.argv[3]
W = sys.argv[4]

LOEPS = sys.argv[5]
HIEPS = sys.argv[6]

NumtBins = int(sys.argv[7])
NumPhiBins = int(sys.argv[8])

kinematics = sys.argv[9]
OutFilename = sys.argv[10]

if int(POL) == 1:
    pol_str = "pl"
elif int(POL) == -1:
    pol_str = "mn"
else:
    print("ERROR: Invalid polarity...must be +1 or -1")
    sys.exit(2)
    
###############################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_Prod")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

foutname = OUTPATH+"/" + OutFilename + ".root"
fouttxt  = OUTPATH+"/" + OutFilename + ".txt"
outputpdf  = OUTPATH+"/" + OutFilename + ".pdf"

################################################################################################################################################
'''
Import separated xsects model
'''

sys.path.append("../models")
if pol_str == "pl" and ParticleType == "kaon":
    from sep_xsect_kaon_pl import import_model

###############################################################################################################################################
'''
Import model parameterization
'''

param_file = '{}/src/{}/parameters/par.{}_Q{}W{}.dat'.format(LTANAPATH, ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""))

param_arr = []
with open(param_file, 'r') as f:
    for i, line in enumerate(f):
        columns = line.split()
        param_arr.append(str(columns[0]))    

###############################################################################################################################################

def are_within_tolerance(num1, num2, tolerance=0.1):
    return abs(num1 - num2) <= tolerance

def file_to_df(f_name, columns, filter_conditions=None):
    '''
    Read in file and convert to dataframe with custom column names
    Optionally filter the dataframe based on specified conditions
    '''
    lineskip = False
    
    # Open the file for reading
    with open(f_name, 'r') as file:
        lines = file.readlines()
        if '1.000000\n' in lines:
            lineskip = True

    if lineskip:
        df = pd.read_csv(f_name, header=None, sep=' ', skiprows=1, skipfooter=1, engine='python')
    else:
        df = pd.read_csv(f_name, header=None, sep=' ')    
    df.columns = columns
    
    # Filter the dataframe based on specified conditions
    if filter_conditions:
        for column, value in filter_conditions.items():
            df = df[df[column] == value]
    
    return df

################################################################################################################################################

def fix_spacing(f_name):
    '''
    Fortran created files are bad with spacing. This fixes it.
    '''

    # Open the file for reading
    with open(f_name, 'r') as file:

        # Read the lines of the file and split based on whitespace
        lines = file.readlines()
        lines = [re.split(r'\s+', line.strip()) for line in lines]

        # Join the split lines with a single space
        lines = [' '.join(line) for line in lines]

        # Write the lines back to the file
        with open(f_name, 'w') as output:
            output.write('\n'.join(lines))

# Fix file spacing to work in pandas
fix_spacing(LTANAPATH+"/src/{}/averages/avek.Q{}W{}.dat".format(ParticleType, Q2.replace("p",""), W.replace("p","")))
fix_spacing(LTANAPATH+"/src/{}/averages/aver.{}_Q{}W{}_{:.0f}.dat".format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100))
fix_spacing(LTANAPATH+"/src/{}/averages/aver.{}_Q{}W{}_{:.0f}.dat".format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100))
fix_spacing(LTANAPATH+"/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat".format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100))
fix_spacing(LTANAPATH+"/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat".format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100))
fix_spacing(LTANAPATH+"/src/{}/xsects/x_sep.{}_Q{}W{}.dat".format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p","")))
################################################################################################################################################
# Read in files and convert to dataframes

file_df_dict = {}

setting_file = LTANAPATH+"/src/{}/list.settings".format(ParticleType)
file_df_dict['setting_df'] = file_to_df(setting_file, ['POL', 'Q2', 'W', 'EPSVAL', 'thpq', 'TMIN', 'TMAX', 'NumtBins'],\
                                        filter_conditions={'Q2': float(Q2.replace("p",".")), 'W': float(W.replace("p","."))})

tmin = file_df_dict['setting_df'].iloc[0]['TMIN']
tmax = file_df_dict['setting_df'].iloc[0]['TMAX']

try:
    with open("{}/src/{}/t_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p","")), "r") as file:
        # Read all lines from the file into a list
        all_lines = file.readlines()
        # Check if the file has at least two lines
        if len(all_lines) >= 2:
            # Extract the second line and remove leading/trailing whitespace
            t_bins = all_lines[1].split("\t")
            del t_bins[0]
            t_bins = np.array([float(element) for element in t_bins])
except FileNotFoundError:
    print("{} not found...".format("{}/src/{}/t_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p",""))))
except IOError:
    print("Error reading {}...".format("{}/src/{}/t_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p",""))))    

t_bin_centers = (t_bins[:-1] + t_bins[1:]) / 2

try:
    with open("{}/src/{}/phi_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p","")), "r") as file:
        # Read all lines from the file into a list
        all_lines = file.readlines()
        # Check if the file has at least two lines
        if len(all_lines) >= 2:
            # Extract the second line and remove leading/trailing whitespace
            phi_bins = all_lines[1].split("\t")
            del phi_bins[0]
            phi_bins = np.array([float(element) for element in phi_bins])
except FileNotFoundError:
    print("{} not found...".format("{}/src/{}/phi_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p",""))))
except IOError:
    print("Error reading {}...".format("{}/src/{}/phi_bin_interval_Q{}W{}".format(LTANAPATH, ParticleType, Q2.replace("p",""), W.replace("p",""))))    

#phi_bin_centers = (phi_bins[:-1] + phi_bins[1:]) / 2
phi_bin_centers = phi_bins
    
for i,row in file_df_dict['setting_df'].iterrows():
    if row['Q2'] == float(Q2.replace("p",".")):
        file_df_dict['beam_file'] = file_to_df(LTANAPATH+"/src/{}/beam/Eb_KLT.dat".format(ParticleType), ['ebeam', 'Q2', 'W', 'EPSVAL'])
        file_df_dict['avek_file'] = file_to_df(LTANAPATH+"/src/{}/averages/avek.Q{}W{}.dat".format(ParticleType, Q2.replace("p",""), W.replace("p","")) \
                                               , ['W', 'dW', 'Q2', 'dQ2', 't', 'dt', 'th_pos', "tbin"])        
        
        if row['EPSVAL'] == float(LOEPS):
            file_df_dict['aver_loeps'] = file_to_df( \
                                                     LTANAPATH+"/src/{}/averages/aver.{}_Q{}W{}_{:.0f}.dat" \
                                                     .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100) \
                                                     , ['ratio', 'dratio', 'phibin', 'tbin']).sort_values(by='tbin')
            if row['thpq'] < 0.0:
                file_df_dict['kindata_loeps_{}'.format('right')] = file_to_df( \
                                                                               LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_{}.dat" \
                                                                               .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                       float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                               , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] > 0.0:
                file_df_dict['kindata_loeps_{}'.format('left')] = file_to_df( \
                                                                              LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+{}.dat" \
                                                                              .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                      float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                              , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] == 0.0:
                file_df_dict['kindata_loeps_{}'.format('center')] = file_to_df( \
                                                                                LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+0000.dat" \
                                                                                .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                        float(LOEPS)*100) \
                                                                                , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            file_df_dict['unsep_file_loeps'] = file_to_df( \
                                                            LTANAPATH+"/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat" \
                                                            .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100) \
                                                            , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 't', 't_min', 'W', 'Q2']).sort_values(by='t')

        if row['EPSVAL'] == float(HIEPS):
            file_df_dict['aver_hieps'] = file_to_df( \
                                                     LTANAPATH+"/src/{}/averages/aver.{}_Q{}W{}_{:.0f}.dat" \
                                                     .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100) \
                                                     , ['ratio', 'dratio', 'phibin', 'tbin']).sort_values(by='tbin')
            if row['thpq'] < 0.0:
                file_df_dict['kindata_hieps_{}'.format('right')] = file_to_df( \
                                                                               LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_{}.dat" \
                                                                               .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                       float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                               , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] > 0.0:
                file_df_dict['kindata_hieps_{}'.format('left')] = file_to_df( \
                                                                              LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+{}.dat" \
                                                                              .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                      float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                              , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            if row['thpq'] == 0.0:
                file_df_dict['kindata_hieps_{}'.format('center')] = file_to_df( \
                                                                                LTANAPATH+"/src/{}/kindata/kindata.{}_Q{}W{}_{:.0f}_+0000.dat" \
                                                                                .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                        float(HIEPS)*100) \
                                                                                , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
            file_df_dict['unsep_file_hieps'] = file_to_df( \
                                                            LTANAPATH+"/src/{}/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat" \
                                                            .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100) \
                                                            , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 't', 't_min', 'W', 'Q2']).sort_values(by='t')
        file_df_dict['sep_file'] = file_to_df( \
                                               LTANAPATH+"/src/{}/xsects/x_sep.{}_Q{}W{}.dat" \
                                               .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p","")) \
                                               , ['sigL', 'dsigL', 'sigT', 'dsigT', 'sigLT', 'dsigLT', 'sigTT', 'dsigTT', 'chisq', 't', 'tm', 'W', 'Q2', 'th_cm'])
            
################################################################################################################################################

# Create a PdfPages object to manage the PDF file
with PdfPages(outputpdf) as pdf:
    
    # Create a figure and axis objects
    fig, axes = plt.subplots(NumtBins, 1, figsize=(8, 6 * NumtBins), sharex=True)

    # Define markers and colors
    markers = ['o', 's'] # 'o'->circle, 's'->square
    colors = ['blue', 'orange']

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t = {:.2f}".format(t_bin_centers[k]))

        for i, df_key in enumerate(['aver_loeps', 'aver_hieps']):
            df = file_df_dict[df_key]
            mask = (df['tbin'] == (k+1))
            ax.errorbar(phi_bin_centers[df['phibin'][mask]], df['ratio'][mask], yerr=df['dratio'][mask], marker=markers[i], linestyle='None', label=df_key, color=colors[i])

        ax.axhline(1.0, color='gray', linestyle='--')

        ax.set_xlabel('$\phi$')
        ax.set_ylabel('Ratio')
        ax.set_ylim(0.0, 2.0)
        ax.set_xlim(0, 360)

        ax.legend()
        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')

    # Loop through t bins and plot data
    for k in range(NumtBins):

        # Create a figure and axis objects
        fig, axes = plt.subplots(1, 1, figsize=(8, 6), sharex=True)

        ax = axes
        ax.set_title("t = {:.2f}".format(t_bin_centers[k]))

        for i, df_key in enumerate(['aver_loeps', 'aver_hieps']):
            df = file_df_dict[df_key]
            mask = (df['tbin'] == (k+1))
            ax.errorbar(phi_bin_centers[df['phibin'][mask]], df['ratio'][mask], yerr=df['dratio'][mask], marker=markers[i], linestyle='None', label=df_key, color=colors[i])

        ax.axhline(1.0, color='gray', linestyle='--')

        ax.set_xlabel('$\phi$')
        ax.set_ylabel('Ratio')
        ax.set_ylim(0.0, 2.0)
        ax.set_xlim(0, 360)

        ax.legend()
        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')

    # Create a figure and axis objects
    fig, axes = plt.subplots(NumtBins, 1, figsize=(8, 6 * NumtBins), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t = {:.2f}".format(t_bin_centers[k]))

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])
            ax.scatter(df['phi'][mask], df['Q2'][mask], marker=markers[i], linestyle='None', label=df_key, color=colors[i])

        ax.set_xlabel('$\phi$')
        ax.set_ylabel('$Q^2$')
        ax.set_xlim(0, 360)
        ax.legend()
        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    
    # Create a figure and axis objects
    fig, axes = plt.subplots(NumtBins, 1, figsize=(8, 6 * NumtBins), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t = {:.2f}".format(t_bin_centers[k]))

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])
            ax.scatter(df['phi'][mask], df['W'][mask], marker=markers[i], linestyle='None', label=df_key, color=colors[i])

        ax.set_xlabel('$\phi$')
        ax.set_ylabel('W')
        ax.set_xlim(0, 360)
        ax.legend()
        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')

    # Create a figure and axis objects
    fig, axes = plt.subplots(NumtBins, 1, figsize=(8, 6 * NumtBins), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t = {:.2f}".format(t_bin_centers[k]))

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])
            ax.errorbar(df['th_cm'][mask], df['x_real'][mask], yerr=df['dx_real'][mask], marker=markers[i], linestyle='None', label=df_key, color=colors[i])

        ax.set_xlabel('$\theta_{cm}$')
        ax.set_ylabel('x_real')
        ax.legend()
        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')

    fig, axes = plt.subplots(NumtBins, 1, figsize=(8, 6 * NumtBins), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t = {:.2f}".format(t_bin_centers[k]))

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])
            ax.scatter(df['th_cm'][mask], df['x_mod'][mask], marker=markers[i], linestyle='None', label=df_key, facecolors='none' , edgecolors=colors[i])

        ax.set_xlabel('$\theta_{cm}$')
        ax.set_ylabel('x_mod')
        ax.legend()
        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    
    fig, axes = plt.subplots(NumtBins, 1, figsize=(8, 6 * NumtBins), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t = {:.2f}".format(t_bin_centers[k]))

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])
            ax.errorbar(df['phi'][mask], df['x_real'][mask], yerr=df['dx_real'][mask], marker=markers[i], linestyle='None', label=df_key, color=colors[i])

        ax.set_xlabel('$\phi$')
        ax.set_ylabel('x_real')
        ax.set_xlim(0, 360)
        ax.legend()
        # Add grid to subplot
        ax.grid(True, linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    
    fig, axes = plt.subplots(NumtBins, 1, figsize=(8, 6 * NumtBins), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t = {:.2f}".format(t_bin_centers[k]))

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])
            ax.scatter(df['phi'][mask], df['x_mod'][mask], marker=markers[i], linestyle='None', label=df_key, facecolors='none' , edgecolors=colors[i])

        ax.set_xlabel('$\phi$')
        ax.set_ylabel('x_mod')
        ax.set_xlim(0, 360)
        ax.legend()
        # Add grid to subplot
        ax.grid(True, linestyle='--', linewidth=0.5)
        
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    
    # Loop through t bins and plot data
    for k in range(NumtBins):

        # Create a figure and axis objects
        fig, axes = plt.subplots(1, 1, figsize=(8, 6), sharex=True)

        ax = axes
        ax.set_title("t = {:.2f}".format(t_bin_centers[k]))

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])
            ax.errorbar(df['phi'][mask], df['x_real'][mask], yerr=df['dx_real'][mask], marker=markers[i], linestyle='None', label=df_key, color=colors[i])
            ax.scatter(df['phi'][mask], df['x_mod'][mask], marker=markers[i], linestyle='None', label=df_key+" Model", facecolors='none' , edgecolors=colors[i])

        ax.set_xlabel('$\phi$')
        ax.set_ylabel('x_real')
        ax.set_xlim(0, 360)
        ax.legend()
        # Add grid to subplot
        ax.grid(True, linestyle='--', linewidth=0.5)
        # Set y-axis to logarithmic scale
        ax.set_yscale('log')

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
    
    fig, axes = plt.subplots(NumtBins, 1, figsize=(8, 6 * NumtBins), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t = {:.2f}".format(t_bin_centers[k]))

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])
            ax.errorbar(df['phi'][mask], df['x_real'][mask], yerr=df['dx_real'][mask], marker=markers[i], linestyle='None', label=df_key, color=colors[i])
            ax.scatter(df['phi'][mask], df['x_mod'][mask], marker=markers[i], linestyle='None', label=df_key+" Model", facecolors='none' , edgecolors=colors[i])

        ax.set_xlabel('$\phi$')
        ax.set_ylabel('x_real')
        ax.set_xlim(0, 360)
        ax.legend()
        # Add grid to subplot
        ax.grid(True, linestyle='--', linewidth=0.5)
        
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    
    fig, axes = plt.subplots(2, 2, figsize=(8, 6), sharex=True)

    for k, sig in enumerate(['sigL','sigT','sigLT','sigTT']):
        
        # Use integer division to get the correct subplot position
        ax = axes[k // 2, k % 2]
        formatted_sig = sig.replace("sig", "\sigma_{") + "}"
        ax.set_title("${}$".format(formatted_sig))
        for i, df_key in enumerate(['sep_file']):
            df = file_df_dict[df_key]
            print("="*50)
            model = []
            # Generate model for comparison
            for j, row in df.iterrows():
                print("-"*50)
                print("Data {} = {:.4e}".format(sig, row[sig]))
                inp_param = '{} {} {} {} {} '.format(Q2.replace("p","."), row['th_cm'], row['t'], row['Q2'], row['W'])+' '.join(param_arr)
                model.append(import_model(sig, inp_param))
            # Check that model sig is not all zeros
            if not all(element == 0 for element in model):
                ax.plot(df['t'], model, linestyle='-.', color='red', label='Model Fit')
            ax.errorbar(df['t'], df['{}'.format(sig)], yerr=df['d{}'.format(sig)], marker=markers[i], linestyle='None', label='Data', color=colors[i])
        ax.set_xlabel('t')
        ax.set_ylabel("${}$".format(formatted_sig))
        ax.set_xlim(tmin, tmax)
        ax.legend()
        # Add grid to subplot
        ax.grid(True, linestyle='--', linewidth=0.5)
    print("="*50)
        
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
