#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-05-23 00:25:30 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#

import pandas as pd
import root_pandas as rpd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit
import re
import sys, os

###############################################################################################################################################
'''
ltsep package import and pathing definitions
'''

OutFilename = "Q2_W_dep"

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
CACHEPATH=lt.CACHEPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH

foutname = OUTPATH+"/" + OutFilename + ".root"
fouttxt  = OUTPATH+"/" + OutFilename + ".txt"
outputpdf  = OUTPATH+"/" + OutFilename + ".pdf"

################################################################################################################################################
'''
Import separated xsects model
'''

pol_str = "pl"
ParticleType = "kaon"

sys.path.append("../models")
if pol_str == "pl" and ParticleType == "kaon":
    from sep_xsect_kaon_pl import import_model

################################################################################################################################################

Q2_lst = ["2p1", "3p0"]
W_lst = ["2p95", "3p14"]
LOEPS_lst = [0.2477, 0.3935]
HIEPS_lst = [0.7864, 0.6668]

comb_dict = {}

for Q2, W, LOEPS, HIEPS in zip(Q2_lst,W_lst, LOEPS_lst, HIEPS_lst):

    if Q2 == "2p1" and W == "2p95":
        inp_dir = CACHEPATH+"/{}/{}/Q{}W{}/".format(USER,ParticleType,Q2,W)+"2024May22_H11M45S25"
    elif Q2 == "3p0" and W == "3p14":
        inp_dir = CACHEPATH+"/{}/{}/Q{}W{}/".format(USER,ParticleType,Q2,W)+"trial_23/2024May06_H23M03S07"
    else:
        print("Error: No valid path for Q{}W{}...".format(Q2,W))
        sys.exit(2)
        
    ###############################################################################################################################################
    '''
    Import model parameterization
    '''

    param_file = '{}/parameters/par.{}_Q{}W{}.dat'.format(inp_dir, pol_str, Q2.replace("p",""), W.replace("p",""))

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
    fix_spacing(inp_dir+"/averages/avek.Q{}W{}.dat".format(Q2.replace("p",""), W.replace("p","")))
    fix_spacing(inp_dir+"/averages/aver.{}_Q{}W{}_{:.0f}.dat".format(pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100))
    fix_spacing(inp_dir+"/averages/aver.{}_Q{}W{}_{:.0f}.dat".format(pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100))
    fix_spacing(inp_dir+"/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat".format(pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100))
    fix_spacing(inp_dir+"/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat".format(pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100))
    fix_spacing(inp_dir+"/xsects/x_sep.{}_Q{}W{}.dat".format(pol_str, Q2.replace("p",""), W.replace("p","")))
    ################################################################################################################################################
    # Read in files and convert to dataframes

    file_df_dict = {}

    setting_file = inp_dir+"/list.settings".format(ParticleType)
    file_df_dict['setting_df'] = file_to_df(setting_file, ['POL', 'Q2', 'W', 'EPSVAL', 'thpq', 'TMIN', 'TMAX', 'NumtBins'],\
                                            filter_conditions={'Q2': float(Q2.replace("p",".")), 'W': float(W.replace("p","."))})

    tmin = file_df_dict['setting_df'].iloc[0]['TMIN']
    tmax = file_df_dict['setting_df'].iloc[0]['TMAX']

    try:
        with open("{}/t_bin_interval_Q{}W{}".format(inp_dir, Q2.replace("p",""), W.replace("p","")), "r") as file:
            # Read all lines from the file into a list
            all_lines = file.readlines()
            # Check if the file has at least two lines
            if len(all_lines) >= 2:
                # Extract the second line and remove leading/trailing whitespace
                t_bins = all_lines[1].split("\t")
                del t_bins[0]
                t_bins = np.array([float(element) for element in t_bins])
    except FileNotFoundError:
        print("{} not found...".format("{}/t_bin_interval_Q{}W{}".format(inp_dir, Q2.replace("p",""), W.replace("p",""))))
    except IOError:
        print("Error reading {}...".format("{}/t_bin_interval_Q{}W{}".format(inp_dir, Q2.replace("p",""), W.replace("p",""))))    

    t_bin_centers = (t_bins[:-1] + t_bins[1:]) / 2        
    file_df_dict['t_bin_centers'] = pd.DataFrame(t_bin_centers, columns=['t_bin_centers'])
    
    try:
        with open("{}/phi_bin_interval_Q{}W{}".format(inp_dir, Q2.replace("p",""), W.replace("p","")), "r") as file:
            # Read all lines from the file into a list
            all_lines = file.readlines()
            # Check if the file has at least two lines
            if len(all_lines) >= 2:
                # Extract the second line and remove leading/trailing whitespace
                phi_bins = all_lines[1].split("\t")
                del phi_bins[0]
                phi_bins = np.array([float(element) for element in phi_bins])
    except FileNotFoundError:
        print("{} not found...".format("{}/phi_bin_interval_Q{}W{}".format(inp_dir, Q2.replace("p",""), W.replace("p",""))))
    except IOError:
        print("Error reading {}...".format("{}/phi_bin_interval_Q{}W{}".format(inp_dir, Q2.replace("p",""), W.replace("p",""))))    

    phi_bin_centers = phi_bins
    file_df_dict['phi_bin_centers'] = pd.DataFrame(phi_bin_centers, columns=['phi_bin_centers'])

    for i,row in file_df_dict['setting_df'].iterrows():
        if row['Q2'] == float(Q2.replace("p",".")):
            file_df_dict['beam_file'] = file_to_df(inp_dir+"/beam/Eb_KLT.dat".format(ParticleType), ['ebeam', 'Q2', 'W', 'EPSVAL'])
            file_df_dict['avek_file'] = file_to_df(inp_dir+"/averages/avek.Q{}W{}.dat".format(Q2.replace("p",""), W.replace("p","")) \
                                                   , ['W', 'dW', 'Q2', 'dQ2', 't', 'dt', 'th_pos', "tbin"])        

            if row['EPSVAL'] == float(LOEPS):
                file_df_dict['aver_loeps'] = file_to_df( \
                                                         inp_dir+"/averages/aver.{}_Q{}W{}_{:.0f}.dat" \
                                                         .format(pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100) \
                                                         , ['ratio', 'dratio', 'phibin', 'tbin']).sort_values(by='tbin')
                if row['thpq'] < 0.0:
                    file_df_dict['kindata_loeps_{}'.format('right')] = file_to_df( \
                                                                                   inp_dir+"/kindata/kindata.{}_Q{}W{}_{:.0f}_{}.dat" \
                                                                                   .format(pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                           float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                                   , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
                if row['thpq'] > 0.0:
                    file_df_dict['kindata_loeps_{}'.format('left')] = file_to_df( \
                                                                                  inp_dir+"/kindata/kindata.{}_Q{}W{}_{:.0f}_+{}.dat" \
                                                                                  .format(pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                          float(LOEPS)*100, int(row['thpq']*1000)) \
                                                                                  , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
                if row['thpq'] == 0.0:
                    file_df_dict['kindata_loeps_{}'.format('center')] = file_to_df( \
                                                                                    inp_dir+"/kindata/kindata.{}_Q{}W{}_{:.0f}_+0000.dat" \
                                                                                    .format(pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                            float(LOEPS)*100) \
                                                                                    , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
                file_df_dict['unsep_file_loeps'] = file_to_df( \
                                                                inp_dir+"/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat" \
                                                                .format(pol_str, Q2.replace("p",""), W.replace("p",""), float(LOEPS)*100) \
                                                                , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 't', 'W', 'Q2']).sort_values(by='t')

            if row['EPSVAL'] == float(HIEPS):
                file_df_dict['aver_hieps'] = file_to_df( \
                                                         inp_dir+"/averages/aver.{}_Q{}W{}_{:.0f}.dat" \
                                                         .format(pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100) \
                                                         , ['ratio', 'dratio', 'phibin', 'tbin']).sort_values(by='tbin')
                if row['thpq'] < 0.0:
                    file_df_dict['kindata_hieps_{}'.format('right')] = file_to_df( \
                                                                                   inp_dir+"/kindata/kindata.{}_Q{}W{}_{:.0f}_{}.dat" \
                                                                                   .format(pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                           float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                                   , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
                if row['thpq'] > 0.0:
                    file_df_dict['kindata_hieps_{}'.format('left')] = file_to_df( \
                                                                                  inp_dir+"/kindata/kindata.{}_Q{}W{}_{:.0f}_+{}.dat" \
                                                                                  .format(pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                          float(HIEPS)*100, int(row['thpq']*1000)) \
                                                                                  , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
                if row['thpq'] == 0.0:
                    file_df_dict['kindata_hieps_{}'.format('center')] = file_to_df( \
                                                                                    inp_dir+"/kindata/kindata.{}_Q{}W{}_{:.0f}_+0000.dat" \
                                                                                    .format(pol_str, Q2.replace("p",""), W.replace("p",""), \
                                                                                            float(HIEPS)*100) \
                                                                                    , ['Q2', 'dQ2', 'W', 'dW', 't', 'dt'])
                file_df_dict['unsep_file_hieps'] = file_to_df( \
                                                                inp_dir+"/xsects/x_unsep.{}_Q{}W{}_{:.0f}.dat" \
                                                                .format(pol_str, Q2.replace("p",""), W.replace("p",""), float(HIEPS)*100) \
                                                                , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 't', 'W', 'Q2']).sort_values(by='t')
            file_df_dict['sep_file'] = file_to_df( \
                                                   inp_dir+"/xsects/x_sep.{}_Q{}W{}.dat" \
                                                   .format(pol_str, Q2.replace("p",""), W.replace("p","")) \
                                                   , ['sigL', 'dsigL', 'sigT', 'dsigT', 'sigLT', 'dsigLT', 'sigTT', 'dsigTT', 'chisq', 't', 'W', 'Q2', 'th_cm'])

    comb_dict["Q{}W{}".format(Q2,W)] = file_df_dict
    
print("\n\ncomb_dict")
print(comb_dict)
print("\n\n")

# Initialize the merged dictionary
merged_dict = {}

# Iterate over all subdictionaries
for subdict in comb_dict.values():
    for key, value in subdict.items():
        if key not in merged_dict:
            merged_dict[key] = []
        merged_dict[key].append(value)
        print("!!!!!!!", key, merged_dict[key])

# Flatten the merged dictionary
for key in merged_dict.keys():
    merged_dict[key] = pd.concat(merged_dict[key], ignore_index=True)
    print("$$$$$$$", key, merged_dict[key])
    
print("\n\nmerged_dict")
#print(merged_dict)
print("\n\n")

# Create a PdfPages object to manage the PDF file
with PdfPages(outputpdf) as pdf:

    ###

    # Define exponential function
    def exp_func(t, a, b):
        return a * np.exp(b * t)

    # Create a figure and axis objects for Q2 plot
    fig, axes = plt.subplots(1, 1, figsize=(12, 8), sharex=True)

    # Define markers and colors
    markers = ['x', 'o', '*', 'D'] # 'x'->x, 'o'->circle, '*'->star, 'D'->diamond
    colors = ['black', 'red']
    
    ax = axes
    #ax.set_title("$Q^2$={:.1f}, W={:.2f}".format(float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

    for i, df_key in enumerate(['kindata_loeps_{}'.format('center'), 'kindata_hieps_{}'.format('center')]):
        df = merged_dict[df_key]
        if "hi" in df_key:
            df_key = "High $\epsilon$"
        else:
            df_key = "Low $\epsilon$"

        print("\n\n",df_key,"\nt_bin_centers\n",merged_dict['t_bin_centers']['t_bin_centers'], "\nQ2\n", df['Q2'])
        ax.errorbar(merged_dict['t_bin_centers']['t_bin_centers'], df['Q2'], yerr=df['dQ2'], marker=markers[i], linestyle='None', label=df_key, color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)

        # Fit the data using exponential function
        popt, _ = curve_fit(exp_func, merged_dict['t_bin_centers']['t_bin_centers'], df['Q2'])
        fit_line = exp_func(merged_dict['t_bin_centers']['t_bin_centers'], *popt)
        ax.plot(merged_dict['t_bin_centers']['t_bin_centers'], fit_line, linestyle='-', color=colors[i], label="{0} Fit: Q(t) = {1:.2f}e^({2:.2f}t)".format(df_key, popt[0], popt[1]))

    ax.set_xlabel('-t', fontsize=24)
    ax.set_ylabel('$Q^2$', fontsize=24)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)        
    ax.set_xlim(tmin, tmax)
    ax.legend(fontsize=16)
    # Add grid
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')

    # Create a figure and axis objects for W plot
    fig, axes = plt.subplots(1, 1, figsize=(12, 8), sharex=True)

    ax = axes
    #ax.set_title("$Q^2$={:.1f}, W={:.2f}".format(float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

    for i, df_key in enumerate(['kindata_loeps_{}'.format('center'), 'kindata_hieps_{}'.format('center')]):
        df = merged_dict[df_key]
        if "hi" in df_key:
            df_key = "High $\epsilon$"
        else:
            df_key = "Low $\epsilon$"

        ax.errorbar(merged_dict['t_bin_centers']['t_bin_centers'], df['W'], yerr=df['dW'], marker=markers[i], linestyle='None', label=df_key, color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)

        # Fit the data using exponential function
        popt, _ = curve_fit(exp_func, merged_dict['t_bin_centers']['t_bin_centers'], df['W'])
        fit_line = exp_func(merged_dict['t_bin_centers']['t_bin_centers'], *popt)
        ax.plot(merged_dict['t_bin_centers']['t_bin_centers'], fit_line, linestyle='-', color=colors[i], label="{0} Fit: W(t) = {1:.2f}e^({2:.2f}t)".format(df_key, popt[0], popt[1]))

    ax.set_xlabel('-t', fontsize=24)
    ax.set_ylabel('W', fontsize=24)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)        
    ax.set_xlim(tmin, tmax)
    ax.legend(fontsize=16)
    # Add grid
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)

    for k, sig in enumerate(['sigL','sigT','sigLT','sigTT']):
        
        # Use integer division to get the correct subplot position
        ax = axes[k // 2, k % 2]
        formatted_sig = sig.replace("sig", "\sigma_{") + "}"
        ax.set_title("${}$".format(formatted_sig), fontsize=24)
        for i, df_key in enumerate(['sep_file']):
            df = file_df_dict[df_key]
            if "hi" in df_key:
                df_key = "High $\epsilon$"
            else:
                df_key = "Low $\epsilon$"
                
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
                ax.plot(merged_dict['t_bin_centers']['t_bin_centers'], model, linestyle='-.', color='red', label='Model Fit')
            ax.errorbar(merged_dict['t_bin_centers']['t_bin_centers'], df['{}'.format(sig)], yerr=df['d{}'.format(sig)], marker=markers[i], linestyle='None', label='Data', color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)
        ax.set_xlabel('t')
        ax.set_ylabel("${}$".format(formatted_sig))
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.set_xlim(tmin, tmax)
        ax.legend(fontsize=24)
        # Add grid to subplot
        ax.grid(True, linestyle='--', linewidth=0.5)
    print("="*50)
        
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    
    ###

