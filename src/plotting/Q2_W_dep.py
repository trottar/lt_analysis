#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-11 05:13:41 trottar"
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
import math
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

TEMP_CACHEPATH=f"{OUTPATH}/cache_transfer"

foutname = OUTPATH+"/" + OutFilename + ".root"
fouttxt  = OUTPATH+"/" + OutFilename + ".txt"
outputpdf  = OUTPATH+"/" + OutFilename + ".pdf"

##################################################################################################################################################
# Importing utility functions

sys.path.append("../utility")
from utility import show_pdf_with_evince

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

# Define constants
PI = math.pi
m_p = 0.93827231
m_n = 0.93956541
mkpl = 0.493677

settings = {
    'set_1': {
        'Q2': '2p1',
        'W': '2p95',
        'LOEPS': 0.2477,
        'HIEPS': 0.7864,
        #'inp_dir': "trial_9/2024July25_H17M05S03" # i=1
        #'inp_dir': "trial_9/2024July26_H05M28S47" # i=12
        #'inp_dir': "2024September11_H10M39S26" # i=2, 4p4 parameterization
        'inp_dir': "2024October12_H03M42S23" # i=10, w_fac->(n=2.45)
    },
    'set_2': {
        'Q2': '3p0',
        'W': '2p32',
        'LOEPS': 0.5736,
        'HIEPS': 0.8791,
        #'inp_dir': "2024September11_H12M56S06" # i=2, 4p4 parameterization
        #'inp_dir': "2024September12_H14M27S29" # i=5, 4p4 parameterization
        'inp_dir': "2024October12_H08M08S16" # i=10, w_fac->(n=3.40)
    },    
    'set_3': {
        'Q2': '3p0',
        'W': '3p14',
        'LOEPS': 0.3935,
        'HIEPS': 0.6668,
        #'inp_dir': "trial_30/2024July25_H17M19S51" # i=1
        #'inp_dir': "trial_30/2024July26_H08M04S15" # i=12
        #'inp_dir': "2024September11_H16M37S05" # i=2, 4p4 parameterization
        #'inp_dir': "2024September12_H15M13S17" # i=5, 4p4 parameterization
        'inp_dir': "trial_32/2024October15_H20M22S46" # i=10, w_fac->(n=2.25)
    },    
    'set_4': {
        'Q2': '4p4',
        'W': '2p74',
        'LOEPS': 0.4805,
        'HIEPS': 0.7148,
        #'inp_dir': "2024August13_H17M30S09" # i=1
        #'inp_dir': "trial_1/2024August13_H22M17S21" # i=12
        #'inp_dir': "2024September11_H17M50S25" # i=2, 4p4 parameterization
        #'inp_dir': "2024September12_H15M58S25" # i=5, 4p4 parameterization
        'inp_dir': "2024October12_H00M05S21" # i=10, w_fac->(n=2.65)
    },    
    'set_5': {
        'Q2': '5p5',
        'W': '3p02',
        'LOEPS': 0.1838,
        'HIEPS': 0.5291,
        #'inp_dir': "trial_14/2024July25_H17M34S49" # i=1
        #'inp_dir': "trial_14/2024July26_H10M45S46" # i=12
        #'inp_dir': "2024September11_H19M05S14" # i=2, 4p4 parameterization
        #'inp_dir': "2024September12_H16M39S54" # i=5, 4p4 parameterization
        'inp_dir': "" # i=10, w_fac->(n=2.45)
    }
}

comb_dict = {}

for key, values in settings.items():
    #if key in ('set_1', 'set_2', 'set_3', 'set_4', 'set_5'):
    if key in ('set_1', 'set_2', 'set_3', 'set_4'):
        Q2, W, LOEPS, HIEPS, inp_dir = values.values()
    else:
        continue

    inp_dir = TEMP_CACHEPATH+"/{}/{}/Q{}W{}/".format(USER,ParticleType,Q2,W)+inp_dir
    
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

# Flatten the merged dictionary
for key in merged_dict.keys():
    merged_dict[key] = pd.concat(merged_dict[key], ignore_index=True)
    print("-"*10, key, "-"*10, "\n", merged_dict[key])
    
print("\n\n")

# Redefine tmin and tmax
# Calculate the first tmin and tmax based on the data
tmin_initial = merged_dict['sep_file']['t'].min() - 0.1
tmax_initial = merged_dict['sep_file']['t'].max() + 0.1

# Define the list of (tmin, tmax) pairs, including the initial pair
tmin_tmax_pairs = [
    (tmin_initial, tmax_initial),  # Calculated from merged_dict['sep_file']['t']
    (0.15, 0.2),                   # Q2=2.115+3.0, 1
    (0.2, 0.25),                   # Q2=2.115+3.0, 2
    (0.3, 0.35),                   # Q2=2.115+3.0, 3
    (0.425, 0.5),                  # Q2=3.0+4.4+5.5
    (0.85, 0.9)                    # Q2=4.4+5.5
]

# Loop through the list and set tmin and tmax for each iteration
for tmin, tmax in tmin_tmax_pairs:
    # Use tmin and tmax in your calculations here
    print(f"tmin: {tmin:.3f}, tmax: {tmax:.3f}")

    tmin  = float(f"{tmin:.3f}")
    tmax  = float(f"{tmax:.3f}")    

    outputpdf = outputpdf.replace(".pdf", f"_tmin{str(tmin).replace('.','p')}-tmax{str(tmax).replace('.','p')}.pdf")
    
    # Create a PdfPages object to manage the PDF file
    with PdfPages(outputpdf) as pdf:

        ###

        # Define exponential function
        def exp_func(t, a, b):
            return a * np.exp(b * t)

        def sigl_func(data, p1, p2):
            q2, t = data
            ft = abs(t) / (abs(t) + mkpl**2)**2 # pole term
            Qdep_L=q2/(1.0+(1.77*q2)+0.12*(q2**2))
            sigl=(p1*Qdep_L*ft)*np.exp(-p2*(abs(t)))
            return sigl

        def sigt_func(data, p5, p6, p7, p8):
            q2, t = data            
            Qdep_T=(np.exp(-q2**2))/q2
            sigt=(p5*np.exp(-p6*(abs(t)))+p7*(abs(t)))*(Qdep_T**p8)
            ##
            ##ft = abs(t) / (abs(t) + mkpl**2)**2 # pole term
            ##sigt=(p5*np.exp(-p6*(abs(t)))+p7)*(ft*(abs(t)))*(Qdep_T**p8) # Testing
            return sigt

        def siglt_func(data, p9, p10):
            q2, t, theta = data
            siglt=(p9/(1+q2))*np.sin(theta*(PI/180))*np.exp(-p10*(abs(t)))
            ##
            ##ft = abs(t) / (abs(t) + mkpl**2)**2 # pole term
            ##siglt=(p9/(1+q2))*np.sin(theta*(PI/180))*ft*np.exp(-p10*(abs(t))) # Testing
            return siglt

        def sigtt_func(data, p13, p14):
            q2, t, theta = data
            ft = abs(t) / (abs(t) + mkpl**2)**2 # pole term
            sigtt=(p13/(1+q2))*(np.sin(theta*(PI/180))**2)*ft*np.exp(-p14*(q2))
            return sigtt

        # Create a figure and axis objects for Q2 plot
        fig, axes = plt.subplots(1, 1, figsize=(12, 8), sharex=True)

        # Define markers and colors
        markers = ['x', 'o', '*', 'D', '^', '+'] # 'x'->x, 'o'->circle, '*'->star, 'D'->diamond
        colors = ['black', 'red', 'blue', 'green', 'purple']

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
            #ax.plot(merged_dict['t_bin_centers']['t_bin_centers'], fit_line, linestyle='-', color=colors[i], label="{0} Fit: Q(t) = {1:.2f}e^({2:.2f}t)".format(df_key, popt[0], popt[1]))

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
            #ax.plot(merged_dict['t_bin_centers']['t_bin_centers'], fit_line, linestyle='-', color=colors[i], label="{0} Fit: W(t) = {1:.2f}e^({2:.2f}t)".format(df_key, popt[0], popt[1]))

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
            df = merged_dict["sep_file"]

            W_ref = 3.0 # Scale data to W=3.0 GeV
            n = 2.0
            #n = 2.25 # From Phys. Rev. C 85, 018202 (2012), kaon
            #n = 2.41 # From Phys. Rev. C 85, 018202 (2012), pion
            #w_scale_factor = ((W_ref**2-mkpl**2)**n)/((df['W']**2-mkpl**2)**n)
            w_scale_factor = ((W_ref**2-mkpl**2)**n)

            scaled_sig = df['{}'.format(sig)]*w_scale_factor
            d_scaled_sig = df['d{}'.format(sig)]*w_scale_factor

            tolerance = 0.5

            if (abs(df['Q2'] - 2.115) < tolerance).any():
                mask = abs(df['Q2'] - 2.115) < tolerance
                ax.errorbar(df.loc[mask, 't'], scaled_sig[mask], yerr=d_scaled_sig[mask], 
                            marker=markers[0], linestyle='', label='$Q^2$=2.115, W=2.95', 
                            color=colors[0], markeredgecolor=colors[0], 
                            markerfacecolor='none', capsize=2)

            if (abs(df['Q2'] - 3.0) < tolerance).any():
                if (abs(df['W'] - 2.32) < tolerance).any():
                    mask = (abs(df['Q2'] - 3.0) < tolerance) & (abs(df['W'] - 2.32) < tolerance)
                    ax.errorbar(df.loc[mask, 't'], scaled_sig[mask], yerr=d_scaled_sig[mask], 
                                marker=markers[1], linestyle='', label='$Q^2$=3.0, W=2.32', 
                                color=colors[1], markeredgecolor=colors[1], 
                                markerfacecolor='none', capsize=2)
                if (abs(df['W'] - 3.14) < tolerance).any():
                    mask = (abs(df['Q2'] - 3.0) < tolerance) & (abs(df['W'] - 3.14) < tolerance)
                    ax.errorbar(df.loc[mask, 't'], scaled_sig[mask], yerr=d_scaled_sig[mask], 
                                marker=markers[2], linestyle='', label='$Q^2$=3.0, W=3.14', 
                                color=colors[2], markeredgecolor=colors[2], 
                                markerfacecolor='none', capsize=2)                    

            if (abs(df['Q2'] - 4.4) < tolerance).any():
                mask = abs(df['Q2'] - 4.4) < tolerance
                ax.errorbar(df.loc[mask, 't'], scaled_sig[mask], yerr=d_scaled_sig[mask], 
                            marker=markers[3], linestyle='', label='$Q^2$=4.4, W=2.74', 
                            color=colors[3], markeredgecolor=colors[3], 
                            markerfacecolor='none', capsize=2)

            if (abs(df['Q2'] - 5.5) < tolerance).any():
                mask = abs(df['Q2'] - 5.5) < tolerance
                ax.errorbar(df.loc[mask, 't'], scaled_sig[mask], yerr=d_scaled_sig[mask], 
                            marker=markers[4], linestyle='', label='$Q^2$=5.5, W=3.02', 
                            color=colors[4], markeredgecolor=colors[4], 
                            markerfacecolor='none', capsize=2)

            # Generate points for smooth curve
            #q2_fit = np.linspace(df['Q2'].min(), df['Q2'].max(), 100)
            #t_fit = np.linspace(tmin, tmax, 100)
            #theta_fit = np.linspace(df['th_cm'].min(), df['th_cm'].max(), 100)
            q2_fit = np.linspace(df['Q2'].min(), df['Q2'].max(), 10000)
            t_fit = np.linspace(0.0001, 2.0, 10000)
            theta_fit = np.linspace(0.0, 360.0, 10000)
                
            if sig == "sigL":
                # Perform exponential fit
                popt, _ = curve_fit(sigl_func, (df['Q2'], df['t']), scaled_sig, sigma=d_scaled_sig, absolute_sigma=True, maxfev = 100000)
                p1, p2 = popt
                param_str = f"{p1:.3e}, {p2:.3e}"
                y_fit = sigl_func((q2_fit, t_fit), p1, p2)
            if sig == "sigT":
                # Perform exponential fit
                popt, _ = curve_fit(sigt_func, (df['Q2'], df['t']), scaled_sig, sigma=d_scaled_sig, absolute_sigma=True, maxfev = 100000)
                p5, p6, p7, p8 = popt
                param_str = f"{p5:.3e}, {p6:.3e}, {p7:.3e}, {p8:.3e}"
                y_fit = sigt_func((q2_fit, t_fit), p5, p6, p7, p8)
            if sig == "sigLT":
                # Perform exponential fit
                popt, _ = curve_fit(siglt_func, (df['Q2'], df['t'], df['th_cm']), scaled_sig, sigma=d_scaled_sig, absolute_sigma=True, maxfev = 100000)
                p9, p10 = popt
                param_str = f"{p9:.3e}, {p10:.3e}"
                y_fit = siglt_func((q2_fit, t_fit, theta_fit), p9, p10)
            if sig == "sigTT":
                # Perform exponential fit
                popt, _ = curve_fit(sigtt_func, (df['Q2'], df['t'], df['th_cm']), scaled_sig, sigma=d_scaled_sig, absolute_sigma=True, maxfev = 100000)
                p13, p14 = popt
                param_str = f"{p13:.3e}, {p14:.3e}"
                y_fit = sigtt_func((q2_fit, t_fit, theta_fit), p13, p14)

            # Plot the fit
            ax.plot(t_fit, y_fit, '-', color='xkcd:light blue', label=f'Fit: {param_str}')

            ax.set_xlabel('t')
            ax.set_ylabel("${}$".format(formatted_sig))
            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=16)        
            ax.set_xlim(tmin-0.01, tmax+0.01)
            ax.legend(fontsize=8)
            # Add grid to subplot
            ax.grid(True, linestyle='--', linewidth=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')

        fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

        for k, eps_val in enumerate(['hi','lo']):

            # Use integer division to get the correct subplot position
            ax = axes[k]
            eps_str = "High" if eps_val == "hi" else "Low"
            ax.set_title(f"{eps_str} Ratio", fontsize=24)
            df = merged_dict[f"unsep_file_{eps_val}eps"]

            ratios = df['x_real']/df['x_mod']
            errors = df['dx_real']/df['x_mod']
            non_zero_mask = (ratios != 0) & (errors != 0)
            ratios = ratios[non_zero_mask]
            errors = errors[non_zero_mask]

            # Use x_increment for x-axis values
            x_values = np.arange(len(ratios[non_zero_mask]))
            
            tolerance = 0.5
            
            if (abs(df[non_zero_mask]['Q2'] - 2.115) < tolerance).any():
                mask = abs(df[non_zero_mask]['Q2'] - 2.115) < tolerance
                ax.errorbar(x_values[mask], ratios[mask], yerr=errors[mask], 
                            marker=markers[0], linestyle='', label='$Q^2$=2.115, W=2.95', 
                            color=colors[0], markeredgecolor=colors[0], 
                            markerfacecolor='none', capsize=2)

            if (abs(df[non_zero_mask]['Q2'] - 3.0) < tolerance).any():
                if (abs(df[non_zero_mask]['W'] - 2.32) < tolerance).any():
                    mask = (abs(df[non_zero_mask]['Q2'] - 3.0) < tolerance) & (abs(df[non_zero_mask]['W'] - 2.32) < tolerance)
                    ax.errorbar(x_values[mask], ratios[mask], yerr=errors[mask], 
                                marker=markers[1], linestyle='', label='$Q^2$=3.0, W=2.32', 
                                color=colors[1], markeredgecolor=colors[1], 
                                markerfacecolor='none', capsize=2)
                if (abs(df[non_zero_mask]['W'] - 3.14) < tolerance).any():
                    mask = (abs(df[non_zero_mask]['Q2'] - 3.0) < tolerance) & (abs(df[non_zero_mask]['W'] - 3.14) < tolerance)
                    ax.errorbar(x_values[mask], ratios[mask], yerr=errors[mask], 
                                marker=markers[2], linestyle='', label='$Q^2$=3.0, W=3.14', 
                                color=colors[2], markeredgecolor=colors[2], 
                                markerfacecolor='none', capsize=2)                    

            if (abs(df[non_zero_mask]['Q2'] - 4.4) < tolerance).any():
                mask = abs(df[non_zero_mask]['Q2'] - 4.4) < tolerance
                ax.errorbar(x_values[mask], ratios[mask], yerr=errors[mask], 
                            marker=markers[3], linestyle='', label='$Q^2$=4.4, W=2.74', 
                            color=colors[3], markeredgecolor=colors[3], 
                            markerfacecolor='none', capsize=2)

            if (abs(df[non_zero_mask]['Q2'] - 5.5) < tolerance).any():
                mask = abs(df[non_zero_mask]['Q2'] - 5.5) < tolerance
                ax.errorbar(x_values[mask], ratios[mask], yerr=errors[mask], 
                            marker=markers[4], linestyle='', label='$Q^2$=5.5, W=3.02', 
                            color=colors[4], markeredgecolor=colors[4], 
                            markerfacecolor='none', capsize=2)

            ax.set_xlabel('t-$\phi$ bin')
            ax.set_ylabel("Ratio")
            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=16)
            ax.legend(fontsize=8)
            # Add grid to subplot
            ax.grid(True, linestyle='--', linewidth=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
        
        fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
        for k, sig in enumerate(['sigL','sigT','sigLT','sigTT']):

            # Use integer division to get the correct subplot position
            ax = axes[k // 2, k % 2]
            formatted_sig = sig.replace("sig", "\sigma_{") + "}"
            ax.set_title("${}$".format(formatted_sig), fontsize=24)

            df = merged_dict["sep_file"]
            df = df[(df['t'] >= tmin) & (df['t'] <= tmax)]
            cut_str = f"t = [{tmin:.3f}, {tmax:.3f}]"

            W_ref = 3.0 # Scale data to W=3.0 GeV
            n = 2.0
            #n = 2.25 # From Phys. Rev. C 85, 018202 (2012), kaon
            #n = 2.41 # From Phys. Rev. C 85, 018202 (2012), pion
            #w_scale_factor = ((W_ref**2-mkpl**2)**n)/((df['W']**2-mkpl**2)**n)
            w_scale_factor = ((W_ref**2-mkpl**2)**n)

            scaled_sig = df['{}'.format(sig)]*w_scale_factor
            d_scaled_sig = df['d{}'.format(sig)]*w_scale_factor

            print("\n\n",df[['t', 'Q2', '{}'.format(sig), 'd{}'.format(sig)]])
            ax.errorbar(df['Q2'], scaled_sig, yerr=d_scaled_sig, marker=markers[0], linestyle='None', label=cut_str, color=colors[0], markeredgecolor=colors[0], markerfacecolor='none', capsize=2)

            # Perform exponential fit
            popt, _ = curve_fit(exp_func, df['Q2'], scaled_sig, sigma=d_scaled_sig, absolute_sigma=True, maxfev = 100000)

            # Generate points for smooth curve
            x_fit = np.linspace(df['Q2'].min(), df['Q2'].max(), 100)
            y_fit = exp_func(x_fit, *popt)

            # Plot the fit
            ax.plot(x_fit, y_fit, '-', color='xkcd:light blue', label=f'Fit: {popt[0]:.2e}*exp({popt[1]:.2f}*Q^2)')

            ax.set_xlabel('$Q^2$')
            ax.set_ylabel("${}$".format(formatted_sig))
            ax.tick_params(axis='x', labelsize=16)
            ax.tick_params(axis='y', labelsize=16)        
            ax.set_xlim(df['Q2'].min()-0.1, df['Q2'].max()+0.1)
            ax.legend(fontsize=12)
            # Add grid to subplot
            ax.grid(True, linestyle='--', linewidth=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')    
        ###

    show_pdf_with_evince(outputpdf)
