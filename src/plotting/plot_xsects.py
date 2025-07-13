#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2025-05-20 15:23:59 trottar"
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
import warnings
from collections import defaultdict
from scipy.optimize import curve_fit
import math
import re
import sys, os

##################################################################################################################################################
# Check the number of arguments provided to the script

if len(sys.argv)-1!=10:
    print("!!!!! ERROR !!!!!\n Expected 10 arguments\n Usage is with - ParticleType POL Q2 W LOEPS HIEPS NumtBins NumPhiBins KIN OutUnsepxsectsFilename\n!!!!! ERROR !!!!!")
    sys.exit(1)

###############################################################################################################################################    
# Suppress the OptimizeWarning
warnings.filterwarnings("ignore", category=RuntimeWarning)    
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

lt=Root(os.path.realpath(__file__),"Plot_LTSep")

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
#t_bin_centers = t_bins
#print("!!!!!!!!!!NumtBins",NumtBins)
#print("!!!!!!!!!!1",len(t_bins),t_bins)
#print("!!!!!!!!!!2",len(t_bin_centers),t_bin_centers)

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
                                                            , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 't', 'W', 'Q2']).sort_values(by='t')

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
                                                            , ['x_real', 'dx_real', 'x_mod', 'eps', 'th_cm', 'phi', 't', 'W', 'Q2']).sort_values(by='t')
        file_df_dict['sep_file'] = file_to_df( \
                                               LTANAPATH+"/src/{}/xsects/x_sep.{}_Q{}W{}.dat" \
                                               .format(ParticleType, pol_str, Q2.replace("p",""), W.replace("p","")) \
                                               , ['sigL', 'dsigL', 'sigT', 'dsigT', 'sigLT', 'dsigLT', 'sigTT', 'dsigTT', 'chisq', 't', 'W', 'Q2', 'th_cm'])
            
################################################################################################################################################

# Create a PdfPages object to manage the PDF file
with PdfPages(outputpdf) as pdf:
    
    # Create a figure and axis objects
    fig, axes = plt.subplots(NumtBins, 1, figsize=(12, 8), sharex=True)

    # Define markers and colors
    markers = ['x', 'o', '*', 'D'] # 'x'->x, 'o'->circle, '*'->star, 'D'->diamond
    colors = ['black', 'red']

    # Loop through t bins and plot data
    for k in range(NumtBins):
        
        ax = axes[k]
        ax.set_title("t={:.3f}, $Q^2$={:.1f}, W={:.2f}".format(t_bin_centers[k], float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

        for i, df_key in enumerate(['aver_loeps', 'aver_hieps']):
            df = file_df_dict[df_key]
            if "hi" in df_key:
                df_key = "High $\epsilon$"
            else:
                df_key = "Low $\epsilon$"
                
            mask = (df['tbin'] == (k+1))
            ax.errorbar(phi_bin_centers[df['phibin'][mask]], df['ratio'][mask], yerr=np.maximum(df['dratio'][mask], 0.0), marker=markers[i], linestyle='None', label=df_key, color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)

        ax.axhline(1.0, color='gray', linestyle='--')

        ax.set_xlabel('$\phi$', fontsize=24)
        ax.set_ylabel('Ratio', fontsize=24)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.set_ylim(0.0, 2.0)
        ax.set_xlim(-185, 185)

        ax.legend(fontsize=14)
        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')

    # Loop through t bins and plot data
    for k in range(NumtBins):
        # Create a figure and axis objects
        fig, axes = plt.subplots(1, 1, figsize=(12, 8), sharex=True)
        ax = axes
        ax.set_title("t={:.3f}, $Q^2$={:.1f}, W={:.2f}".format(t_bin_centers[k], float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

        for i, df_key in enumerate(['aver_loeps', 'aver_hieps']):
            df = file_df_dict[df_key]
            if "hi" in df_key:
                epsilon_label = "High $\epsilon$"
            else:
                epsilon_label = "Low $\epsilon$"

            mask = (df['tbin'] == (k+1))

            # Filter out zero values
            ratios = df['ratio'][mask]
            errors = df['dratio'][mask]
            non_zero_mask = (ratios != 0) & (errors != 0)
            ratios = ratios[non_zero_mask]
            errors = errors[non_zero_mask]

            if len(ratios) > 0:  # Check if we have any non-zero data points
                weights = 1 / (errors ** 2)

                weighted_average = np.average(ratios, weights=weights)
                weighted_error = np.sqrt(1 / np.sum(weights))

                # Update the legend label to include the weighted average and its error
                legend_label = "{} (Avg: {:.3e} ± {:.3e})".format(epsilon_label, weighted_average, weighted_error)
            else:
                legend_label = "{} (No valid data)".format(epsilon_label)

            ax.errorbar(phi_bin_centers[df['phibin'][mask][non_zero_mask]], ratios, yerr=np.maximum(errors, 0.0),
                        marker=markers[i], linestyle='None', label=legend_label, 
                        color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)

        ax.axhline(1.0, color='gray', linestyle='--')
        ax.set_xlabel('$\phi$', fontsize=24)
        ax.set_ylabel('Ratio', fontsize=24)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.set_xlim(-185, 185)
        ax.set_ylim(0.0, 2.0)
        ax.legend(fontsize=14)

        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')

    ##########
    # Fix ratio data vs Phi

    # functional form in φ [degrees]
    def R_model(phi_deg, A, B, C):
        phi = np.deg2rad(phi_deg)
        return A + B * np.cos(phi) + C * np.cos(2*phi)
    
    # containers for the coefficients we’ll harvest inside your φ-loop
    coeffs = defaultdict(lambda: {'t': [], 'A': [], 'B': [], 'C': []})   

    # Loop through t bins and plot data + fits
    for k in range(NumtBins):
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title(
            "t={:.3f}, $Q^2$={:.1f}, W={:.2f}".format(
                t_bin_centers[k],
                float(Q2.replace("p",".")),
                float(W.replace("p","."))
            ),
            fontsize=24
        )

        for i, df_key in enumerate(['aver_loeps','aver_hieps']):
            df = file_df_dict[df_key]
            epsilon_label = "High $\epsilon$" if "hi" in df_key else "Low $\epsilon$"

            mask = (df['tbin'] == (k+1))
            φ  = phi_bin_centers[df['phibin'][mask]]
            R  = df['ratio'][mask]
            dR = df['dratio'][mask]
            good = (R != 0) & (dR != 0)
            φ, R, dR = φ[good], R[good], dR[good]

            # plot the data points
            ax.errorbar(
                φ, R, yerr=dR,
                marker=markers[i], linestyle='None',
                label=epsilon_label,
                color=colors[i], markeredgecolor=colors[i],
                markerfacecolor='none', capsize=2
            )

            if len(R) >= 3:   # need at least 3 points to fit A,B,C
                # do the weighted fit
                popt, pcov = curve_fit(
                    R_model, φ, R,
                    sigma=dR, absolute_sigma=True,
                    p0=[1.0, 0.1, 0.1]
                )
                A, B, C = popt
                coeffs[df_key]['t'].append(t_bin_centers[k])
                coeffs[df_key]['A'].append(A)
                coeffs[df_key]['B'].append(B)
                coeffs[df_key]['C'].append(C)                
                # smooth curve for plotting
                φ_smooth = np.linspace(-180, 180, 361)
                R_smooth = R_model(φ_smooth, A, B, C)

                # plot the fit curve
                fit_label = f"Fit {epsilon_label}: A={A:.3f}, B={B:.3f}, C={C:.3f}"
                ax.plot(
                    φ_smooth, R_smooth,
                    linestyle='-',
                    color=colors[i],
                    label=fit_label
                )

        # rest of your formatting
        ax.axhline(1.0, color='gray', linestyle='--')
        ax.set_xlabel('$\phi$ (deg)', fontsize=24)
        ax.set_ylabel('Ratio', fontsize=24)
        ax.set_xlim(-185, 185)
        ax.set_ylim(0.0, 2.0)
        ax.tick_params(axis='both', labelsize=16)
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        ax.legend(fontsize=14)
        plt.tight_layout(rect=[0,0,1,0.96])
        pdf.savefig(fig, bbox_inches='tight')

    ##########
    # Plot the fit parameters A, B, C vs t and fit them with three models

    # analytic t-dependence models
    def f_linear(t, a, b):      return a + b * t
    def f_exponential(t, a, b): return a * np.exp(-abs(b * t))
    def f_inverse(t, a, b):     return a + b / t
    #def f_custom1(t, a, b, c):      return a * (t/(t + 0.493677**2)**2 + 1/abs(t**c)) * np.exp(-abs(b * t))
    def f_custom1(t, a, b, c):      return a * (t/(t + 0.493677**2)**2 + t) * np.exp(-abs(b * t))
    #def f_custom2(t, a, b):      return (a/t) * np.exp(-abs(b * t))
    def f_custom2(t, a, b):      return (a * t) * np.exp(-abs(b * t**2))

    component_labels = ['A', 'B', 'C']
    fit_styles = {
        #'Linear':      (f_linear,      '-',  1.00),
        #'Exponential': (f_exponential, '--', 0.95),
        #'1/t':         (f_inverse,     ':',  0.95),
        'Custom1':      (f_custom1,      '-',  1.00),
        'Custom2':      (f_custom2,      '--',  0.95) 
    }

    # first high-ε, then low-ε so the PDF order is exactly as requested
    page_order = [('aver_hieps', 'High ε', colors[1]),
                ('aver_loeps', 'Low ε',  colors[0])]

    t_dense = np.linspace(min(t_bin_centers), max(t_bin_centers), 400)

    for df_key, eps_label, base_color in page_order:
        if not coeffs[df_key]['t']:         # skip if that ε-sample had no fits
            continue

        fig, axes = plt.subplots(
            3, 1, figsize=(11, 13), sharex=True,
            gridspec_kw=dict(hspace=0.30)
        )

        for ax, comp in zip(axes, component_labels):
            t_arr = np.array(coeffs[df_key]['t'])
            y_arr = np.array(coeffs[df_key][comp])

            # scatter points
            ax.scatter(t_arr, y_arr, marker='o', s=70,
                    facecolors='none', edgecolors=base_color,
                    label='data')

            if len(t_arr) >= 2:             # need ≥2 points to fit in t
                for name, (func, style, alpha) in fit_styles.items():
                    try:
                        p, _ = curve_fit(func, t_arr, y_arr, maxfev=6000)
                        ax.plot(t_dense, func(t_dense, *p),
                                linestyle=style, linewidth=2,
                                color=base_color, alpha=alpha,
                                label=f'{name}: {p[0]:.3g}, {p[1]:.3g}')
                    except RuntimeError:
                        pass    # silently ignore failed model

            ax.set_ylabel(comp, fontsize=15)
            ax.grid(True, linestyle='--', linewidth=0.6, alpha=0.7)
            ax.legend(fontsize=9)

        axes[0].set_title(f'{eps_label} coefficients vs $t$', fontsize=19, pad=10)
        axes[-1].set_xlabel('$t$', fontsize=15)

        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

    ###########

    ### HERE 0

    # Create a single figure and axis object for all phi bins
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_title(f"$Q^2$={float(Q2.replace('p', '.'))}, W={float(W.replace('p', '.'))}", fontsize=24)
    
    # Loop through t bins and plot data
    for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
        df = file_df_dict[df_key]
        if "hi" in df_key:
            epsilon_label = "High $\epsilon$" if k == 0 else ""
            epsilon_fit_color = "r-"
        else:
            epsilon_label = "Low $\epsilon$" if k == 0 else ""
            epsilon_fit_color = "-"

        ratios = df['x_real']/df['x_mod']
        errors = df['dx_real']/df['x_mod']
        non_zero_mask = (ratios != 0) & (errors != 0)
        ratios = ratios[non_zero_mask]
        errors = errors[non_zero_mask]        

        # Use x_increment for x-axis values
        x_values = np.arange(0, len(ratios))
        
        ax.errorbar(x_values, ratios, yerr=np.maximum(errors, 0.0), marker=markers[i], linestyle='None', 
                    label=epsilon_label, color=colors[i], markeredgecolor=colors[i], 
                    markerfacecolor='none', capsize=2)

        x_len = len(x_values)
        
    # Add vertical lines every NumPhiBins
    for x in range(0, x_len, NumPhiBins):
        ax.axvline(x, color='blue', linestyle='-', linewidth=0.75, alpha=0.5)        
        
    ax.axhline(1.0, color='gray', linestyle='--')
    ax.set_xlabel('t-$\phi$ bin', fontsize=24)
    ax.set_ylabel('Ratio', fontsize=24)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.set_ylim(0.0, 2.0)
    ax.legend(fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')

    # Set integer ticks on x-axis
    ax.set_xticks(range(0, x_len, 2))
    ax.set_xticklabels(range(1, x_len + 1, 2))  # Start from 1 instead of 0

    # Add grid
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')

    ### HERE 1

    a_hi_lst  = []
    b_hi_lst  = []
    c_hi_lst  = []
    d_hi_lst  = []
    
    a_lo_lst  = []
    b_lo_lst  = []
    c_lo_lst  = []
    d_lo_lst  = []
    
    def fit_function(Wset, Q2set, a, b, c, d):
        Wval = np.linspace(float(W.replace("p","."))-0.5, float(W.replace("p","."))+0.5, len(Wset))
        Q2val = np.linspace(float(Q2.replace("p","."))-0.5, float(Q2.replace("p","."))+0.5, len(Q2set))
        return 1 + b*(Wval-Wset) + c*(Q2val-Q2set) + d*(Wval-Wset)*(Q2val-Q2set)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        j=0
        # Create a single figure and axis object for all phi bins
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title("t={:.3f}, $Q^2$={:.1f}, W={:.2f}".format(t_bin_centers[k], float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)        
        # Loop through t bins and plot data
        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            if "hi" in df_key:
                epsilon_label = "High $\epsilon$" if k == 0 else ""
                epsilon_fit_color = "r-"
            else:
                epsilon_label = "Low $\epsilon$" if k == 0 else ""
                epsilon_fit_color = "-"

            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])
                
            ratios = df['x_real'][mask]/df['x_mod'][mask]
            errors = df['dx_real'][mask]/df['x_mod'][mask]
            non_zero_mask = (ratios != 0) & (errors != 0)
            ratios = ratios[non_zero_mask]
            errors = errors[non_zero_mask]

            x_increment = j+k*NumPhiBins
            
            # Use x_increment for x-axis values
            x_values = np.arange(x_increment, x_increment+len(ratios))

            ax.errorbar(x_values, ratios, yerr=np.maximum(errors, 0.0), marker=markers[i], linestyle='None', 
                        label=epsilon_label, color=colors[i], markeredgecolor=colors[i], 
                        markerfacecolor='none', capsize=2)

            def fit_func(data, a, b, c, d):
                Wval, Q2val = data
                return fit_function(Wval, Q2val, a, b, c, d)

            try:
                popt, pcov = curve_fit(fit_func, (df['W'][mask][non_zero_mask], df['Q2'][mask][non_zero_mask]), ratios, sigma=errors, absolute_sigma=True)
                
                a_fit, b_fit, c_fit, d_fit = popt

                fitted_values = fit_function(df['W'][mask][non_zero_mask], df['Q2'][mask][non_zero_mask], a_fit, b_fit, c_fit, d_fit)

                # Plot fitted function
                ax.plot(range(x_increment, x_increment+len(ratios)), fitted_values, epsilon_fit_color, label=f'a = {a_fit:.4f}\nb = {b_fit:.4f}\nc = {c_fit:.4f}\nd = {d_fit:.4f}')
            except (TypeError, ValueError) as e:
                print("Fit not found!")
                continue

            x_len = x_increment+len(x_values)

            if "hi" in df_key:
                a_hi_lst.append((df['t'][x_increment], a_fit))
                b_hi_lst.append((df['t'][x_increment], b_fit))
                c_hi_lst.append((df['t'][x_increment], c_fit))
                d_hi_lst.append((df['t'][x_increment], d_fit))
            else:
                a_lo_lst.append((df['t'][x_increment], a_fit))
                b_lo_lst.append((df['t'][x_increment], b_fit))
                c_lo_lst.append((df['t'][x_increment], c_fit))
                d_lo_lst.append((df['t'][x_increment], d_fit))            
                
        # Add the equation as text above the legend
        equation = r'$a + b\cdot(W - W_{\text{c}}) + c\cdot(Q^2 - Q^2_{\text{c}}) + d\cdot(W - W_{\text{c}}) (Q^2 - Q^2_{\text{c}})$'
        ax.text(1.05, 1.02, equation, transform=ax.transAxes, fontsize=10, verticalalignment='bottom')

        ax.axhline(1.0, color='gray', linestyle='--')
        ax.set_xlabel('t-$\phi$ bin', fontsize=24)
        ax.set_ylabel('Ratio', fontsize=24)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.legend(fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')

        # Set integer ticks on x-axis
        ax.set_xticks(range(x_increment, x_len, 2))
        ax.set_xticklabels(range(x_increment+1, x_len + 1, 2))  # Start from 1 instead of 0

        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')

        j+=1

    '''
    # Create a figure with 4 subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    # Scatter data for 'a'
    axs[0, 0].scatter(*zip(*a_hi_lst), marker=markers[0], linestyle='None', label='High $\epsilon$', color=colors[0])
    axs[0, 0].scatter(*zip(*a_lo_lst), marker=markers[1], linestyle='None', label='Low $\epsilon$', color=colors[1])
    axs[0, 0].set_title('a')
    axs[0, 0].legend()

    # Scatter data for 'b'
    axs[0, 1].scatter(*zip(*b_hi_lst), marker=markers[0], linestyle='None', label='High $\epsilon$', color=colors[0])
    axs[0, 1].scatter(*zip(*b_lo_lst), marker=markers[1], linestyle='None', label='Low $\epsilon$', color=colors[1])
    axs[0, 1].set_title('$b\cdot(W - W_{\text{c}})$')
    #axs[0, 1].legend()

    # Scatter data for 'c'
    axs[1, 0].scatter(*zip(*c_hi_lst), marker=markers[0], linestyle='None', label='High $\epsilon$', color=colors[0])
    axs[1, 0].scatter(*zip(*c_lo_lst), marker=markers[1], linestyle='None', label='Low $\epsilon$', color=colors[1])
    axs[1, 0].set_title('$c\cdot(Q^2 - Q^2_{\text{c}})$')
    #axs[1, 0].legend()

    # Scatter data for 'd'
    axs[1, 1].scatter(*zip(*d_hi_lst), marker=markers[0], linestyle='None', label='High $\epsilon$', color=colors[0])
    axs[1, 1].scatter(*zip(*d_lo_lst), marker=markers[1], linestyle='None', label='Low $\epsilon$', color=colors[1])
    axs[1, 1].set_title('$d\cdot(W - W_{\text{c}}) (Q^2 - Q^2_{\text{c}})$')
    #axs[1, 1].legend()

    # Adjust layout
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')

    ### HERE 2

    a_hi_lst  = []
    b_hi_lst  = []
    c_hi_lst  = []
    d_hi_lst  = []
    
    a_lo_lst  = []
    b_lo_lst  = []
    c_lo_lst  = []
    d_lo_lst  = []
    
    def fit_function(phival, thetaval, a, b, c, d):
        #phival = np.linspace(0.0, 360, len(thetaval)) 
        return 1 + b*(np.sin(thetaval)**2) + c*(np.sin(thetaval)*np.cos(phival)) + d*((np.sin(thetaval)**2)*np.cos(2*phival))

    # Loop through t bins and plot data
    for k in range(NumtBins):
        j=0
        # Create a single figure and axis object for all phi bins
        fig, ax = plt.subplots(figsize=(12, 8))
        ax.set_title("t={:.3f}, $Q^2$={:.1f}, W={:.2f}".format(t_bin_centers[k], float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)
        # Loop through t bins and plot data
        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            if "hi" in df_key:
                epsilon_label = "High $\epsilon$" if k == 0 else ""
                epsilon_fit_color = "r-"
            else:
                epsilon_label = "Low $\epsilon$" if k == 0 else ""
                epsilon_fit_color = "-"

            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])

            ratios = df['x_real'][mask]/df['x_mod'][mask]
            errors = df['dx_real'][mask]/df['x_mod'][mask]
            non_zero_mask = (ratios != 0) & (errors != 0)
            ratios = ratios[non_zero_mask]
            errors = errors[non_zero_mask]        

            x_increment = j+k*NumPhiBins
            
            # Use x_increment for x-axis values
            x_values = np.arange(x_increment, x_increment+len(ratios))
            
            ax.errorbar(x_values, ratios, yerr=np.maximum(errors, 0.0), marker=markers[i], linestyle='None', 
                        label=epsilon_label, color=colors[i], markeredgecolor=colors[i], 
                        markerfacecolor='none', capsize=2)

            def fit_func(data, a, b, c, d):
                phival, thetaval = data
                return fit_function(phival, thetaval, a, b, c, d)

            try:
                popt, pcov = curve_fit(fit_func, (df['phi'][mask][non_zero_mask].to_numpy(), df['th_cm'][mask][non_zero_mask].to_numpy()), ratios, sigma=errors, absolute_sigma=True)

                a_fit, b_fit, c_fit, d_fit = popt
                
                fitted_values = fit_function(df['phi'][mask][non_zero_mask], df['th_cm'][mask][non_zero_mask], a_fit, b_fit, c_fit, d_fit)

                # Plot fitted function
                ax.plot(range(x_increment, x_increment+len(ratios)), fitted_values, epsilon_fit_color, label=f'a = {a_fit:.4f}\nb = {b_fit:.4f}\nc = {c_fit:.4f}\nd = {d_fit:.4f}')
                
            except (TypeError, ValueError) as e:
                print("Fit not found!")
                continue
            
            x_len = x_increment+len(x_values)

            if "hi" in df_key:
                a_hi_lst.append((df['t'][x_increment], a_fit))
                b_hi_lst.append((df['t'][x_increment], b_fit))
                c_hi_lst.append((df['t'][x_increment], c_fit))
                d_hi_lst.append((df['t'][x_increment], d_fit))
            else:
                a_lo_lst.append((df['t'][x_increment], a_fit))
                b_lo_lst.append((df['t'][x_increment], b_fit))
                c_lo_lst.append((df['t'][x_increment], c_fit))
                d_lo_lst.append((df['t'][x_increment], d_fit))
            
        # Add the equation as text above the legend
        equation = r'$a + b\cdot\sin^2(\theta) + c\cdot\sin(\theta) \cos(\phi) + d\cdot\sin^2(\theta) \cos(2\phi)$'
        ax.text(1.05, 1.02, equation, transform=ax.transAxes, fontsize=10, verticalalignment='bottom')

        ax.axhline(1.0, color='gray', linestyle='--')
        ax.set_xlabel('t-$\phi$ bin', fontsize=24)
        ax.set_ylabel('Ratio', fontsize=24)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.legend(fontsize=10, bbox_to_anchor=(1.05, 1), loc='upper left')

        # Set integer ticks on x-axis
        ax.set_xticks(range(x_increment, x_len, 2))
        ax.set_xticklabels(range(x_increment+1, x_len + 1, 2))  # Start from 1 instead of 0

        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')

        j+=1

    # Create a figure with 4 subplots
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))

    # Scatter data for 'a'
    axs[0, 0].scatter(*zip(*a_hi_lst), marker=markers[0], linestyle='None', label='High $\epsilon$', color=colors[0])
    axs[0, 0].scatter(*zip(*a_lo_lst), marker=markers[1], linestyle='None', label='Low $\epsilon$', color=colors[1])
    axs[0, 0].set_title('a')
    axs[0, 0].legend()

    # Scatter data for 'b'
    axs[0, 1].scatter(*zip(*b_hi_lst), marker=markers[0], linestyle='None', label='High $\epsilon$', color=colors[0])
    axs[0, 1].scatter(*zip(*b_lo_lst), marker=markers[1], linestyle='None', label='Low $\epsilon$', color=colors[1])
    axs[0, 1].set_title('$b\cdot\sin^2(\theta)$')
    #axs[0, 1].legend()

    # Scatter data for 'c'
    axs[1, 0].scatter(*zip(*c_hi_lst), marker=markers[0], linestyle='None', label='High $\epsilon$', color=colors[0])
    axs[1, 0].scatter(*zip(*c_lo_lst), marker=markers[1], linestyle='None', label='Low $\epsilon$', color=colors[1])
    axs[1, 0].set_title('$c\cdot\sin(\theta) \cos(\phi)$')
    #axs[1, 0].legend()

    # Scatter data for 'd'
    axs[1, 1].scatter(*zip(*d_hi_lst), marker=markers[0], linestyle='None', label='High $\epsilon$', color=colors[0])
    axs[1, 1].scatter(*zip(*d_lo_lst), marker=markers[1], linestyle='None', label='Low $\epsilon$', color=colors[1])
    axs[1, 1].set_title('$d\cdot\sin^2(\theta) \cos(2\phi)$')
    #axs[1, 1].legend()

    # Adjust layout
    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
        
    ##
                
    
        
    ###

    # Define exponential function
    def exp_func(t, a, b):
        return a * np.exp(b * t)

    # Create a figure and axis objects for Q2 plot
    fig, axes = plt.subplots(1, 1, figsize=(12, 8), sharex=True)

    ax = axes
    ax.set_title("$Q^2$={:.1f}, W={:.2f}".format(float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

    for i, df_key in enumerate(['kindata_loeps_{}'.format('center'), 'kindata_hieps_{}'.format('center')]):
        df = file_df_dict[df_key]
        if "hi" in df_key:
            df_key = "High $\epsilon$"
        else:
            df_key = "Low $\epsilon$"

        ax.errorbar(t_bin_centers, df['Q2'], yerr=np.maximum(df['dQ2'], 0.0), marker=markers[i], linestyle='None', label=df_key, color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)

        try:
            # Fit the data using exponential function
            popt, _ = curve_fit(exp_func, t_bin_centers, df['Q2'])
            fit_line = exp_func(t_bin_centers, *popt)
            ax.plot(t_bin_centers, fit_line, linestyle='-', color=colors[i], label="{0} Fit: Q(t) = {1:.2f}e^({2:.2f}t)".format(df_key, popt[0], popt[1]))
        except (TypeError, ValueError) as e:
            print("Fit not found!")
            continue

    ax.set_xlabel('-t', fontsize=24)
    ax.set_ylabel('$Q^2$', fontsize=24)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)        
    ax.set_xlim(tmin-0.1, tmax+0.1)
    ax.legend(fontsize=16)
    # Add grid
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')

    # Create a figure and axis objects for W plot
    fig, axes = plt.subplots(1, 1, figsize=(12, 8), sharex=True)

    ax = axes
    ax.set_title("$Q^2$={:.1f}, W={:.2f}".format(float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

    for i, df_key in enumerate(['kindata_loeps_{}'.format('center'), 'kindata_hieps_{}'.format('center')]):
        df = file_df_dict[df_key]
        if "hi" in df_key:
            df_key = "High $\epsilon$"
        else:
            df_key = "Low $\epsilon$"

        ax.errorbar(t_bin_centers, df['W'], yerr=np.maximum(df['dW'], 0.0), marker=markers[i], linestyle='None', label=df_key, color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)

        try:
            # Fit the data using exponential function
            popt, _ = curve_fit(exp_func, t_bin_centers, df['W'])
            fit_line = exp_func(t_bin_centers, *popt)
            ax.plot(t_bin_centers, fit_line, linestyle='-', color=colors[i], label="{0} Fit: W(t) = {1:.2f}e^({2:.2f}t)".format(df_key, popt[0], popt[1]))
        except (TypeError, ValueError) as e:
            print("Fit not found!")
            continue
            
    ax.set_xlabel('-t', fontsize=24)
    ax.set_ylabel('W', fontsize=24)
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)        
    ax.set_xlim(tmin-0.1, tmax+0.1)
    ax.legend(fontsize=16)
    # Add grid
    ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout()
    pdf.savefig(fig, bbox_inches='tight')
        
    ###
        
    # Create a figure and axis objects
    fig, axes = plt.subplots(NumtBins, 1, figsize=(12, 8), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t={:.3f}, $Q^2$={:.1f}, W={:.2f}".format(t_bin_centers[k], float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            if "hi" in df_key:
                df_key = "High $\epsilon$"
            else:
                df_key = "Low $\epsilon$"
                
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t']) & (df['x_real'] > 0.0)  & (df['dx_real'] > 0.0)
            ax.errorbar(df['th_cm'][mask], df['x_real'][mask], yerr=np.maximum(df['dx_real'][mask], 0.0), marker=markers[i], linestyle='None', label=df_key, color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)

        ax.set_xlabel('$\theta_{cm}$')
        ax.set_ylabel('x_real')
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.legend(fontsize=14)
        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')

    fig, axes = plt.subplots(NumtBins, 1, figsize=(12, 8), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t={:.3f}, $Q^2$={:.1f}, W={:.2f}".format(t_bin_centers[k], float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            if "hi" in df_key:
                df_key = "High $\epsilon$"
            else:
                df_key = "Low $\epsilon$"
                
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])
            ax.scatter(df['th_cm'][mask], df['x_mod'][mask], marker=markers[i+2], linestyle='None', label=df_key, color=colors[i])

        ax.set_xlabel('$\theta_{cm}$')
        ax.set_ylabel('x_mod')
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.legend(fontsize=14)
        # Add grid
        ax.grid(True, which='both', linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    
    fig, axes = plt.subplots(NumtBins, 1, figsize=(12, 8), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t={:.3f}, $Q^2$={:.1f}, W={:.2f}".format(t_bin_centers[k], float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            if "hi" in df_key:
                df_key = "High $\epsilon$"
            else:
                df_key = "Low $\epsilon$"
                
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t']) & (df['x_real'] > 0.0)  & (df['dx_real'] > 0.0)
            ax.errorbar(df['phi'][mask], df['x_real'][mask], yerr=np.maximum(df['dx_real'][mask], 0.0), marker=markers[i], linestyle='None', label=df_key, color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)

        ax.set_xlabel('$\phi$', fontsize=24)
        ax.set_ylabel('x_real', fontsize=24)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.set_xlim(-5, 365)
        ax.legend(fontsize=14)
        # Add grid to subplot
        ax.grid(True, linestyle='--', linewidth=0.5)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    
    fig, axes = plt.subplots(NumtBins, 1, figsize=(12, 8), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t={:.3f}, $Q^2$={:.1f}, W={:.2f}".format(t_bin_centers[k], float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            if "hi" in df_key:
                df_key = "High $\epsilon$"
            else:
                df_key = "Low $\epsilon$"
                
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t'])
            ax.scatter(df['phi'][mask], df['x_mod'][mask], marker=markers[i+2], linestyle='None', label=df_key, color=colors[i])

        ax.set_xlabel('$\phi$', fontsize=24)
        ax.set_ylabel('x_mod', fontsize=24)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.set_xlim(-5, 365)
        ax.legend(fontsize=14)
        # Add grid to subplot
        ax.grid(True, linestyle='--', linewidth=0.5)
        
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')
    
    # Loop through t bins and plot data
    for k in range(NumtBins):

        # Create a figure and axis objects
        fig, axes = plt.subplots(1, 1, figsize=(12, 8), sharex=True)

        ax = axes
        ax.set_title("t={:.3f}, $Q^2$={:.1f}, W={:.2f}".format(t_bin_centers[k], float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            if "hi" in df_key:
                df_key = "High $\epsilon$"
            else:
                df_key = "Low $\epsilon$"
                
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t']) & (df['x_real'] > 0.0)  & (df['dx_real'] > 0.0)
            ax.errorbar(df['phi'][mask], df['x_real'][mask], yerr=np.maximum(df['dx_real'][mask], 0.0), marker=markers[i], linestyle='None', label=df_key, color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)
            ax.scatter(df['phi'][mask], df['x_mod'][mask], marker=markers[i+2], linestyle='None', label=df_key+" Model")

        ax.set_xlabel('$\phi$', fontsize=24)
        ax.set_ylabel('x_real', fontsize=24)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.set_xlim(-5, 365)
        ax.legend(fontsize=14)
        # Add grid to subplot
        ax.grid(True, linestyle='--', linewidth=0.5)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        pdf.savefig(fig, bbox_inches='tight')
    
    fig, axes = plt.subplots(NumtBins, 1, figsize=(12, 8), sharex=True)

    # Loop through t bins and plot data
    for k in range(NumtBins):
        ax = axes[k]
        ax.set_title("t={:.3f}, $Q^2$={:.1f}, W={:.2f}".format(t_bin_centers[k], float(Q2.replace("p",".")), float(W.replace("p","."))), fontsize=24)

        for i, df_key in enumerate(['unsep_file_loeps', 'unsep_file_hieps']):
            df = file_df_dict[df_key]
            if "hi" in df_key:
                df_key = "High $\epsilon$"
            else:
                df_key = "Low $\epsilon$"
                
            mask =  (df['t'][k*NumPhiBins+int(i/NumPhiBins)] == df['t']) & (df['x_real'] > 0.0)  & (df['dx_real'] > 0.0)
            ax.errorbar(df['phi'][mask], df['x_real'][mask], yerr=np.maximum(df['dx_real'][mask], 0.0), marker=markers[i], linestyle='None', label=df_key, color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)
            ax.scatter(df['phi'][mask], df['x_mod'][mask], marker=markers[i+2], linestyle='None', label=df_key+" Model")

        ax.set_xlabel('$\phi$', fontsize=24)
        ax.set_ylabel('x_real', fontsize=24)
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.set_xlim(-5, 365)
        ax.legend(fontsize=14)
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
                inp_param = '{} {} {} {} {} {} '.format(Q2.replace("p","."), W.replace("p","."), row['th_cm'], row['t'], row['Q2'], row['W'])+' '.join(param_arr)
                model.append(import_model(sig, inp_param))
            # Check that model sig is not all zeros
            if not all(element == 0 for element in model):
                ax.plot(df['t'], model, linestyle='-.', color='red', label='Model Fit')
            ax.errorbar(df['t'], df['{}'.format(sig)], yerr=np.maximum(df['d{}'.format(sig)], 0.0), marker=markers[i], linestyle='None', label='Data', color=colors[i], markeredgecolor=colors[i], markerfacecolor='none', capsize=2)
        ax.set_xlabel('t')
        ax.set_ylabel("${}$".format(formatted_sig))
        ax.tick_params(axis='x', labelsize=16)
        ax.tick_params(axis='y', labelsize=16)        
        ax.set_xlim(tmin-0.1, tmax+0.1)
        ax.legend(fontsize=14)
        # Add grid to subplot
        ax.grid(True, linestyle='--', linewidth=0.5)
    print("="*50)
        
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches='tight')

    '''
