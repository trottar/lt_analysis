#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-11 06:20:02 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import os, sys, re

'''
# Define the function f(x; params)
def f(x, params):
    """
    Generalized function for multiple parameters.
    params: List of parameters [a1, a2, ..., an].
    """
    return sum(p * x**i for i, p in enumerate(params))

# Mock analysis procedure A
def analysis_procedure(data, params):
    """
    Update parameters based on the data and current parameterization.
    """
    new_params = [p + 0.01 * (-1)**i * np.mean(data) for i, p in enumerate(params)]
    return new_params

# Iterative process
def iterative_process(initial_params, initial_data, iterations=10):
    params = initial_params
    data = initial_data
    
    params_values, data_values = [params], [data]
    
    for _ in range(iterations):
        # Apply the analysis procedure to update parameters
        params = analysis_procedure(data, params)
        params_values.append(params)
        
        # Generate new data using the updated function
        data = f(data, params)
        data_values.append(data)
    
    return params_values, data_values
'''

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

##################################################################################################################################################
# Importing utility functions

sys.path.append("../utility")
from utility import show_pdf_with_evince

###############################################################################################################################################

def compare_iters(pol_str, ParticleType, Q2, W, LOEPS, HIEPS):


    ################################################################################################################################################
    '''
    Import separated xsects model
    '''
    
    sys.path.append("models")
    if pol_str == "pl" and ParticleType == "kaon":
        from sep_xsect_kaon_pl import import_model

    ################################################################################################################################################        
    
    f_iter = "{}/{}_Q{}W{}_iter.dat".format(LTANAPATH,ParticleType,Q2,W)

    inp_dir = f_iter.replace(f"{ParticleType}_Q{Q2}W{W}_iter.dat", f"src/{ParticleType}")

    ################################################################################################################################################
    # Import iteration directory...
    # Grab functional forms and parameterization for each iteration directory
    # 

    # Open the file, read the lines, and filter them
    with open(f_iter, "r") as infile:
        lines = infile.readlines()    

    iter_arr = [f.rstrip('\n') for f in lines]
    iter_start = iter_arr[0]
    iter_end = iter_arr[-1]
    iterations = len(lines)-1

    OutFilename = f"compare_iters_{iter_start}-{iter_end}"
    outputpdf  = OUTPATH+"/" + OutFilename + ".pdf"

    # Create settings dictionary using a loop
    settings = {}

    for i, date in enumerate(iter_arr):
        ###############################################################################################################################################
        '''
        Import model parameterization
        '''

        inp_dir = f"{TEMP_CACHEPATH}/{ParticleType}/Q{Q2}W{W}/{date}"
        
        param_file = '{}/parameters/par.{}_Q{}W{}.dat'.format(inp_dir, pol_str, Q2.replace("p",""), W.replace("p",""))

        ###############################################################################################################################################
        
        param_arr = []
        
        with open(param_file, 'r') as f:
            for j, line in enumerate(f):
                columns = line.split()
                param_arr.append(str(columns[0]))    
                
        settings[f'set_{i+1}'] = {
            'Q2': Q2,
            'W': W,
            'LOEPS': LOEPS,
            'HIEPS': HIEPS,
            'params' : param_arr,
            'date': date
        }

    comb_dict = {}

    file_df_dict = {}
    
    # Extract values for specific sets
    for key, values in settings.items():

        # Unpack values into variables
        Q2, W, LOEPS, HIEPS, param_arr, date = values.values()

        inp_dir = f"{TEMP_CACHEPATH}/{ParticleType}/Q{Q2}W{W}/{date}"

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
                
                file_df_dict["params"] = pd.DataFrame(param_arr, columns=["params"])
                file_df_dict["date"] = pd.DataFrame([date], columns=["date"])
                
            comb_dict[f'Q{Q2}W{W}'] = file_df_dict

    # Create a structured dictionary to track iterations explicitly
    iteration_data = {}

    # Iterate through the settings to populate iteration_data
    for key, values in settings.items():
        # Unpack values into variables
        Q2, W, LOEPS, HIEPS, param_arr, date = values.values()

        # Create a structured entry for each iteration
        iteration_data[date] = {
            'params': param_arr,
            'Q2': Q2,
            'W': W,
            'sep_file': comb_dict[f'Q{Q2}W{W}']['sep_file'],
            'sigL': comb_dict[f'Q{Q2}W{W}']['sep_file']['sigL'].tolist(),
            'sigT': comb_dict[f'Q{Q2}W{W}']['sep_file']['sigT'].tolist(),
            'sigLT': comb_dict[f'Q{Q2}W{W}']['sep_file']['sigLT'].tolist(),
            'sigTT': comb_dict[f'Q{Q2}W{W}']['sep_file']['sigTT'].tolist()
        }

    print("\n\nIteration Data:")
    print(iteration_data)            
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
    
    # Redefine the plotting section
    for k, sig in enumerate(['sigL','sigT','sigLT','sigTT']):
        
        # Lists to store iteration-specific data
        params_values = []
        data_values = []
        dates = []
        
        tmp_file_name = outputpdf.replace(OutFilename, f"{OutFilename}_{sig}")

        # Collect data from iteration_data
        for date, data in iteration_data.items():
            dates.append(date)
            params_values.append(data['params'])
            data_values.append(data[sig])
        
        # Create a PdfPages object to manage the PDF file
        with PdfPages(tmp_file_name) as pdf:    
            # Collect data and parameters for each iteration
            for iter_key, iter_data in comb_dict.items():

                df = iter_data["sep_file"]

                date_df = iter_data["date"]
                param_df  = iter_data["params"]
                                        
                # 1. Parameter Evolution Plot
                fig = plt.figure(figsize=(12, 6))
                for i in range(len(params_values[0][:k+4])): # Maximum of 4 parameters
                    param_evolution = [params[i] for params in params_values]
                    plt.plot(dates, param_evolution, label=f'Parameter {i}', marker='o')

                plt.xlabel('Iteration Date')
                plt.ylabel('Parameter Value')
                plt.title(f'Parameter Evolution Across Iterations for {sig}')
                plt.xticks(rotation=45, ha='right')
                plt.legend()
                plt.tight_layout()
                plt.grid(True)
                pdf.savefig(fig, bbox_inches='tight')

                # 2. Data Distribution Plot
                fig = plt.figure(figsize=(12, 6))
                for i, (data, date) in enumerate(zip(data_values, dates)):
                    plt.hist(data, bins=20, alpha=0.5, label=f'Iteration {date}', density=True)

                plt.xlabel(f'{sig} Values')
                plt.ylabel('Density')
                plt.title(f'Data Distribution Across Iterations for {sig}')
                plt.legend()
                plt.tight_layout()
                plt.grid(True)
                pdf.savefig(fig, bbox_inches='tight')

                # 3. Model Comparison Plot
                fig = plt.figure(figsize=(12, 6))
                model_values = []

                for j, (df_row, date, params) in enumerate(zip(comb_dict.values(), dates, params_values)):
                    df = df_row["sep_file"]
                    model_for_iteration = []

                    for _, row in df.iterrows():
                        inp_param = '{} {} {} {} {} {} '.format(
                            Q2.replace("p","."), 
                            W.replace("p","."), 
                            row['th_cm'], 
                            row['t'], 
                            row['Q2'], 
                            row['W']
                        ) + ' '.join(map(str, params))

                        model_for_iteration.append(import_model(sig, inp_param))

                    model_values.append(model_for_iteration)

                    plt.plot(df['t'], model_for_iteration, label=f'Iteration {date}', marker='o')

                plt.xlabel('t')
                plt.ylabel(f'{sig} Model Value')
                plt.title(f'Model Evolution Across Iterations for {sig}')
                plt.legend()
                plt.tight_layout()
                plt.grid(True)
                pdf.savefig(fig, bbox_inches='tight')

                # 4. Residuals and Convergence Plot
                fig = plt.figure(figsize=(12, 6))
                residuals = []

                for data, model in zip(data_values, model_values):
                    # Calculate residuals (difference between data and model)
                    res = np.abs(data - model)
                    residuals.append(np.mean(res))

                plt.plot(dates, residuals, marker='o', color='purple')
                plt.xlabel('Iteration Date')
                plt.ylabel('Mean Absolute Residual')
                plt.title(f'Convergence of Residuals Across Iterations for {sig}')
                plt.xticks(rotation=45, ha='right')
                plt.tight_layout()
                plt.grid(True)
                pdf.savefig(fig, bbox_inches='tight')

                # 5. Scatter Plot of Data vs Model
                fig = plt.figure(figsize=(12, 6))
                for i, (data, model, date) in enumerate(zip(data_values, model_values, dates)):
                    plt.scatter(data, model, label=f'Iteration {date}', alpha=0.7)

                plt.plot([data.min(), data.max()], [data.min(), data.max()], 'r--', label='Ideal Fit')
                plt.xlabel('Observed Data')
                plt.ylabel('Model Prediction')
                plt.title(f'Data vs Model Predictions for {sig}')
                plt.legend()
                plt.tight_layout()
                plt.grid(True)
                pdf.savefig(fig, bbox_inches='tight')

                # Open the PDF
                show_pdf_with_evince(tmp_file_name)
