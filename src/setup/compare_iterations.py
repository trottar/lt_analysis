#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-12-11 04:15:38 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os, sys

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

def compare_iters(ParticleType,Q2,W):

    f_iter = "{}/{}_Q{}W{}_iter.dat".format(LTANAPATH,ParticleType,Q2,W)

    inp_dir = f_iter.replace(f"{ParticleType}_Q{Q2}W{W}_iter.dat", "src/{ParticleType}")

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

    ################################################################################################################################################
    # Import iteration directory...
    # Grab functional forms and parameterization for each iteration directory
    # 

    # Open the file, read the lines, and filter them
    with open(f_iter, "r") as infile:
        lines = infile.readlines()    

    iter_arr = [f for f in f_iter]
    iterations = len(lines)

    OutFilename = f"compare_iters_{iter_start}-{iter_end}"
    outputpdf  = OUTPATH+"/" + OutFilename + ".pdf"

    # Create settings dictionary using a loop
    settings = {}

    for i, date in enumerate(iterations):
        settings[f'set_{i+1}'] = {
            'inp_dir': date
        }

    comb_dict = {}

    # Extract values for specific sets
    for key, values in settings.items():

        # Unpack values into variables
        Q2, W, LOEPS, HIEPS, inp_dir = values.values()

        # Add the extracted values into comb_dict
        comb_dict[key] = {
            'inp_dir': inp_dir
        }

        inp_dir = TEMP_CACHEPATH+"/{}/Q{}W{}/".format(ParticleType,Q2,W)+inp_dir

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

    for val in file_df_dict['sep_file']:

        for index, data in val.iterrows():

            tmp_file_name = outputpdf.replace(OutFilename, f"{OutFilename}_{index}")
            # Create a PdfPages object to manage the PDF file
            with PdfPages(tmp_file_name) as pdf:

                params_values = param_arr
                data_values = data

                # 1. Parameter Evolution Plot
                # Purpose: Track how parameters change across iterations.
                fig = plt.figure(figsize=(10, 6))
                for i in range(len(initial_params)):
                    plt.plot(range(iterations + 1), [p[i] for p in params_values], label=f'Parameter {i}', marker='o')
                plt.xlabel('Iteration')
                plt.ylabel('Parameter Value')
                plt.title('Parameter Evolution Across Iterations')
                plt.legend()
                plt.grid()
                pdf.savefig(fig, bbox_inches='tight')

                # 2. Function Evolution Plot
                # Purpose: Visualize how f(x) changes over iterations.
                x_values = np.linspace(0, 2, 100)  # Example x range
                fig = plt.figure(figsize=(10, 6))
                for i, params in enumerate(params_values):
                    plt.plot(x_values, f(x_values, params), label=f'Iteration {i}')
                plt.xlabel('x')
                plt.ylabel('f(x)')
                plt.title('Function Evolution Across Iterations')
                plt.legend()
                plt.grid()
                pdf.savefig(fig, bbox_inches='tight')

                # 3. Data Distribution Plot
                # Purpose: Observe how the data distribution changes with updated parameterization.
                fig = plt.figure(figsize=(10, 6))
                for i, data in enumerate(data_values):
                    plt.hist(data, bins=20, alpha=0.5, label=f'Iteration {i}', density=True)
                plt.xlabel('Data Values')
                plt.ylabel('Density')
                plt.title('Data Distribution Across Iterations')
                plt.legend()
                plt.grid()
                pdf.savefig(fig, bbox_inches='tight')

                # 4. Parameter vs. Output Data Plot
                # Purpose: Correlate parameters with properties of the generated data.
                data_means = [np.mean(data) for data in data_values]
                fig = plt.figure(figsize=(10, 6))
                for i in range(len(initial_params)):
                    plt.scatter([p[i] for p in params_values], data_means, label=f'Param {i} vs Data Mean', marker='o')
                plt.xlabel('Parameter Value')
                plt.ylabel('Data Mean')
                plt.title('Parameter vs. Data Property')
                plt.legend()
                plt.grid()
                pdf.savefig(fig, bbox_inches='tight')

                # 5. Convergence Plot
                # Purpose: Show how a chosen metric (e.g., error or residuals) evolves toward convergence.
                residuals = [np.sum((data - f(data, params))**2) for data, params in zip(data_values, params_values)]
                fig = plt.figure(figsize=(10, 6))
                plt.plot(range(iterations + 1), residuals, marker='o', color='purple')
                plt.xlabel('Iteration')
                plt.ylabel('Residual Sum of Squares')
                plt.title('Convergence of Residuals Across Iterations')
                plt.grid()
                pdf.savefig(fig, bbox_inches='tight')

                # 6. Heatmap of Function Changes
                # Purpose: Provide a comprehensive view of f(x) across x-values and iterations.
                z_values = np.array([f(x_values, params) for params in params_values])
                fig = plt.figure(figsize=(10, 6))
                plt.imshow(z_values, aspect='auto', cmap='viridis', extent=[0, 2, 0, iterations])
                plt.colorbar(label='f(x)')
                plt.xlabel('x')
                plt.ylabel('Iteration')
                plt.title('Heatmap of Function Changes Across Iterations')
                pdf.savefig(fig, bbox_inches='tight')

                show_pdf_with_evince(tmp_file_name)
