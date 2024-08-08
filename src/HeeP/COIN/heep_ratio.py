#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-08-08 16:25:55 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

# Settings
Q2 = [7.52, 6.66, 4.92, 3.40, 3.40]
Ebeam = [10.59070, 8.21290, 6.18720, 4.92950, 3.8338]

# Data W(0.7,1.1)
#y_data_nooffset = [5.15e1, 9.15e1, 2.91e2, 1.35e3, 8.75e2]
#y_simc_nooffset = [5.56e1, 8.32e1, 2.62e2, 1.30e3, 9.04e2]
#rel_yield_nooffset = [9.25e-1, 1.10, 1.11, 1.04, 9.68e-1]
#y_data = [6.50e1, 9.00e1, 2.89e2, 1.34e3, 8.71e2]
#y_simc = [5.01e1, 8.32e1, 2.62e2, 1.30e3 9.04e2]
#rel_yield = [1.17, 1.08, 1.10, 1.03, 9.64e-1]
# Data W(0.5,2.0)
#y_data_nooffset = [6.50e1, 1.07e2, 3.23e2, 1.46e3, 9.11e2]
#y_simc_nooffset = [5.56e1, 8.32e1, 3.05e2, 1.37e3, 9.04e2]
#rel_yield_nooffset = [1.17, 1.29, 1.06, 1.06, 1.01]
#y_data = [6.50e1, 1.07e2, 3.23e2, 1.46e3, 9.11e2]
#y_simc = [5.57e1, 8.32e1, 3.06e2, 1.41e3, 9.04e2]
#rel_yield = [1.17, 1.29, 1.05, 1.03, 1.01]
# RLT (7/29/2024): Integrated fit
y_data = [1.45e+03, 9.08e+02]
y_simc = [1.41e+03, 9.04e+02]
rel_yield = y_simc/y_data

# Calculate linear fit for y_data_nooffset
#slope_data_nooffset, intercept_data_nooffset, _, _, _ = linregress(Q2, y_data_nooffset)

# Calculate linear fit for y_simc_nooffset
#slope_simc_nooffset, intercept_simc_nooffset, _, _, _ = linregress(Q2, y_simc_nooffset)

# Calculate linear fit for y_data
slope_data, intercept_data, _, _, _ = linregress(Q2, y_data)

# Calculate linear fit for y_simc
slope_simc, intercept_simc, _, _, _ = linregress(Q2, y_simc)

#slope_relyield_nooffset, intercept_relyield_nooffset, _, _, _ = linregress(Q2, rel_yield_nooffset)

slope_relyield, intercept_relyield, _, _, _ = linregress(Q2, rel_yield)

# Plotting
plt.figure(figsize=(12,8))

'''
# Plot y_data_nooffset and y_simc_nooffset on the same plot with linear fits
plt.subplot(221)
plt.plot(Q2, y_data_nooffset, 'bo', label='y_data_nooffset')
plt.plot(Q2, y_simc_nooffset, 'ro', label='y_simc_nooffset')
plt.plot(Q2, slope_data_nooffset * np.array(Q2) + intercept_data_nooffset, 'b--', label='m={:.2f}, b={:.2f}'.format(slope_data_nooffset, intercept_data_nooffset))
plt.plot(Q2, slope_simc_nooffset * np.array(Q2) + intercept_simc_nooffset, 'r--', label='m={:.2f}, b={:.2f}'.format(slope_simc_nooffset, intercept_simc_nooffset))
plt.title("No Offset")
plt.xlabel('Q2')
plt.ylabel('Yield')
plt.legend()
'''

# Plot y_data and y_simc on the same plot with linear fits
plt.subplot(211)
plt.plot(Q2, y_data, 'bo', label='y_data')
plt.plot(Q2, y_simc, 'ro', label='y_simc')
plt.plot(Q2, slope_data * np.array(Q2) + intercept_data, 'b--', label='m={:.2f}, b={:.2f}'.format(slope_data, intercept_data))
plt.plot(Q2, slope_simc * np.array(Q2) + intercept_simc, 'r--', label='m={:.2f}, b={:.2f}'.format(slope_simc, intercept_simc))
plt.title("Offset")
plt.xlabel('Q2')
plt.ylabel('Yield')
plt.legend()

'''
# Plot rel_yield on a different plot with a horizontal line at y=1.0
plt.subplot(223)
plt.plot(Q2, rel_yield, 'go')
plt.plot(Q2, slope_relyield_nooffset * np.array(Q2) + intercept_relyield_nooffset, 'g--', label='m={:.2f}, b={:.2f}'.format(slope_relyield_nooffset, intercept_relyield_nooffset))
plt.axhline(y=1.0, color='gray')
plt.xlabel('Q2')
plt.ylabel('Rel. Yield')
plt.legend()
'''

# Duplicate the third plot for better visualization
plt.subplot(212)
plt.plot(Q2, rel_yield, 'go')
plt.plot(Q2, slope_relyield * np.array(Q2) + intercept_relyield, 'g--', label='m={:.2f}, b={:.2f}'.format(slope_relyield, intercept_relyield))
plt.axhline(y=1.0, color='gray')
plt.xlabel('Q2')
plt.ylabel('Rel. Yield')
plt.legend()

# Adjust layout for better spacing
plt.tight_layout()

# Calculate linear fit for y_data_nooffset
#slope_data_nooffset, intercept_data_nooffset, _, _, _ = linregress(Ebeam, y_data_nooffset)

# Calculate linear fit for y_simc_nooffset
#slope_simc_nooffset, intercept_simc_nooffset, _, _, _ = linregress(Ebeam, y_simc_nooffset)

# Calculate linear fit for y_data
slope_data, intercept_data, _, _, _ = linregress(Ebeam, y_data)

# Calculate linear fit for y_simc
slope_simc, intercept_simc, _, _, _ = linregress(Ebeam, y_simc)

#slope_relyield_nooffset, intercept_relyield_nooffset, _, _, _ = linregress(Ebeam, rel_yield_nooffset)

slope_relyield, intercept_relyield, _, _, _ = linregress(Ebeam, rel_yield)

# Plotting
plt.figure(figsize=(12,8))

'''
# Plot y_data_nooffset and y_simc_nooffset on the same plot with linear fits
plt.subplot(221)
plt.plot(Ebeam, y_data_nooffset, 'bo', label='y_data_nooffset')
plt.plot(Ebeam, y_simc_nooffset, 'ro', label='y_simc_nooffset')
plt.plot(Ebeam, slope_data_nooffset * np.array(Ebeam) + intercept_data_nooffset, 'b--', label='m={:.2f}, b={:.2f}'.format(slope_data_nooffset, intercept_data_nooffset))
plt.plot(Ebeam, slope_simc_nooffset * np.array(Ebeam) + intercept_simc_nooffset, 'r--', label='m={:.2f}, b={:.2f}'.format(slope_simc_nooffset, intercept_simc_nooffset))
plt.title("No Offset")
plt.xlabel('Ebeam')
plt.ylabel('Yield')
plt.legend()
'''

# Plot y_data and y_simc on the same plot with linear fits
plt.subplot(221)
plt.plot(Ebeam, y_data, 'bo', label='y_data')
plt.plot(Ebeam, y_simc, 'ro', label='y_simc')
plt.plot(Ebeam, slope_data * np.array(Ebeam) + intercept_data, 'b--', label='m={:.2f}, b={:.2f}'.format(slope_data, intercept_data))
plt.plot(Ebeam, slope_simc * np.array(Ebeam) + intercept_simc, 'r--', label='m={:.2f}, b={:.2f}'.format(slope_simc, intercept_simc))
plt.title("Offset")
plt.xlabel('Ebeam')
plt.ylabel('Yield')
plt.legend()

'''
# Plot rel_yield on a different plot with a horizontal line at y=1.0
plt.subplot(223)
plt.plot(Ebeam, rel_yield, 'go')
plt.plot(Ebeam, slope_relyield_nooffset * np.array(Ebeam) + intercept_relyield_nooffset, 'g--', label='m={:.2f}, b={:.2f}'.format(slope_relyield_nooffset, intercept_relyield_nooffset))
plt.axhline(y=1.0, color='gray')
plt.xlabel('Ebeam')
plt.ylabel('Rel. Yield')
plt.legend()
'''

# Duplicate the third plot for better visualization
plt.subplot(222)
plt.plot(Ebeam, rel_yield, 'go')
plt.plot(Ebeam, slope_relyield * np.array(Ebeam) + intercept_relyield, 'g--', label='m={:.2f}, b={:.2f}'.format(slope_relyield, intercept_relyield))
plt.axhline(y=1.0, color='gray')
plt.xlabel('Ebeam')
plt.ylabel('Rel. Yield')
plt.legend()

# Duplicate the third plot for better visualization
plt.subplot(223)
plt.plot(Q2, rel_yield, 'go')
plt.plot(Q2, slope_relyield * np.array(Q2) + intercept_relyield, 'g--', label='m={:.2f}, b={:.2f}'.format(slope_relyield, intercept_relyield))
plt.axhline(y=1.0, color='gray')
plt.xlabel('Q2')
plt.ylabel('Rel. Yield')
plt.legend()

# Adjust layout for better spacing
plt.tight_layout()

# Show plots
plt.show()
