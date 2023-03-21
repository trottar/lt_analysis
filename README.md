![LT_Analyisis Workflow Part 1](docs/LT_Analysis_Workflow_1.png)
![LT_Analyisis Workflow Part 2](docs/LT_Analysis_Workflow_2.png)
![LT_Analyisis Workflow Part 3](docs/LT_Analysis_Workflow_3.png)

- **Note:** I need to update the above diagram as it is no longer correct!

# HeeP Analyisis

- **Note:** I will elaborate on this when I, eventually, update the scripts

1. First run SIMC...

```
./set_HeepInput.sh -a <Kinematics>
```

- **Note:** Flag '-a' is required for analysis, see '-h' for all options

- If you want to check offsets running...

```
./check_Offsets.sh -c <Kinematics>
```

- Again, see '-h' for further explanations

2. Next cuts needs to be applied to data...

```
run_HeeP_Analysis.sh -a <Kinematics>
```

3. Compare data and SIMC...

```
run_HeeP_Analysis.sh <Kinematics>
```

# Production Analyisis

1. First run SIMC...

```
set_ProdInput.sh <epsilon> <Q2> <W>
```

2. Apply PID & CT cuts to the data (see diagram above)...

```
applyCuts_Prod.sh <epsilon> <phi_setting> <Q2> <W> <target> <run_number>
```

3. Combine run numbers for each phi setting, apply diamond cuts, bin data in t & phi, and produce input files for cross-section fortran scripts...

```
run_Prod_Analysis.sh -at <epsilon> <Q2> <W>
```

4. Calculate average kinematic values per t, phi bin and the unseparated cross-section...

```
run_xsect.sh <Q2> <W>
```
