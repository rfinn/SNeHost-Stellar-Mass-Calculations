# GBT-Stellar-Mass-Calculations
Python pipeline for creating stellar mass calculations for galaxies observed in the infrared via the Wide-field Infrared Survey Explorer. Multiple outside repositories are referenced in the tutorial, plus instructions for cloning and working with the repositories are detailed in the tutorial.

**Overview**:
- Setup a conda environment and launch the `TUTORIAL.ipynb` notebook with the environment kernel.
- Steup home directory, current working directory, subdirectories for data and plots, and append a path to the `halphagui` repository.
- Read in your table of galaxy with the *at least* the following data columns:
  - `'RA'`, `'DEC'`, `'AGCnr'`, `'mu'`, `'e_mu'`
- Run the cells containing all of the functions
- Call the `getMasses` function with your table, let it run, and voilà!
  - The output table `galTable_withMasses.fits` automatically is downloaded into your `data` folder in your current working directory.

# TUTORIAL.ipynb
- Main python notebook where the calculations are done. All of the code is set, just download a copy, and follow the directions & change the what is necessary!

# UATGBT_dataTab.fits
- This `.fits` file is the table that I used for the calculations, including the data for 191 galaxies in the GBT SNe host galaxies. This is the format that will be needed for your table! You can include other columns for further analysis, but refer to the **Overview** section for the minimum requirements.

# galTable_withMasses.fits
- This is the output table after calling the `getMasses` function. It is a duplicate of the input table (for this specific table, it is a duplicate of the `UATGBT_dataTab.fits`), but now it has new columns for the calculated stellar masses.

# StMassCalcs-functions.py
- (***WIP***) This will be a python file containing all of the functions needed to calculate the stellar masses.
