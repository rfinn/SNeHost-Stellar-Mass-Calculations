# WISE-Stellar-Mass-Calculations
Python pipeline for creating stellar mass calculations for galaxies observed in the infrared via the Wide-field Infrared Survey Explorer. Multiple outside repositories are referenced in the tutorial, plus instructions for cloning and working with the repositories are detailed in the tutorial.

**Overview**:
- First, designate a folder/directory that you want to work out of, this directory will be where all of the files and images that are generated will go. This directory can be anywhere on your device, just make sure you can get the path for it (and remember where it is)!
- Second, go to [this link](https://github.com/rfinn/halphagui) and clone the repository to your device. You will then need to go into that cloned repo and create a python virtual environment with Anaconda. Here is a [link](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment) on how to create the virtual environment. One created, activate the environment.
- Next, you will need to create a python environment kernel, install the `requirements.txt` for all of the dependencies of the `halphagui` repo, and install `jupyterlab` via these steps:
  1. `pip install ipykernel`
  2. `pip install -r requirements.txt`
  3. `python -m ipykernel install --user --name=[insert name here]`
  4. `conda install -c conda-forge jupyterlab`
    - If Jupyter Lab does not start, you may have to install `chardet`:
      - `conda install chardet`
- You will also need to install `SWarp` via MacPorts onto *your device*.
  - `SWarp` is used to co-add multiple `.fits` images using any astrometric projection.
  - Go to [this link](https://www.macports.org/install.php) for the MacPorts installation guide. Once MacPorts is successfully installed onto your deivce, in a terminal, install `SWarp` with this line:
    - `port install swarp`
- Once all of these steps are complete, open the `TUTORIAL.ipynb` notebook and select the kernel that you just created.

- Now you can start running the notebook by following the instructions layed out! 
  - Setup the home directory, current working directory, subdirectories for data and plots, and append a path to the `halphagui` repository.
  - Read in your table of galaxy with the *at least* the following data columns:
    - `'RA'`, `'DEC'`, `'AGCnr'`, `'mu'`, `'e_mu'`
    - Note: The table that I provided in this repo does not *have* to be the table that you run, so you can have your own table!
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
