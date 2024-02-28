# WISE Stellar Mass Calculations
Python pipeline for creating stellar mass calculations for galaxies observed in the infrared via the Wide-field Infrared Survey Explorer. Multiple outside repositories are referenced in the tutorial, plus instructions for cloning and working with the repositories are detailed in the tutorial.

**Requirements**:
- First, designate a folder/directory that you want to work out of, this directory will be where all of the files and images that are generated will go. This directory can be anywhere on your device, just make sure you can get the path for it (and remember where it is)!
- Second, you will need to clone two repositories: This respository (GBT-Stellar-Mass-Calculations) and the halphagui repository.
  - Go to [this link](https://github.com/rfinn/halphagui) and clone the repository to your device. Do the same for this repository. Make sure you clone them to somewhere you remember!
- You will then need to go into that cloned repo and create a python virtual environment with Anaconda. Here is a [link](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment) on how to create the virtual environment. One created, activate the environment.
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
- Once all of these steps are complete, launch Jupyter Lab, create a notebook with the `halphagui` kernel, and run the `StMassCalcs_functions.py` file! Steps for this are below.

- For running the `StMassCalcs_functions.py` file:
  1. Make sure the selected kernel for your notebook is the kernel you created earlier!
  2. Use this line of code to run the python file:
    - `%run StMassCalcs_functions.py --cwd [full path to desired working directory] --tablepath [full path to your compiled galaxy table]`
      - Note: You might have to provide the path to `StMassCalcs_functions.py`.
  3. You will then have access to your table under the variable name `mytable` and all of the functions in the python file.
  - Call the `getMasses` function with your table, let it run, and voilà!
    - The output table `galTable_withMasses.fits` automatically is downloaded into your `data` folder in your current working directory.

# UATGBT_dataTab.fits
- This `.fits` file is the table that I used for the calculations, including the data for 191 galaxies in the GBT SNe host galaxies. This is the format that will be needed for your table! You can include other columns for further analysis, but refer to the **Overview** section for the minimum requirements.

# galTable_withMasses.fits
- This is the output table after calling the `getMasses` function. It is a duplicate of the input table (for this specific table, it is a duplicate of the `UATGBT_dataTab.fits`), but now it has new columns for the calculated stellar masses.

# StMassCalcs-functions.py
- (***WIP***) This will be a python file containing all of the functions needed to calculate the stellar masses.
