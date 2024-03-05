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

**Folders**:
## dataTables 
#### MissingAGC_fromSNeHostCatalog_20230425_GBTsample.csv
  - Comma separated list of the 16 SNe host galaxies that had missing AGC numbers from the `UATSNeHostGalaxyCatalogue_20230425.csv`.
  - We are working to verify which AGC numbers are the correct numbers based on their SN ID and position. Most of these galaxies have been verified, but there are a few that require the full AGC allsky catalog.
#### SGA_matchedwith_GBT.fits
  - File containing 191 GBT+Archival galaxies that were matched with the Siena Galaxy Atlas. These are the galaxies that I will run bagpipes on to calculate logarithmic stellar masses from photometry.
#### SNe_Hosts_Combined_Sample_20220922.csv
  - The original 334 SNe host galaxies from the GBT+Archival sample.
  - Contains column names: AGCnr, HIMANGA, ALFALFA, SPRINGOB,	radeg,	decdeg,	hiflux,	width,	widtherr,	snr,	D,	vhelagc,	vhelio,	vopt,	v21,	Ty,	a,	b,	b/a,	zmag,	rms,	hisrc,	logMH,	logMH_Limit, and Type.
#### SNe_inUAT_withTypes.fits
  - 
#### UATSNeHostGalaxyCatalogue_20230425.csv
  - 
## galdata
  - Will contain all of the data for each galaxy from the SNe host sample. Each galaxy will have a folder named with its AGC number, and subfolders `figures` containing images/cutouts for the unWISE images and photometry plots, and `imdat` contains all of the .fits and .csv files generated in `StMassCalcs_functions.py`.
## notebooks
  - Nothing in here yet, but I will include sandbox notebooks that I have worked through that created the output stellar mass calculations for the GBT+Archival sample.
## outputs
#### galTable_withMasses.fits
