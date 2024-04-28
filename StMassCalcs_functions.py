import os;import sys

import warnings
warnings.filterwarnings('ignore')

import glob
import argparse

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('Agg')

from astropy.table import Table,Column
from astropy.io import fits

from PIL import Image

parser = argparse.ArgumentParser()
parser.add_argument('--cwd',dest='cwd',default=None,
        help='Set the current working directory for the tutorial.')
parser.add_argument('--tablepath',dest='tablepath',default=None,
        help='Set the path to where your data table is located.')
parser.add_argument('--overwrite',dest='overwrite',default=False,
        help='Overwrites check to see if mass has already been calculated, set to True if you need to recalculate the masses.')
args = parser.parse_args()

homedir = os.getenv('HOME')

if args.cwd is not None:
    cwd=args.cwd
    os.chdir(cwd)

if os.path.exists(cwd+'data/'):
    print('\nData Table (/data) subfolder already exists.\n')
else:
    print('Data Table (/data) subfolder does not exists, creating one now.\n')
    os.mkdir(cwd+'data/')

sys.path.append(homedir+'/github/halphagui/')

import photwrapper
from image_functions import *

if args.tablepath is not None:
    mytable=Table.read(args.tablepath)
    print('Your Data Table was read in as `mytable`.\n')

#################################################################################

def getLW1(Fenc,mu,e_mu,m):
    '''
    Takes in ptab and matched table values for distance modulus plus error, then calculates the luminosity two
    different ways; one in relation to the distance modulus and absolute magnitude, and the other through the flux-
    luminosity relation using the max enclosed flux from photutils.

    Inputs:
        Fenc = max enclosed flux of target galaxy; ergs/cm^2 s
        mu = target galaxy distance modulus
        e_mu = target galaxy error in the distance modulus
        m = source apparent magnitude of target galaxy

    RETURNS:
        LW1_mag = numpy array including the upper, median, and lower W1 Luminosities via distance modulus, mu, 
                     and absolute magnitude of galaxy; erg/s.
        LW1_flux = numpy array including the upper, median, and lower W1 Luminosities via flux-luminosity relation;
                     erg/s.
    '''
    
    cmConversion = 3.086e18
    d_Upper=(10**((mu+e_mu)/5+1))*cmConversion # gets d in pc, converts to cm; includes upper error for mu
    d_Med=(10**(mu/5+1))*cmConversion # gets d in pc, converts to cm
    d_Lower=(10**((mu-e_mu)/5+1))*cmConversion # gets d in pc, converts to cm; includes lower error for mu
    d=np.array([d_Upper,d_Med,d_Lower])

    M_Upper=m-(mu+e_mu)
    M_Med=m-mu
    M_Lower=m-(mu-e_mu)
    
    M=np.array([M_Upper,M_Med,M_Lower])
    
    LW1_mag_Upper=10**(-0.4*(M_Upper-3.24))
    LW1_mag_Med=10**(-0.4*(M_Med-3.24))
    LW1_mag_Lower=10**(-0.4*(M_Lower-3.24))
    
    LW1_mag=np.array([LW1_mag_Upper,LW1_mag_Med,LW1_mag_Lower])

    ### Luminosity via enclosed flux
    Lsun_ergs = 3.846e33
    LW1_flux_Upper=Fenc*(4*np.pi*d_Upper**2)/Lsun_ergs
    LW1_flux_Med=Fenc*(4*np.pi*d_Med**2)/Lsun_ergs
    LW1_flux_Lower=Fenc*(4*np.pi*d_Lower**2)/Lsun_ergs
    
    LW1_flux=np.array([LW1_flux_Upper,LW1_flux_Med,LW1_flux_Lower])

    return LW1_mag,LW1_flux

#################################################################################

def getPlots(ptab,galID,image_names):
    '''
    Generates plots for the enclosed flux and source magnitude versus semi-major axis from the
    photometry table via photwrapper. These plots are saved under their own subfolder.
    
    Inputs:
        ptab = photwrapper output table for a given galaxy containing data for flux and magnitude for
                for intervals of semi-major axis.
        galID = used only for the title of the plot, has both the AGC prefix and identification number
                 out to 6 digits.
        image_names = list of the image file names for WISE W1-4.
    '''
    
    SMA = ptab['sma_arcsec'] # semi-major axis measured in arcseconds
    FLUX = ptab['flux_erg'] # enclosed flux measured in ergs
    MAG = ptab['mag'] # source magnitude through each fitted ellipse
    
    plt.figure(figsize=(12,5))
    
    plt.subplot(1,2,1)
    plt.plot(SMA,FLUX)
    plt.title(f'{galID} $F_{{enc}}$ Vs. SMA',fontsize=16)
    plt.xlabel("Semi-major axis (arcsec)",fontsize=14)
    plt.ylabel("$F_{enc}$ (erg / $cm^2$ s)",fontsize=14)
    
    plt.subplot(1,2,2)
    plt.plot(SMA,MAG)
    plt.title(f'{galID} Magnitude Vs. SMA',fontsize=16)
    plt.xlabel("Semi-major axis (arcsec)",fontsize=14)
    plt.ylabel("Magnitude",fontsize=14)    
    
    ptabimpath = args.cwd+f'{galID}/figures/'+f'{galID}-ptab-plots.png'
    if os.path.exists(ptabimpath):
        os.remove(ptabimpath)
    plt.savefig(ptabimpath)
    plt.close()

    plt.figure(figsize=(12,6.5))
    # plot WISE images
    image_names.sort()
    imname = image_names[0]
    imnames = ['W1','W2','W3','W4']
    for i,im in enumerate(image_names):
        plt.subplot(1,4,i+1)
        data = fits.getdata(im)
        display_image(data,percent=92)
        plt.title(imnames[i],fontsize=14)
            
    hdu = fits.open(imname)
    data = hdu[0].data
    hdu.close()
    
    W1impath = args.cwd+f'{galID}/figures/'+f'{galID}-W1cutout.png'
    if os.path.exists(W1impath):
        os.remove(W1impath)
    plt.savefig(W1impath)
    plt.close()
    
#################################################################################

def getPhot(imname,galID,ngrow=1):
    '''
    Uses photwrapper, the halphagui repo photometry calculator for a detected galaxy through various 
    photutils methods/functions, to fit concentric ellipses on both the detected galaxy and the background sources.
    The background sources are masked out and the flux is measured at each ring; based on SMA.
    
    Inputs:
        imname = the selected image that you want to measure the flux of. For this these stellar mass calculations,
                     the W1.fits image for each galaxy will be used for the calculations, however you can input other files.
        galID = AGC prefix and identification number for target galaxy out to 6 digits.
        ngrow = default value of 1 for WISE images, can be changed depending on what images you are analyzing.

    Returns:
        ptab = photwrapper output table for a given galaxy containing data for flux and magnitude for
                for intervals of semi-major axis. 
        e = photwrapper ellipse fitting/analysis for target galaxy.
    '''
        
    # check if mask exists
    maskname = imname.replace('.fits','-mask.fits')
    maskpath = args.cwd+f'{galID}/imdat/{maskname}'
    
    # looking for mask_wrapper file
    if os.path.exists(maskpath):
        maskfile = maskpath
        print(f'\nMaskfile for {galID} found! No need to generate one.\n')

        e = photwrapper.ellipse(imname,mask=maskpath)

    else:
        print(f'\nNo maskfile for {galID} found! Using maskwrapper.py to generate one now.\n')
        
        # call the mask wrapper from within your script
        os.system(f"python ~/github/halphagui/maskwrapper.py --image {imname} --ngrow {ngrow} --sesnr {2} --minarea {5} --auto")
        maskfile = maskpath
    
        e = photwrapper.ellipse(imname,mask=maskname)
    
    e.detect_objects()
    e.find_central_object()    
    e.get_ellipse_guess()    
    e.measure_phot()     
    e.calc_sb()
    e.convert_units()
    e.write_phot_tables()
    e.write_phot_fits_tables()
    
    allApertures_filename = f'{galID}-allApertures.png'
    e.show_seg_aperture(plotname=allApertures_filename)
    plt.savefig(args.cwd+f'{galID}/figures/{allApertures_filename}')
    plt.close()
    
    photApertures_filename = f'{galID}-photApertures.png'
    e.draw_phot_apertures(plotname=photApertures_filename)
    plt.savefig(args.cwd+f'{galID}/figures/{photApertures_filename}')
    plt.close()
    
    ptabName = imname.replace('.fits','_phot.fits')
    ptab = Table.read(ptabName)
    
    return ptab,e

#################################################################################

def getLogMass(LW1_array):
    '''
    Reads in the W1 luminosities calculations with the +/- error on the distance modulus, including the 
    calculation without error propagation, and converts those luminosities into stellar masses, via 
    equation from Jarrett 2023.

    Inputs:
        LW1_array = contains the upper, median, and lower luminosities as a (3,) shape array:
            LW1_upper = W1 luminosity including the +error on the distance modulus; distance modulus will 
                            be the largest source of error; erg/s
            LW1_med =  W1 luminosity NOT including the +/-error on the distance modulus; erg/s
            LW1_lower = W1 luminosity including the -error on the distance modulus; distance modulus will 
                            be the largest source of error; erg/s
        
    RETURNS:
        logMstar = numpy array of including the upper, median, and lower W1 stellar mass calculations.
    '''
    
    A0= -12.62185; A1= 5.00155; A2= -0.43857; A3= 0.01593
    
    logMstar_upper=A0+(A1*np.log10(LW1_array[0]))+(A2*(np.log10(LW1_array[0]))**2)+(A3*(np.log10(LW1_array[0]))**3)
    logMstar_med=A0+(A1*np.log10(LW1_array[1]))+(A2*(np.log10(LW1_array[1]))**2)+(A3*(np.log10(LW1_array[1]))**3)
    logMstar_lower=A0+(A1*np.log10(LW1_array[2]))+(A2*(np.log10(LW1_array[2]))**2)+(A3*(np.log10(LW1_array[2]))**3)
    
    logMstar=np.array([logMstar_upper,logMstar_med,logMstar_lower])
    
    return logMstar

#################################################################################

def getMasses(galTab,verbose=False):
    '''
    Reads in a table of galaxy coordinates to get wise/legacy images to calculate photometry; W1 phot will be converted
    into a flux, then converted into a stellar mass.
    
    Inputs:
        galTab = table that has galaxies that need to have stellar masses calculated, was read in at the beginning of the notebook.
                    Must have these columns: 'RA' , 'DEC' , 'AGCnr' , 'Imsize' , 'mu' , and 'e_mu' !!!
        verbose = conditional to have function talk to you through the processes within; be verbose!
    
    RETURNS:
        galTable_withMasses = duplicate table of galTab but includes new columns for the calculated stellar masses
                                for magnitude and flux, plus their upper and lower calculations, too.
    '''

    massPath = args.cwd+'data/galTable_withMasses.fits'
    if os.path.isfile(massPath):
        # reads in the masses table
        massTab = Table.read(massPath)
        if verbose:
            print('Table with calculated masses was found and read in.\n')
    else:
        # Create a duplicate table of input table.
        massTab = Table.read(args.tablepath)          

        # Create and add columns with 0.0 entries for the log stellar masses and their upper and lower bounds.
        flux_col_Median = Column(0.0, name='logMstar_flux_Median'); massTab.add_column(flux_col_Median)
        flux_col_Upper = Column(0.0, name='logMstar_flux_Upper'); massTab.add_column(flux_col_Upper)
        flux_col_Lower = Column(0.0, name='logMstar_flux_Lower'); massTab.add_column(flux_col_Lower)
        
    RaCol=input('Enter RA column header:'); 
    DecCol=input('Enter DEC column header:'); 
    GalIDCol=input('Enter AGC number column header:'); print('\n') 

    for row in galTab:
        
        # get sky coords and galaxy AGC ID from input table
        ind = row.index
        ra = float(row[RaCol])
        dec = float(row[DecCol])
        galID = f"AGC{row[GalIDCol]:06d}"; galNum = row[GalIDCol]
        imgsize = 120

        # checks to see if overwrite is set to True.
        # if False: check to see if mass has been calculated and to continue to next galaxy if there is one.
        
        if not args.overwrite:
            if massTab['logMstar_flux_Median'][ind]>0.0:
                if verbose:
                    print(f'{galID} mass already calculated, moving to next galaxy.')
                continue

        print(f'\n------BEGIN CALCULATION: {galID}------\n')
    
        # Checks to see if galaxy subfolder exists, along with its own image data and figure subfolders.
        if not os.path.exists(args.cwd+galID):
            os.mkdir(args.cwd+galID)
            os.mkdir(args.cwd+f'{galID}/imdat/')
            os.mkdir(args.cwd+f'{galID}/figures/')
            if verbose:
                print(f'Galaxy {galID} subfolders were not found, created them now.')
        if verbose:
            print(f'Galaxy {galID} subfolders already created.\n')
        
        # setting pixel scaling
        UNWISE_PIXSCALE = 2.75
        LEGACY_PIXSCALE = 1

        ###########################
        ### Getting WISE Images ###
        ###########################
        
        try:
            image_names,weight_names,multiframe = get_unwise_image(ra, dec, galid=galID,
                                                    pixscale=UNWISE_PIXSCALE, imsize=imgsize, bands='1234',
                                                    makeplots=False,subfolder=None,verbose=False)
        except:
            print(f'\nWarning: Problem with getting UNWISE images for {galID}!\n')
            continue

        # We need to sort the images because sometimes the images will be out of order in the array.
        # This ensures that we always grab the W1 image for the calculations, but it can be changed if needed.

        image_names.sort()
        imname = image_names[0]
        if verbose:
            print(f'\nunWISE images: {image_names}')

        ##########################################
        ### Maskwrapper & Photwrapper Analysis ###
        ##########################################
        
        try:
            ptab,e=getPhot(imname,galID,ngrow=1)
        except:
            print(f'\nWarning: Problem with photometry for {galID}!\n')
            continue

        # Plot F_enc vs SMA and Mag vs SMA, then save the plot as .png in the galaxy folder.
        
        getPlots(ptab,galID,image_names)
        
        #################################
        ### Stellar Mass Calculations ###
        #################################
    
        maxFenc=np.max(ptab['flux_erg']) # max enclosed flux in ergs/cm^2 s
        mu=row['mu'] # distance modulus
        mu_err=row['e_mu'] # distance modulus error
        mmag=np.min(ptab['mag']) # source magnitude via photutils

        LW1_mag,LW1_flux=getLW1(maxFenc,mu,mu_err,mmag)

        logMstar_mag=getLogMass(LW1_mag).round(5)
        logMstar_flux=getLogMass(LW1_flux).round(5)
        if verbose:
            print(f'\nLW1_mag = {LW1_mag}'+f'\nLW1_flux = {LW1_flux}\n')
            
        massTab['logMstar_flux_Median'][ind]=logMstar_flux[1]
        massTab['logMstar_flux_Upper'][ind]=logMstar_flux[0]
        massTab['logMstar_flux_Lower'][ind]=logMstar_flux[2]
        
        if verbose:
            print(f'logMstar_flux_Median = {logMstar_flux[1]}')
            print(f'\nlogMstar_flux_Upper = {logMstar_flux[0]}')
            print(f'\nlogMstar_flux_Lower = {logMstar_flux[2]}')

        # this rewrites the galTable_withMasses.fits file in the data folder every time a new mass is calculated,
        # which will keep track of all the masses that have been done without altering the primary fits file!
        massTab.write(massPath,format='fits',overwrite=True)

        # gathering all of the plots made and removing/renaming them to fall
        # into respective galaxy image folder
        galplots=glob.glob(galID+'*.png')
        for img in galplots:
            os.remove(img)
        
        imdata=glob.glob(galID+'-*')
        for img in imdata:
            if not os.path.exists(args.cwd+f'{galID}/imdat/{img}'):
                os.rename(img,args.cwd+f'{galID}/imdat/{img}')
            else:
                os.remove(img)
        if verbose:
            print(f'Removed {len(galplots)} plots for {galID}.\n'+f'Renamed/Removed {len(imdata)} data files for {galID}.')
        print(f'\n------END CALCULATION: {galID}------\n')
            
    print('\nCalculations complete, view galTable_withMasses.fits for masses!')
    print('\nReading in galTable_withMasses.fits now...')
    
    return massTab
    