# imports
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from astropy.table import Table
from astropy.io import fits
from PIL import Image
import photutils
from image_functions import *
from urllib.parse import urlencode
from urllib.request import urlretrieve

# setting home directory & appending github
homedir = os.getenv('HOME')
sys.path.append(homedir+'/github/APPSS/python/')

cwd = homedir+'/Desktop/gbtresearch/'
os.chdir(cwd)

# set scales
UNWISE_PIXSCALE = 2.75
LEGACY_PIXSCALE = 1


########################################
### READ IN & CONVERT DATA TO TABLE  ###
########################################   

def gettable(path):
    '''
    Input the path for required data (without CWD path).
    Reads in .fits data.
    Converts .fits to a table and returns the table.
    '''

    fitsdata=path
    datatable=Table.read(fitsdata)
    
    return datatable


########################################
### GETTING RA, DEC, NAMES FROM DATA ###
########################################

def getgaldata(mydata):
    '''
    Takes dataset as input.
    Creates lists for RA, DEC, and galaxy name.
    Appends RA, DEC, and galaxy name.
    Returns lists.
    '''
    
    ra=[]
    dec=[]
    galname=[]
    
    for i in range(len(mydata)):
        ra.append(mydata['RA'][i])
        dec.append(mydata['DEC'][i])
        galname.append(f"AGC{mydata['name'][i]:06d}")
    
    return ra,dec,galname


#############################
### GETTING LEGACY IMAGES ###
#############################

def legacyimage(ra,dec,galname,imsize=120,pixscale=LEGACY_PIXSCALE):
    '''
    Set input RA, DEC, Galaxy Name, with optional imsize and pixscale inputs.
    Bands set in the g,r,z.
    Function will generate .fits files for each band, as well as legacy image.
    Appends these files and images to respective master lists.
    '''
    
    legimfiles = []
    legjpgfiles = []
    bands = ['g','r','z']
    
    imsize_pixels_legacy = round(imsize/LEGACY_PIXSCALE)
    
    try:
        
        for i,b in enumerate(bands):
            url='http://legacysurvey.org/viewer/jpeg-cutout?ra='+str(ra)+'&dec='+str(dec)+'&layer=dr8&size='+str(imsize)+'&pixscale='+str(pixscale)
            urlretrieve(url,'test.jpg')
    
            if i == 0:
                # only need to download this once
                getjpg = True
            else:
                getjpg = False
                
            t = get_legacy_images(ra,dec,galid=galname,band=b,makeplots=False,imsize=str(imsize_pixels_legacy))
            
            if i == 0:
                legimfiles.append(t[0])
                legjpgfiles.append(t[1])
            else:
                legimfiles.append(t[0])
                
    except:
        
        print(f'Warning: {galname} outside Legacy Survey!')
    
    return [legimfiles],[legjpgfiles]


###########################
### GETTING WISE IMAGES ###
###########################

def wiseimage(ra,dec,galname,imsize=120,pixscale=UNWISE_PIXSCALE):
    '''
    Set input RA, DEC, Galaxy Name, with optional imsize and pixscale inputs.
    Function will generate wise image files and noise files for input galaxies.
    Appends these files and images to respective master lists.
    '''
    
    wiseimagefiles = []
    noisefiles = []
    
    try:
        
        url='http://legacysurvey.org/viewer/jpeg-cutout?ra='+str(ra)+'&dec='+str(dec)+'&layer=dr8&size='+str(imsize)+'&pixscale='+str(LEGACY_PIXSCALE)
        urlretrieve(url,'test.jpg')
        
        imsize_pixels_unwise = round(imsize/UNWISE_PIXSCALE)

        t = get_unwise_image(ra,dec,galid=galname,makeplots=False,imsize=str(imsize_pixels_unwise))
    
        wiseimagefiles = t[0]
        noisefiles = t[1]
        wiseimagefiles.sort()
        noisefiles.sort()
    
    except:
        
        print(f'Warning: {galname} outside Legacy Survey!')
    
    return [wiseimagefiles],[noisefiles]


#######################
### PLOTTING IMAGES ###
#######################

def plotting(jpgfile,imfile,wise,galname):
    
    plt.figure(figsize=(12,6.5))

    # concatinate lists
    legacy_images = jpgfile+imfile
    imnames = [f"{galname} : grz",'g','r','z']
    
    # plot legacy images in top row
    for i,im in enumerate(legacy_images):
        
        plt.subplot(2,4,i+1)
        
        if i == 0:
            # display jpg
            t = Image.open(im)
            plt.imshow(t,origin='upper')
            
        else:
            data = fits.getdata(im)
            display_image(data,lowrange=False,percent=95)
            
        plt.title(imnames[i],fontsize=14)

    # plot WISE images
    imnames = ['W1','W2','W3','W4']
    
    for i,im in enumerate(wise):
        plt.subplot(2,4,4+i+1)
        data = fits.getdata(im)
        display_image(data,percent=92)
        plt.title(imnames[i],fontsize=14)
    
    
###########################
### PUT IT ALL TOGETHER ###
###########################    

def getall(ra,dec,galname):
    '''
    Takes in RA, DEC, and galaxy names as input arguments.
    Creates new lists for all legacy, IM, and WISE images/files.
    Uses legacyimage and wiseimage functions to get these files and append them to their lists.
    Returns the four lists of files.
    '''
    
    allimfiles=[]
    alljpgfiles=[]

    allwise=[]
    allnoise=[]
    
    for i in range(len(ra)):
        fitslist,legjpgfiles = legacyimage(ra[i],dec[i],galname[i])
    
        allimfiles += fitslist
        alljpgfiles += legjpgfiles 
    
        wiseimagefiles,noisefiles = wiseimage(ra[i],dec[i],galname[i])
        allwise += wiseimagefiles
        allnoise += noisefiles
    
        # plotting isn't totally needed, made the function to test to see if wise images were
        # being created
        # for greater inputs, comment out plotting line
        # plotting(alljpgfiles[i],allimfiles[i],allwise[i],galname[i])
        
    return allimfiles,alljpgfiles,allwise,allnoise


###############################
### PHOTOMETRY FOR GALAXIES ###
###############################

#%run ~/github/APPSS/python/image_functions.py
#sys.path.append(homedir+'/github/halphagui')

#import photutils

#from photutils.isophote import EllipseGeometry, Ellipse
#from astropy.visualization import ZScaleInterval
#from astropy.visualization import AsymmetricPercentileInterval

#from astropy.table import Table
#from PIL import Image

#import photwrapper

#%load_ext autoreload
#%autoreload 2


#def photometry():
    