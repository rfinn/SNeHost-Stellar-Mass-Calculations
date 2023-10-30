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

def getPlots(ptab,galID):
    '''
    Generates plots for the enclosed flux and source magnitude versus semi-major axis from the
    photometry table via photwrapper. These plots are saved under the galPlots subfolder.
    
    Inputs:
        ptab = photwrapper output table for a given galaxy containing data for flux and magnitude for
                for intervals of semi-major axis.
        galID = used only for the title of the plot, has both the AGC prefix and identification number
                 out to 6 digits.
    '''
    
    if os.path.exists(cwd+f'/galPlots/{galID}'):
        pass
    else:
        os.mkdir(cwd+f'/galPlots/{galID}')
    
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
    
    ptabimpath = cwd+f'galPlots/{galID}/'+f'{galID}-ptab-plots.png'
    if os.path.exists(ptabimpath):
        os.remove(ptabimpath)
    plt.savefig(ptabimpath)

    plt.figure(figsize=(12,6.5))
    # plot WISE images
    imname = wiseImgs[0]
    imnames = ['W1','W2','W3','W4']
    for i,im in enumerate(wiseImgs):
        plt.subplot(2,4,4+i+1)
        data = fits.getdata(im)
        display_image(data,percent=92)
        plt.title(imnames[i],fontsize=14)
            
    hdu = fits.open(imname)
    data = hdu[0].data
    hdu.close()
    
    plt.figure()
    plt.title(f'{galID} W1 cut')
    display_image(data)
    
    W1impath = cwd+f'galPlots/{galID}/'+f'{galID}-W1cutout.png'
    if os.path.exists(W1impath):
        os.remove(W1impath)
    plt.savefig(W1impath)
    
#################################################################################

def getPhot(imname):
    '''
    Uses photwrapper, the halphagui repo photometry calculator for a detected galaxy through various 
    photutils methods/functions, to fit concentric ellipses on both the detected galaxy and the background sources.
    The background sources are masked out and the flux is measured at each ring; based on SMA.
    
    Inputs:
        imname = the selected image that you want to measure the flux of. For this these stellar mass calculations,
                     the W1.fits image for each galaxy will be used for the calculations, however you can input other files.
                     Potential future project: running the r band image through the masker and using that mask for the
                                               W1 image for a better mask.
    '''
    
    e = photwrapper.ellipse(imname)
    
    e.detect_objects()
    e.find_central_object()
    e.get_mask_from_segmentation()
    e.get_ellipse_guess()
    e.measure_phot()     
    e.calc_sb()
    e.convert_units()
    e.write_phot_tables()
    e.write_phot_fits_tables()
    
    ptabName = imname.replace('.fits','-phot.fits')
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

def getMass(myGal,plotting=False,imsize=120,makeMask=True,changeSize=False,verbose=False):
    '''
    Reads in a table of galaxy coordinates to get wise/legacy images to calculate photometry; W1 phot will be converted
    into a flux, then converted into a stellar mass.
    
    Inputs:
        myGal = table entry of galactic coordinate data; needs 'RA', 'DEC', 'AGCnr', 'mu', and 'e_mu' columns
        plotting = conditional to plot each photwrapper image step; default set to False
        imsize = length/width of legacy/WISE image cutout; default set to 120; units: pixels
        makeMask = conditional to mask off subtracted sky objects from segmentation image; default set to True
        changeSize = conditional to search for previously made unwise images, deletes them, and remakes them with given imsize; 
                        default set to False
        verbose = conditional to have function talk to you through the processes within; be verbose!
    
    RETURNS:
        ptab = table created from output photometry file via photwrapper
        e = photwrapper ellipse data
        galTable_withMasses = same as the table you read in, but with new columsn for calculated stellar masses and
                                image size. This table updates every time you call the function!
    '''
    
    # get sky coords and galaxy AGC ID from input table
    ind = myGal.index
    ra = myGal['RA']
    dec = myGal['DEC']
    galID = f"AGC{myGal['AGCnr']:06d}"; galNum = myGal['AGCnr']
    
    # setting pixel scaling
    UNWISE_PIXSCALE = 2.75
    LEGACY_PIXSCALE = 1
    
    if changeSize:
        # checks to see if unwise images were created
        dstring=f'{galID}-unwise*'
        flist=glob.glob(dstring)
        # removes image if it exists in the directory
        for f in flist:
            if os.path.exists(f):
                os.remove(f)
    
    # gets the W1-4 and legacy images for the galaxy
    legacyImgs, wiseImgs = display_legacy_unwise(ra,dec,galID,imsize_arcsec=imsize)
    
    # can be changed to different bands
    imname = wiseImgs[0]
    
    %matplotlib inline
    
    if plotting:
    
        print()
    
        plt.figure(figsize=(12,6.5))

        # concatinate lists
        imnames = ['grz','g','r','z']
        # plot legacy images in top row
        for i,im in enumerate(legacyImgs):
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
        for i,im in enumerate(wiseImgs):
            plt.subplot(2,4,4+i+1)
            data = fits.getdata(im)
            display_image(data,percent=92)
            plt.title(imnames[i],fontsize=14)
            
        hdu = fits.open(imname)
        data = hdu[0].data
        hdu.close()
    
        plt.figure()
        plt.title(f'{imname}')
        display_image(data)
    
    ### PHOTWRAPPER ####
        
    ptab,e=getPhot(imname)
        
    if plotting:
        
        e.show_seg_aperture()
        e.draw_guess_ellipse_mpl()
        e.draw_phot_apertures()
    
    # Plot F_enc vs SMA and Mag vs SMA, then save the plot as .png in galPlots folder.
    
    getPlots(ptab,galID)
    plt.show()

    if verbose:
        print(f"Max enclosed flux: {np.max(ptab['flux']):.2e}")
        print(f"Max SMA: {np.max(ptab['sma_arcsec'])}"); print()

    #################################
    ### Stellar Mass Calculations ###
    #################################
    
    maxFenc=np.max(ptab['flux_erg']) # max enclosed flux in ergs/cm^2 s
    mu=myGal['mu'] # distance modulus
    mu_err=myGal['e_mu'] # distance modulus error
    mmag=np.min(ptab['mag']) # source magnitude via photutils

    
    LW1_mag,LW1_flux=getLW1(maxFenc,mu,mu_err,mmag)        

    logMstar_mag=getLogMass(LW1_mag)
    logMstar_flux=getLogMass(LW1_flux)

    if verbose:
        print(f"LW1_mag array: {LW1_mag}"); print()
        print(f"LW1_flux array: {LW1_flux}"); print()
        
        print(f'Stellar masses based on magnitude: (upper, median, lower) \n{logMstar_mag}'); print()
        print(f'Stellar masses based on flux: (upper, median, lower) \n{logMstar_flux}'); print()

    if os.path.isfile(cwd+'data/galTable_withMasses.fits'):
        if verbose:
            print('Mass-updated table found! Reading in table...')
        # reads in the masses table
        galTable_withMasses = Table.read(cwd+'data/galTable_withMasses.fits')

        # appends the calculated masses and image size to the respective rows
        galTable_withMasses[ind]['logMstar_mag']=logMstar_mag.reshape((3,1))
        galTable_withMasses[ind]['logMstar_flux']=logMstar_flux.reshape((3,1))
        galTable_withMasses[ind]['Imsize']=int(imsize)
        
    else:
        if verbose:
            print('Mass-updated table not found, creating table...')
        # creates a duplicate table
        galTable_withMasses = Table.read(destinationPath)
        
        # creates columns for the stellar masses and image sizes
        logMstar_mag_col = Column([np.zeros((3, 1), dtype=float)], name='logMstar_mag')
        logMstar_flux_col = Column([np.zeros((3, 1), dtype=float)], name='logMstar_flux')
        imsize_col = Column(int(0),name='Imsize')
        
        # adds those columns to the table
        galTable_withMasses.add_column(logMstar_mag_col)
        galTable_withMasses.add_column(logMstar_flux_col)
        galTable_withMasses.add_column(imsize_col)
        
        # appends the calculated masses and image size to the respective rows
        galTable_withMasses[ind]['logMstar_mag']=logMstar_mag.reshape((3,1))
        galTable_withMasses[ind]['logMstar_flux']=logMstar_flux.reshape((3,1))
        galTable_withMasses[ind]['Imsize']=int(imsize)

    # this rewrites the galTable_withMasses.fits file in the data folder every time a new mass is calculated,
    # which will keep track of all the masses that have been done without altering the primary fits file!
    galTable_withMasses.write('data/galTable_withMasses.fits',format='fits',overwrite=True)
    
    return ptab,e,galTable_withMasses

#################################################################################

def getMasses(galTab,verbose=False):
    '''
    Reads in a table of galaxy coordinates to get wise/legacy images to calculate photometry; W1 phot will be converted
    into a flux, then converted into a stellar mass.
    
    Inputs:
        galTab = table of galaxies to get stellar mass calculations, includes the columns:
                    'AGRnr', 'RA', 'DEC', 'mu', 'e_mu', and 'Imsize'
        verbose = conditional to have function talk to you through the processes within; be verbose!
    
    RETURNS:
        galTable_withMasses = duplicated galTab, but with new columns 'logMstar_mag' and 'logMstar_flux' 
                                for the respective calculated masses.
    '''
        
    for row in galTab:
        
        # get sky coords and galaxy AGC ID from input table
        ind = row.index
        ra = row['RA']
        dec = row['DEC']
        galID = f"AGC{row['AGCnr']:06d}"; galNum = row['AGCnr']
        imsize = row['Imsize']
    
        # setting pixel scaling
        UNWISE_PIXSCALE = 2.75
        LEGACY_PIXSCALE = 1
        
        # gets the W1-4 and legacy images for the galaxy
        try:
            legacyImgs, wiseImgs = display_legacy_unwise(ra,dec,galID,imsize_arcsec=imsize)
        except:
            continue
    
        # can be changed to different bands
        imname = wiseImgs[0]
    
        ### PHOTWRAPPER ####
        ptab,e=getPhot(imname)

        # Plot F_enc vs SMA and Mag vs SMA, then save the plot as .png in galPlots folder.
        getPlots(ptab,galID,wiseImgs)
        
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

        if os.path.isfile(cwd+'data/galTable_withMasses.fits'):
            # reads in the masses table
            galTable_withMasses = Table.read(cwd+'data/galTable_withMasses.fits')

            # appends the calculated masses and image size to the respective rows
            galTable_withMasses['logMstar_mag'][ind]=logMstar_mag.reshape((3,1))
            galTable_withMasses['logMstar_flux'][ind]=logMstar_flux.reshape((3,1))
        
        else:
            # creates a duplicate table
            galTable_withMasses = Table.read(myGalTab_path)
        
            # creates columns for the stellar masses and image sizes
            logMstar_mag_col = Column([np.zeros((3, 1), dtype=float)], name='logMstar_mag')
            logMstar_flux_col = Column([np.zeros((3, 1), dtype=float)], name='logMstar_flux')
        
            # adds those columns to the table
            galTable_withMasses.add_column(logMstar_mag_col)
            galTable_withMasses.add_column(logMstar_flux_col)
            
            # appends the calculated masses and image size to the respective rows
            galTable_withMasses['logMstar_mag'][ind]=logMstar_mag.reshape((3,1))
            galTable_withMasses['logMstar_flux'][ind]=logMstar_flux.reshape((3,1))

        # this rewrites the galTable_withMasses.fits file in the data folder every time a new mass is calculated,
        # which will keep track of all the masses that have been done without altering the primary fits file!
        galTable_withMasses.write('data/galTable_withMasses.fits',format='fits',overwrite=True)
            
    return galTable_withMasses

#################################################################################