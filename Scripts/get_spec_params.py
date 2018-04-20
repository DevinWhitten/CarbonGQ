#AUTHOR: Devin Whitten
#DATE: Nov 11, 2016
# We need to determine the stellar params and signal-to-noise for our spectra
# And produce an normalized output version. We'll tackle the analysis later.
# Just open spectra, cross match SPSPEC with segue for params and SN
# then normalize and write out. Easy


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pyfits

#------------------------------------------------------------------
#open spectra

print "get_spec_params.py"
segue = pd.read_csv("/Users/MasterD/Google Drive/Databases/SEGUE_working.csv")
def get_spec_params(spec_name):
    #spec_name = "spec-0269-51910-0353.fits"
    #spec_name = "spec-2187-54270-0564.fits"
    #spec_name = "spec-0323-51615-0110.csv"
    #spec_name = "spec-0323-51615-0110.fits"
    spectra = pyfits.open("/Users/MasterD/Google Drive/CarbonDetection/Contents/SEGUE/"+spec_name)



    plateid = str(spectra[0].header['PLATEID'])
    mjd = str(spectra[0].header['MJD'])
    fiberid = str(spectra[0].header['FIBERID'])

    while True:
        if len(plateid) < 4:            # Just formatting
            plateid = "0"+plateid
            print "plate Formatted"
        elif len(fiberid) < 4:
            fiberid = "0"+fiberid
            print "fiber Formatted"
        else:
            break

    print "Plate ID:   ",plateid
    print "MJD     :   ",mjd
    print "Fiber ID:   ",fiberid

    spspec = plateid + "-" + mjd + "-" +fiberid

    STAR = segue[segue.SPSPEC == spspec]

    TEFF = float(STAR.TEFF_ADOP)
    LOGG = float(STAR.LOGG_ADOP)
    FEH =  float(STAR.FEH_BIW)
    CFE_COR =  float(STAR.CFE_COR)
    CFE_FEH_ERR = float(STAR.CFE_FEH_ERR)
    #-------------------------------------------
    print
    print "SSPP Params:   "
    print "Teff:    ", TEFF
    print "Log(g):  ", LOGG
    print "[Fe/H]:  ", FEH
    print "[C/Fe] Original:   ", float(STAR.CFE)
    print "[C/Fe] Biweight:  ", CFE_COR


    CFE_ARRAY = np.linspace(-1.50, 3.00, 20)
    print
    print "---- Generating carbon_array.csv ----"
    file_out = open("/Users/MasterD/Google Drive/CarbonDetection/Contents/carbon_array.csv", "w+")
    for CFE in CFE_ARRAY:
        file_out.write(str(TEFF) + "," + str(LOGG) + "," +
                        str(FEH) + "," + str('%7.5f' %CFE) + "\n")
    file_out.close()
    return TEFF, LOGG, FEH, CFE_COR, CFE_FEH_ERR

def generate_carbon_array(LIST, number=20):
    '''
    Precondition: Takes in a pandas dataframe of stellar parameters for spectra generation.
    Postcondition: Loads the appropriate files. carbon_array.csv
    '''
    file_out = open("/Users/MasterD/Google Drive/CarbonDetection/Contents/carbon_array.csv", "w+")

    for STAR in LIST:
        CFE_ARRAY = np.linspace(-1.5, 3.0, number)
        for CFE in CFE_ARRAY:
            file_out.write(str(STAR[0]) + ',' + str(STAR[1]) + "," + str(STAR[2]) + "," + str('%7.5f' %CFE) + "\n")

    file_out.close()
    #raw_input("... New Carbon_array.csv")
    return
