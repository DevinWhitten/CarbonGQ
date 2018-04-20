#Author: Devin Whitten

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from norm_functions import Spectra, find_cont, Segment
import pandas as pd
import os
plt.ion()


#-------------------------------------------------------------------------------
def normalize(spec_name, ARG):
    spectra = pyfits.open("/Users/MasterD/Google Drive/CarbonDetection/Contents/SEGUE/"+spec_name)
    print "---------- Spectra: ----------"
    print spec_name

    ################################################################################
    # Preparing Spectrum Components
    Wl = np.array(10**spectra[1].data['loglam'], dtype=np.float64)
    Flux = np.array(spectra[1].data['flux'], dtype=np.float64)
    Ivar = np.array(spectra[1].data['ivar'], dtype=np.float64)
    ################################################################################


    #DATA = pd.DataFrame({"Wave":Wl, "Flux":Flux})
    SN_est = np.median(Flux * np.sqrt(Ivar)) # Global Estimate
    print "S/N Estimation:", SN_est

    ###### This is where we make use of norm_functions.py

    SPECTRUM = Spectra(wave=Wl, flux=Flux)
    SPECTRUM.poly_normalize(nlow=1.0, nhigh=5.0, order=5)
    SPECTRUM.spline_normalize(boost=0.03)

    #-------------------------------------------------------------

    print "Flux Average:  ",np.average(Flux)

    print "Continuum Average:  ",np.average(SPECTRUM.continuum)

    ################################################################################
    #       Plotting
    ################################################################################
    fig = plt.figure(figsize=(10,6))

    ax1 = plt.axes([0.1, 0.5,0.55, 0.40], facecolor=[0.93,0.93,0.93])
    ax2 = plt.axes([0.1, 0.10,0.55, 0.40], facecolor=[0.93,0.93,0.93],sharex=ax1)

    ax4 = plt.axes([0.65, 0.5, 0.3, 0.4], sharey=ax1, facecolor=[0.93,0.93,0.93])
    ax3 = plt.axes([0.65, 0.1, 0.3, 0.4], sharey=ax2, facecolor=[0.93,0.93,0.93], sharex=ax4)

    LABEL =""
    HANDLES = [ax1, ax2, ax3, ax4]

    [labels.set_axisbelow(True) for labels in HANDLES]

    ax1.set_title(spec_name)
    ax1.plot(SPECTRUM.wave, SPECTRUM.flux, color="green", linestyle="-",
             alpha=0.5, zorder=3)
    ax1.plot(SPECTRUM.wave, SPECTRUM.Poly_Continuum, color="purple", linestyle="--", linewidth=2.0, zorder=3)
    ax1.plot(SPECTRUM.wave, SPECTRUM.Spline_Continuum, color="orange", linestyle="--", linewidth=1.5, zorder=3, label="Spline")
    ax1.scatter(SPECTRUM.wl_Bin, SPECTRUM.flux_Bin, color="red", marker="s")

    # These are the continuum points
    #ax1.scatter(SPECTRUM.Poly_WL, SPECTRUM.Poly_FLUX, color="blue", marker="^", s=0.5)

    ax2.plot(SPECTRUM.wave, SPECTRUM.norm_flux, color="black", linewidth=0.50, zorder=3, alpha=0.75)

    # -------------------------       SIDE VIEW
    ax4.plot(SPECTRUM.wave, SPECTRUM.flux, color="green", linestyle="-", alpha=0.5, zorder=3)
    #ax4.scatter(SPECTRUM.Poly_WL, SPECTRUM.Poly_FLUX, color="blue", marker="^", s=0.5)
    ax4.plot(SPECTRUM.wave, SPECTRUM.Poly_Continuum, color="purple", linestyle="--", linewidth=1.0, zorder=3, label="Poly")
    ax4.plot(SPECTRUM.wave, SPECTRUM.Spline_Continuum, color="orange", linestyle="--", linewidth=1.5, zorder=3, label="Spline")
    ax4.axvline(x=4341, ymin=0, ymax = max(SPECTRUM.flux), linewidth=0.5, color='k', linestyle="--", label="$H \gamma$")
    ax3.axvline(x=4341, ymin=0, ymax = max(SPECTRUM.flux), linewidth=0.5, color='k', linestyle="--")

    ax3.plot(SPECTRUM.wave, SPECTRUM.norm_flux, color="black", zorder=3, alpha=0.75)


    ax4.legend()

    ax1.set_ylabel("Flux", fontsize=18, fontname="Times New Roman")
    ax2.set_xlabel("Wavelength $\AA$", fontsize=18, fontname="Times New Roman")

    ax4.set_xlim([4150,4650])  #GP Placco et al 2010.
    ax1.set_ylim([0,max(SPECTRUM.flux)])
    ax2.set_xlim([3500,9500])

    ax1.legend(framealpha=0.75, fontsize=14)

    ax1.tick_params(axis='both', which='major', labelsize=15)
    ax2.tick_params(axis='both', which='major', labelsize=15)


    [label.spines['top'].set_visible(False) for label in [ax1, ax4]]
    [label.spines['bottom'].set_visible(False) for label in [ax2,ax3]]

    [label.spines['left'].set_visible(False) for label in [ax1,ax2]]
    [label.spines['right'].set_visible(False) for label in HANDLES]


    [plt.setp(label.get_xticklabels(), visible=False) for label in [ax1, ax4]]
    [plt.setp(label.get_yticklabels(), visible=False) for label in [ax3, ax4]]

    [label.grid(b=True, which='both', linewidth=2,color='1',linestyle='-', zorder=0) for label in HANDLES]
    plt.show()
    if "p" in ARG:
        raw_input("Press any key to continue ... ")
    plt.close()
    # Write Normalized spectra to FITS file

    ################################################################################
    #       Build normalized spectrum .fits file
    ################################################################################

    print "--------------------------------------------"
    print "... Building Normalized .fits File"

    spectra[0].header['SN_est'] = SN_est
    orig_header = spectra[0].header
    new_header = pyfits.PrimaryHDU(header=orig_header)

    orig_data = spectra[1].data
    orig_cols = orig_data.columns
    cols = []
    cols.append(pyfits.Column(name="norm_flux", format="D", array=SPECTRUM.norm_flux))
    new_cols = pyfits.ColDefs(cols)
    hdu = pyfits.BinTableHDU.from_columns(orig_cols + new_cols)

    thdulist = pyfits.HDUList([new_header, hdu])

    ################################################################################
    # Write to norm_<spec_name>.fits
    ################################################################################

    try:
        thdulist.writeto("/Users/MasterD/Google Drive/CarbonDetection/Contents/Normalized/norm_"+spec_name)
    except:
        #Overwrite = raw_input("norm_"+spec_name + " exists. Overwrite?   (y/n):   ")
        #if Overwrite == "y":
        if True:
            os.remove("/Users/MasterD/Google Drive/CarbonDetection/Contents/Normalized/norm_"+spec_name)
            thdulist.writeto("/Users/MasterD/Google Drive/CarbonDetection/Contents/Normalized/norm_"+spec_name)
            print "... Complete: .fits written"
        else:
            print "... Complete: File not written."

    return
