# Author: Devin Whitten
# Date: Nov 14th 2016
# Here we will actually do the integration
# Assume spectra is already normalized, and interpolation has been run

import pyfits
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd
import os
import scipy.interpolate as interp
from scipy.interpolate import interp1d
from scipy.stats import linregress
from IntegrateFunctions import Spectrum, G_index, LINE, Load_Synthetic

plt.ion()

print "--------------- integrate.py ---------------"
#-------------------------------------------------------------------------------

def integrate(spec_name, PARAMS, iterations=150, ARG = ""):
    spec_name = "norm_" + spec_name
    #-------------------------------------------------------------------------------
    print spec_name
    print "---------------------------------------------------"

    ################################################################################
    # Preparing Spectrum Components
    norm_spec = pyfits.open("/Users/MasterD/Google Drive/CarbonDetection/Contents/Normalized/" + spec_name)
    wavelength = 10**norm_spec[1].data['loglam']
    flux = norm_spec[1].data['flux']
    norm_flux = norm_spec[1].data['norm_flux']
    ivar = norm_spec[1].data['ivar']

    ################################################################################
    #                               For testing only!
    #test = pd.read_csv("Spectra/T6629.0000g4.1010z-2.7820c+0.9375.csv")
    #wavelength=test.wl
    #norm=test.wl
    #Segue = Spectrum(name="Whatever", wl= wavelength, norm_flux = norm, flux=norm, ivar=norm)
    ################################################################################

    Segue = Spectrum(name=spec_name, wl = wavelength, norm_flux = norm_flux, ivar = ivar, flux = flux)

    Segue.compute_sn()
    print "Signal-to-noise:  ", Segue.SN, "+/-", Segue.SN_sigma
    #raw_input()
    SN_430 = np.random.normal(Segue.SN, Segue.SN_sigma,iterations)
    #SN_430 = np.random.normal(100.0, 1.0, 100)

    ################################################################################
    #                       Load in Synthetic Spectra
    ################################################################################

    print "Loading synthetic spectra ..."

    ZERO, BATCH = Load_Synthetic(PARAMS) # Also loads parameters like [C/Fe]



    MY_INDEX = Segue.compute_SEGUE_GQ(ZERO= ZERO)
    print "My Index:  ", MY_INDEX
    #raw_input()
    # Interpolate Zero carbon synth spectrum for subtraction

    #MY_INDEX = np.sum(Segue.norm_flux[(Segue.wave > GMIN)&(Segue.wave<GMAX)] - Interp_0c(Segue.wave)[(Segue.wave > GMIN)&(Segue.wave<GMAX)])
    #----------------------------------------------------------------------

    #GMIN, GMAX = 4297.5, 4312.5 # GPa Placco et al. 2010
    GMIN, GMAX = 4281.0, 4307.0 # GPHES Placco et al. 2010

    CFE_Pred_INT = []
    CFE_Pred_POLY = []
    CFE_Pred_LINE = []

    fig = plt.figure(figsize=(10,5))
    ax1 = plt.axes([0.1,0.1,0.4, 0.8], facecolor=[0.93,0.93,0.93])
    ax2 = plt.axes([0.5, 0.1, 0.4, 0.8], facecolor=[0.93,0.93,0.93], sharey=ax1)
    #ax1.set_xscale("log")
    [labels.set_axisbelow(True) for labels in (ax1,ax2)]

    [label.spines['top'].set_visible(False) for label in [ax1, ax2]]
    [label.spines['bottom'].set_visible(False) for label in [ax1,ax2]]
    [label.grid(b=True, which='both', linewidth=2,color='1',linestyle='-', zorder=0) for label in (ax1, ax2)]
    ax1.spines['left'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    plt.setp(ax2.get_yticklabels(), visible=False)
    #ax1.set_xlim([-30,0])
    #Synth_Indices = ax.scatter([0]*24,[0]*24,marker="^", color="green", s=20)
    #Interp = ax.plot([0]*50,[0]*50,linewidth=2.0, color="red", linestyle="--")
    #Poly = ax.plot([0]*50,[0]*50,linewidth=2.0, color="blue", linestyle="--", label="")
    #Line_Regress = ax.plot([0]*50,[0]*50,linewidth=2.0, color="green", linestyle="--", label="")

    Synth_Indices = Line2D([], [],marker="^", color="green")
    Interp = Line2D([], [],linewidth=1.0, color="red", linestyle="--", label="Interpolation")
    Poly = Line2D([], [], linewidth=1.0, color="blue", linestyle="--", label="Poly Regression")
    Line_Regress = Line2D([], [], linewidth=1.0, color="black", linestyle="--", label="Linear Regression")
    ax1.legend()

    ax2.set_xlabel("Iteration", fontsize=14)
    ax1.tick_params("both", labelsize =14)
    ax1.set_xlabel("GQ index", fontsize=14)
    ax1.set_ylabel("[C/Fe]", fontsize=14)
    ax1.set_ylim([-2,3])
    ax1.set_xlim([0,30])

    ax1.add_line(Synth_Indices)
    ax1.add_line(Interp)
    ax1.add_line(Poly)
    ax1.add_line(Line_Regress)

    POLY_AVG, INTERP_AVG, LINE_AVG = [],[],[]
    #Synth_Indices.set_data()
    ################################################################################
    #                       Main Subtraction Routine
    ################################################################################

    for SN in SN_430:
        # Here we inject noise in all spectra from the Gaussian SN_430,
        # but should keep the zero carbon spectrum noise free

        #print "Flux Average Before:  ", np.average(BATCH[1].noise_flux), np.average(BATCH[1].norm_flux)

        [element.sn_inject(SN) for element in BATCH]

        #print "Following Noise Injection:  ",  np.average(BATCH[1].noise_flux), np.average(BATCH[1].norm_flux)
        #----------------------------------------------------------------------

    ###################### Load GQ_flux for each spectrum ##########################
        [Synthetic.subtract_flux(ZERO) for Synthetic in BATCH]
        [Synthetic.compute_GQ() for Synthetic in BATCH]

    ################################################################################

        indices = np.array([label.GQ_index for label in BATCH])
        carbon = np.array([float(label.cfe) for label in BATCH])
        #print "Max Index: ", max(indices)
        #print "Min Index: ", min(indices)
        #print len(carbon), len(carbon)
        indices = indices[np.isfinite(indices)]
        carbon = carbon[np.isfinite(indices)]
        idx = indices.argsort()


        carbon = carbon[idx]
        indices = indices[idx]

        #print indices, carbon
        #---------------------- Fitting Functions ------------------------------

        InterpFunction = interp1d(indices, carbon)
        PolyFunction = np.poly1d(np.polyfit(indices, carbon, 3))
        Regress = linregress(indices,carbon)

        try:
            INTERP_CFE =  InterpFunction(MY_INDEX)

        except:
            INTERP_CFE = np.nan
            print "Interpolation Error!"
        try:
            POLY_CFE = PolyFunction(MY_INDEX)

        except:
            POLY_CFE = np.nan
            print "Polynomial Regression Error!"
        try:
            LINEAR_CFE = LINE(MY_INDEX, Regress[0], Regress[1])

        except:
            LINEAR_CFE = np.nan
            print "Linear Regression Error!"


        CFE_Pred_INT.append(INTERP_CFE)
        CFE_Pred_POLY.append(POLY_CFE)
        CFE_Pred_LINE.append(LINEAR_CFE)

    #    try:
    #        print "My Index", MY_INDEX
    #        print
    #        print "[C/Fe] Linear:     ", INTERP_CFE
    #        print "[C/Fe] Poly:       ", POLY_CFE
    #        print "[C/Fe] Lin Regress:", LINEAR_CFE

    #    except:
    #        print "Error Creating Fit"
    ################################################################################
    #       Plotting
    ################################################################################

        #ax = fig.add_subplot(111)
        #Synth_Indices, = ax.scatter(indices, carbon, marker="^", color="green", s=20)
        #Interp, = ax.plot(index_range, InterpFunction(index_range), linewidth=2.0, color="red", linestyle="--", label="S/N : "+str(SN))
        #Poly, = ax.plot(index_range, PolyFunction(index_range), linewidth=2.0, color="blue", linestyle="--", label="")
        #Line_Regress, = ax.plot(index_range, LINE(index_range, Regress[0], Regress[1]), linewidth=2.0, color="green", linestyle="--", label="")
        index_range = np.linspace(min(indices), max(indices), 50)

        Synth_Indices.set_data(indices,carbon)

        Interp.set_data(index_range, InterpFunction(index_range))

        Poly.set_data(index_range, PolyFunction(index_range))

        Line_Regress.set_data(index_range,LINE(index_range, Regress[0], Regress[1]))

        ax1.scatter(MY_INDEX, INTERP_CFE, color="red", s = 1.0, alpha=0.75,label="[C/Fe] Linear:  " + str(round(1000*INTERP_CFE)/1000))
        #ax1.scatter(MY_INDEX, POLY_CFE, color="blue",s=1.0,alpha=0.75, label="[C/Fe] Poly:  " + str(round(1000*POLY_CFE)/1000))
        ax1.scatter(MY_INDEX, LINEAR_CFE, color="green", s=1.0, alpha=0.75, label="[C/Fe] Line:  " + str(round(1000*LINEAR_CFE)/1000))
        #print MY_INDEX
        #ax.legend(fontsize=16)
        plt.show()
        if "p" in ARG:
            plt.pause(0.001)
        #raw_input()
        #fig.clf()


        POLY_AVG.append(np.average(np.array(CFE_Pred_POLY)[np.isfinite(np.array(CFE_Pred_POLY))]))
        INTERP_AVG.append(np.average(np.array(CFE_Pred_INT)[np.isfinite(np.array(CFE_Pred_INT))]))
        LINE_AVG.append(np.average(np.array(CFE_Pred_LINE)[np.isfinite(np.array(CFE_Pred_LINE))]))

        ax2.plot(np.arange(len(POLY_AVG)), POLY_AVG, linewidth =1.0, color="blue")
        ax2.plot(np.arange(len(INTERP_AVG)), INTERP_AVG, linewidth=1.0, color="red")
        ax2.plot(np.arange(len(LINE_AVG)), LINE_AVG, linewidth=1.0, color="green")

    CFE_Estimate = np.average((np.median(CFE_Pred_INT), np.median(CFE_Pred_POLY), np.median(CFE_Pred_LINE)))
    CFE_Sigma = max(np.median(CFE_Pred_INT), np.median(CFE_Pred_POLY), np.median(CFE_Pred_LINE)) - min(np.median(CFE_Pred_INT), np.median(CFE_Pred_POLY), np.median(CFE_Pred_LINE))
    #plt.close()
    print
    print "####################################################"
    print "                       Result:                      "
    print "[C/Fe]  Interp: ", np.median(CFE_Pred_INT), "+/-", np.std(CFE_Pred_INT)
    print "[C/Fe]  Poly:   ", np.median(CFE_Pred_POLY), "+/-", np.std(CFE_Pred_POLY)
    print "[C/Fe]  Line:   ", np.median(CFE_Pred_LINE), "+/-", np.std(CFE_Pred_LINE)
    print "[C/Fe] Average: ", np.average((np.median(CFE_Pred_INT), np.median(CFE_Pred_POLY), np.median(CFE_Pred_LINE)))
    print "####################################################"
    print

    if "p" in ARG:
        raw_input("Press any key to continue... ")
    plt.close()
    float_formatter = lambda x: "%7.6f" % x
    print float_formatter(CFE_Estimate)
    print float_formatter(CFE_Sigma)
    print " Writing result to ", spec_name
    norm_spec[0].header['CFE_EST'] = float_formatter(CFE_Estimate)
    norm_spec[0].header['CFE_ERR'] = float_formatter(CFE_Sigma)

    #try:
    norm_spec.writeto("/Users/MasterD/Google Drive/CarbonDetection/Contents/Processed/" + spec_name)
    #except:
        #Overwrite = raw_input(spec_name + " exists. Overwrite?   (y/n):   ")
        #if Overwrite == "y":

#        if True:
#            norm_spec.writeto("/Users/MasterD/Google Drive/CarbonDetection/Contents/Processed/"+spec_name)
#            print "... Complete: .fits written"
#        else:
#            print "... Complete: File not written"

    return CFE_Estimate, CFE_Sigma
