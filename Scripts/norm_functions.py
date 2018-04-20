#Author: Devin Whitten
#Date: Nov 12, 2016
# This is will serve as the interface for the normalization function.
# So just defining some functions in here.

import pyfits
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as interp


class Segment():
    def __init__(self, wl=[], flux=[], mad=0, median=0,midpoint=0):
        self.wl = wl
        self.flux = flux
        self.mad=0
        self.median=0
        self.midpoint = midpoint

    def get_midpoint(self):
        self.midpoint = (self.wl[-1] + self.wl[0])/2.0
        #self.wl.iloc[-1]
    def get_stats(self):
        # Just calculate the robust statistics for the segment
        self.mad = np.median(np.absolute(self.flux - np.median(self.flux)))
        self.median = np.median(self.flux)

    def sigma_clip(self,low_cut, high_cut):
        #Postcondition: Conditions segment by excluding values outside of
        #               high_cut and low_cut tolerance.
        # Now do the clip, rejust points beyond specific sigma
        TOP = self.median + high_cut * self.mad
        BOTTOM = self.median - low_cut * self.mad

        self.wl = self.wl[(self.flux < TOP) & (self.flux > BOTTOM)]
        self.flux = self.flux[(self.flux < TOP) & (self.flux > BOTTOM)]

#-------------------------------------------------------------------------------

def find_cont(WL, FLUX, Kernals=30, low=1.0, high=1.0):
    # Precondition:
    # Postcondition: Divides spectra into <Kernal> segments,
    #               computes flux midpoint of SEGMENTS
    #               returns wl and flux values for continuum determination

    if(len(WL) != len(FLUX)):
        print "Mismatch"
        return

    Len = len(WL)
    wave_out, flux_out  = [], []

    Spectra = pd.DataFrame({"Wave":WL, "Flux":FLUX})

    SEGMENTS = []
    for BIN in np.array_split(Spectra, Kernals):
        SEGMENTS.append(Segment(wl=np.array(BIN.Wave), flux=np.array(BIN.Flux)))

    for Seg in SEGMENTS:
        Seg.get_stats()  #just loads the median and mad values for each segment.
        Seg.sigma_clip(low,high)  # Rejects flux beyond thresholds for each.

    for Seg in SEGMENTS:
        # Here we just rebuild the global wavelength and flux arrays.
        wave_out = np.concatenate((wave_out,Seg.wl), axis=0)
        flux_out = np.concatenate((flux_out,Seg.flux), axis=0)

    print "Continuum Points:   ", len(wave_out), len(flux_out)
    return wave_out, flux_out

#-------------------------------------------------------------------------------

class Spectra():
    # Houses the normalizing functions poly_normalize, spline_normalize
    def __init__(self, wave=[], flux=[], trim_wave=[], trim_flux=[], Continuum=[], Norm=[], Function="none"):
        self.wave = np.array(wave, dtype=np.float64)
        self.flux = np.array(flux, dtype=np.float64)
        self.trim_wave = trim_wave
        self.trim_flux = trim_flux
        self.continuum = Continuum
        self.norm_flux = Norm
        self.Function = Function

    def poly_normalize(self, nlow=3.0, nhigh=3.0, boost=0.05, order=4):
        self.Poly_WL, self.Poly_FLUX = find_cont(self.wave, self.flux, low=nlow, high=nhigh)
        self.Poly_FLUX = np.array(self.Poly_FLUX + boost*self.Poly_FLUX)
        self.Poly_WL = np.array(self.Poly_WL)

        print "Poly Find_cont:  ", len(self.Poly_WL), len(self.Poly_FLUX)
        self.Function = np.poly1d(np.polyfit(self.Poly_WL[(self.Poly_WL < 4200.) | (self.Poly_WL >4400.)], self.Poly_FLUX[(self.Poly_WL < 4200.) | (self.Poly_WL >4400.)], order))
        self.Poly_Continuum = self.Function(self.wave)  # We can get error information here
        self.norm_flux = np.divide(self.flux, self.Poly_Continuum)
        print "... Poly Normalize Complete!"

    def spline_normalize(self, nlow=1.0, nhigh=5.0, BINS=25, boost=0.05, G_exclude=True):

        # Trim Spectra, bin for spline interpolation
        WL, FLUX = find_cont(self.wave, self.flux, low = nlow, high = nhigh)

        Spectra = pd.DataFrame({"Wave":WL, "Flux":FLUX})

#        if G_exclude:
#            LEFT = Spectra[(Spectra.Wave < 4200.0)]
#            RIGHT = Spectra[(Spectra.Wave > 4400.0)]
#            Spectra = pd.concat([LEFT, RIGHT])

        SEGMENTS = []

        for BIN in np.array_split(Spectra, BINS):
            #print np.average(BIN.Wave)
            SEGMENTS.append(Segment(wl=np.array(BIN.Wave), flux=np.array(BIN.Flux)))
            #print "Elements:   ",len(BIN.Wave)
        self.wl_Bin, self.flux_Bin = [], []
        wl_BIN, flux_BIN = [], []

        for kernal in SEGMENTS:
            kernal.get_midpoint()
            kernal.get_stats()
            wl_BIN.append(kernal.midpoint)
            flux_BIN.append(kernal.median + boost*kernal.median)

        wl_BIN = np.array(wl_BIN)
        flux_BIN = np.array(flux_BIN)

        self.wl_Bin =  wl_BIN[(wl_BIN < 4100.0) | (wl_BIN > 4500.0)]
        self.flux_Bin = flux_BIN[(wl_BIN < 4100.0) | (wl_BIN > 4500.0)]
        #flux_BIN[np.all([wl_BIN < 4200.0, wl_BIN >4400.0], axis=0)]

        #print "---------------"
        #self.wl_Bin =  wl_BIN[np.where(wl_BIN < 4200.) & np.where(wl_BIN > 4400.)]
        #self.flux_Bin = flux_BIN[np.where(wl_BIN < 4200.0) & np.where(wl_BIN > 4400.)]

        #print "wl_BIN:  ", self.wl_Bin
        #print 'flux_BIN:  ', self.flux_Bin

        print "Flux Bin:  ",np.average(self.flux_Bin)
        tck = interp.splrep(self.wl_Bin, self.flux_Bin,k=3,s=0.0)
        self.Spline_Continuum = interp.splev(self.wave, tck)
        self.norm_flux = np.divide(self.flux, self.Spline_Continuum)
        print "Spline complete"


    #def linear_normalize(self, )
    #

#-------------------------------------------------------------------------------
