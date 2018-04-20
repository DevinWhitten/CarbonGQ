import numpy as np
import os
import pandas as pd
from scipy.interpolate import interp1d

class Spectrum():
    def __init__(self, name, wl=0, norm_flux=0, ivar=None, flux=None,
                 SN=100, Teff=0,Logg=0, FEH=0,CFE=0, Synthetic=False):
        self.name=name

        self.wave = np.float64(wl[wl <= 5000]) # if wl is not None else list()

        self.norm_flux = np.float64(norm_flux[wl <= 5000]) # if norm_flux is not None else list()

        if not Synthetic:
            self.flux = np.float64(flux[wl <= 5000])
            self.ivar = np.float64(ivar[wl <= 5000])

        self.sn = SN
        self.teff = Teff
        self.logg = Logg
        self.feh = FEH
        self.cfe = CFE
        self.noise_flux = self.norm_flux
        self.GQ_flux = np.nan

    def compute_sn(self):
        # THIS IS INTENDED ONLY FOR THE SEGUE SPECTRA
        # Precondition: norm_flux, and wave have been loaded
        # Postcondition: computes the local S/N, as well as SN sigma
        # (4240 <-> 4254)
        self.SN = np.median(self.flux * np.sqrt(self.ivar))
        self.SN_sigma = np.std(self.flux * np.sqrt(self.ivar))
        return self.SN, self.SN_sigma

    #------------------------------

    def load_params(self):
        '''
        #Precondition: Name has been loaded
        #Postcondition: teff, logg, feh, and cfe loaded
        '''
        self.teff = self.name.split("T")[1].split("g")[0]
        self.logg = self.name.split("g")[1].split("z")[0]
        self.feh = self.name.split("z")[1].split("c")[0]
        self.cfe = self.name.split("c")[1].strip(".")

    def sn_inject(self, new_sn):
        '''
        #Precondition: wave, and flux have been loaded
        #Postcondition: throw in some gaussian noise to flux
        '''
        self.sn = new_sn
            #print sigma, self.norm_flux[I]
        self.noise_flux = np.random.normal(loc=self.norm_flux, scale=abs(self.norm_flux/self.sn)) #abs(self.norm_flux/self.sn)

        #print "Noise Injection Length:  ", len(self.noise_flux)
        return

    def subtract_flux(self, ZERO):
        '''
        Precondition: sn injection must have taken place to produce noise_flux
        Postcondition: subtracts the zero carbon spectrum, stores in GQ_flux
        '''
        if len(self.noise_flux) != len(ZERO.norm_flux):
            print "Length mismatch between current spectra and ZERO.norm_flux"
            print "Current noise flux:  ", len(self.noise_flux)
            print "ZERO.norm_flux:      ", len(ZERO.norm_flux)

        self.GQ_flux =  self.noise_flux - ZERO.norm_flux

        return

    def compute_GQ(self, width="GPE"):
        '''
        Precondition: sn injection and subtract_flux must have been run
        Postcondition: compute the GQ index based off GPE widths in Placco et al. 2010
        '''
        if width=="GPE":
            self.GQ_index = np.sum(self.GQ_flux[(self.wave > 4200.0) & (self.wave < 4400.0)])

        else: print "I haven't implemented that width yet."
        if self.GQ_index > 0:
            self.GQ_index =np.nan
            return

        self.GQ_index = np.abs(self.GQ_index)
        return

    def compute_SEGUE_GQ(self, ZERO, width="GPE"):
        '''
        Precondition: This function is reserved for SEGUE spectra ONLY! Not synthetic.
        Postcondition: performs subract_flux and compute_GQ with interpolated ZERO
        '''
        Interp_0c = interp1d(ZERO.wave, ZERO.norm_flux)
        self.GQ_flux = self.norm_flux[(self.wave > 4200.0)&(self.wave < 4400.0)] - Interp_0c(self.wave[(self.wave > 4200.0)&(self.wave < 4400.0)])

        #self.GQ_flux = self.GQ_flux[(self.wave > 4200.0)&(self.wave < 4400.0)]
        self.GQ_index = np.sum(self.GQ_flux)
        if self.GQ_index > 0.0:
                self.GQ_index = np.nan

        else:
            self.GQ_index = np.abs(self.GQ_index)

        return self.GQ_index
#----------------------------------------------------------------------

class G_index():
    '''
    Just stores the [C/Fe] value associated with the G_index
    '''
    def __init__(self, index=0, cfe=0, dif = []):
        self.index = index
        self.cfe = cfe
        self.dif = dif
    def dif_index(self):    # Just take the sum!
        self.index = np.sum(self.dif)


def LINE(x, m, b):
    return m*x + b

def Load_Synthetic(PARAMS):
    directory = os.listdir("/Users/MasterD/Google Drive/CarbonDetection/Contents/CarbonArray/")
    syn_spectra = []

    for afile in directory:
        if (len(afile.split("T")[0]) == 0) and float(afile.split("T")[1].split("g")[0]) == PARAMS[0] and np.abs(float(afile.split("T")[1].split("g")[1].split("z")[0]) - PARAMS[1]) <0.002 and np.abs(float(afile.split("T")[1].split("g")[1].split("z")[1].split("c")[0]) - PARAMS[2]) <0.002:
            syn_spectra.append(afile)

    print "Number of Synthetic Spectra:  ",len(syn_spectra)

    if len(syn_spectra) == 0:
        print "Error synthetic spectra not found for:  ", PARAMS

    SPECTRA = []
    for FILE in syn_spectra:
        DATA = pd.read_csv("/Users/MasterD/Google Drive/CarbonDetection/Contents/CarbonArray/"+FILE)
        SPECTRA.append( Spectrum(name=FILE, wl=DATA.wl, norm_flux=DATA.flux, Synthetic=True) )

    [label.load_params() for label in SPECTRA]

    ZERO, = [s for s in SPECTRA if s.cfe == "-1.5000"]   #filter(lambda x: x.cfe == "-1.5000", SPECTRA)
    BATCH = [s for s in SPECTRA if s.cfe != "-1.5000"]

    print "--------------------BATCH----------------------------"
    for element in BATCH:
        print element.cfe, np.average(element.norm_flux), element.name, len(element.noise_flux)


    if (1 + len(BATCH) != len(SPECTRA)):
        print "Error with base spectra (Mismatch)"
    else:
        print "Found base spectra:    "
        print ZERO.sn, ZERO.cfe
        print "Zero wavelength"

        return ZERO, BATCH
