#Author: Devin Whitten


# This function will serve as the wrapper for the carbon routine

# First executes normalize.py, to generate normalized fits file.

# Next runs get_spec_params to load the synthetic spectra array

# Finally, runs the interpolate.py routine to compute the GQ index
# and carbon estimate.

import integrate, get_spec_params, normalize, ExecuteFunctions
import subprocess, os
import pandas as pd
import numpy as np
import sys
#normalize.normalize("spec-2187-54270-0564.fits")

#### Let's test.

ARG = ExecuteFunctions.Process_Arguments(sys.argv)

SPEC_LIST = os.listdir("/Users/MasterD/Google Drive/CarbonDetection/Contents/SEGUE")
# First element in SPEC_LIST is always weird thing.
print " ... Found ", len(SPEC_LIST[1:]), " in Contents/SEGUE"

INITIAL_PARAMS = [get_spec_params.get_spec_params(SPECTRUM) for SPECTRUM in SPEC_LIST[1:]]
#Generate one spectra request for entire SEGUE batch
get_spec_params.generate_carbon_array(INITIAL_PARAMS)  # TEFF, LOGG, FEH, CFE_COR

#Push to nuit

subprocess.call("./to_nuit.sh", shell=True)
print "NOT WORKING"
_ = [normalize.normalize(SPECTRUM, ARG) for SPECTRUM in SPEC_LIST[1:]]
CFE_ESTIMATE = [integrate.integrate(SPECTRUM, PARAMS, ARG =ARG) for SPECTRUM, PARAMS in zip(SPEC_LIST[1:], INITIAL_PARAMS)]
#INITIAL_PARAMS, CFE_ESTIMATE = np.matrix(INITIAL_PARAMS), np.matrix(CFE_ESTIMATE)
print "-------------------"

Output_Frame = pd.DataFrame({"SPSPEC":SPEC_LIST[1:],
                            "TEFF":[element[0] for element in INITIAL_PARAMS],
                            "LOGG": [element[1] for element in INITIAL_PARAMS],
                            "FEH": [element[2] for element in INITIAL_PARAMS],
                            "CFE_COR": [element[3] for element in INITIAL_PARAMS],
                            "CFE_ERR": [element[4] for element in INITIAL_PARAMS],
                            "CFE_EST":[CFE[0] for CFE in CFE_ESTIMATE],
                            "CFE_Sigma": [CFE[1] for CFE in CFE_ESTIMATE]
                            })

Output_Frame.to_csv("Result.csv", index=False)
