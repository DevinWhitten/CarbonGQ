
# Just throw the carbon_array.csv to nuit\
# Maybe later I will try to run the idl interp.pro

scp ../Contents/carbon_array.csv dwhitten@nuit.phys.nd.edu:/emc1/home/dwhitten/CarbonDetection/carbon_array.csv

ssh dwhitten@nuit.phys.nd.edu << LETS_TRY

    cd /emc1/home/dwhitten/CarbonDetection
    rm CarbonArray/*
    idl << ANOTHERLEVEL
        interp
        exit
ANOTHERLEVEL
    zip -r CarbonArray.zip CarbonArray/
    exit

LETS_TRY

rm  -r ~/Google\ Drive/CarbonDetection/Contents/CarbonArray

#scp dwhitten@nuit.phys.nd.edu:/emc1/home/dwhitten/CarbonDetection/CarbonArray/* /Users/MasterD/Google\ Drive/CarbonDetection/Contents/Spectra/

scp dwhitten@nuit.phys.nd.edu:/emc1/home/dwhitten/CarbonDetection/CarbonArray.zip /Users/MasterD/Google\ Drive/CarbonDetection/Contents

unzip -o ~/Google\ Drive/CarbonDetection/Contents/CarbonArray.zip -d ~/Google\ Drive/CarbonDetection/Contents/
