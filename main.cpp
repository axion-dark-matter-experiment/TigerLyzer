#include "flatfileinterface.h"
#include "spectrum.h"
#include "spectrumfilter.h"

#include <iostream>
#include <iomanip>      // std::setprecision

#include <chrono>

/*
//idiom for writing spectrum to file
std::ofstream ofs;
ofs.open ("/home/bephillips2/test_data.txt", std::ofstream::out | std::ofstream::app);
ofs << spectra_c;
ofs.close();
*/

//int main(int argc, char *argv[])

int main() {

    auto start = std::chrono::high_resolution_clock::now();

    Spectrum spectra;

    for (int i = 0 ; i < 1 ; i++) {
        auto Reader = FlatFileReader("/home/bephillips2/workspace/Electric_Tiger_Control_Code/data/27_20_00_20.08.2016/");

        #pragma omp parallel for ordered
        for( uint j = 0 ; j < Reader.size() ; j ++) {

            std::cout << "Loading spectrum "<< j << std::endl;
            auto spec = SingleSpectrum( Reader.at(j) ) ;

            //Note that all background subtraction steps should be perfomred -before-
            //initial binning
            if( j == 32 ) {
                plot( spec, "Single Digitized Power Spectrum");
            }

            UnsharpMask( spec, 150 );

            if( j == 32 ) {
                plot( spec, "Background Subtracted Power Spectrum");
            }
            spec.InitialBin( 72 );

            spectra += spec;
        }
    }

    //Note each spectra is implicitly converted from dBm to watts during
    //initialization, so we only need to convert to excess power

    std::cout << "Converting to units of excess power." << std::endl;
    spectra.WattsToExcessPower();

    std::cout << "Weighting spectra by expected axion power." << std::endl;
    spectra.KSVZWeight();

    std::cout << "Building grand spectra." << std::endl;
    auto g_spec = spectra.GrandSpectrum();
    plot ( g_spec, "Grand Spectrum" );

    std::cout << "Building limits." << std::endl;
    auto limits = spectra.Limits();
    plot ( limits, limits.size(), "Limits" );

    std::cout << "Building G^2 prediction." << std::endl;
    auto g_2_prediction = spectra.Limits();
    plot ( g_2_prediction, 100, "G^2 Prediction" );


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = end - start;
    auto time_taken = fp_ms.count();

    std::cout<<"Took "<<time_taken<<" ms."<<std::endl;
}
