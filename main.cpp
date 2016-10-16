#include "flatfileinterface.h"
#include "spectrum.h"
#include "singlespectrum.h"
#include "spectrumfilter.h"
#include "plotter.h"
#include "physicsfunctions.h"

#include <iostream>
#include <iomanip>      // std::setprecision

#include <chrono>

/*! \mainpage Tigerlyzer: An analysis program for the Electric Tiger experiment
 *
 * \section intro_sec Introduction
 *
 * This program is designed to interpret data collected by the Electric Tiger experiment
 * and generate exclusion limits on \f$ g_{a\gamma\gamma} ( GeV^{-1} ) \f$.
 *
 * \section Dependencies
 *      \li boost/algorithm/string
 *      \li OpenMP ( and Compilier that supports OpenMP #pragma's such as GCC)
 *      \li gnuplot-iostream (included with program files)
 *
 */

int main() {

    auto start = std::chrono::high_resolution_clock::now();

    Spectrum spectra;

    for (int i = 0 ; i < 1 ; i++) {
        auto Reader = FlatFileReader("/home/bephillips2/workspace/Electric_Tiger_Control_Code/data/27_20_00_20.08.2016/", "SA_F");
//        auto Reader = FlatFileReader("/home/bephillips2/workspace/Electric_Tiger_Control_Code/data/09_56_11_17.08.2016/");

        for( uint j = 0 ; j < Reader.size() ; j ++) {

            std::cout << "Loading spectrum "<< j << std::endl;
            auto spec = SingleSpectrum( Reader.at(j) ) ;

            //Note that all background subtraction steps should be perfomred -before-
            //initial binning
            if( j == 20 ) {
                plot( spec, "Single Digitized Power Spectrum");
            }

            UnsharpMask( spec, 10 );

            if( j == 20 ) {
                plot( spec, "Background Subtracted Power Spectrum");
            }
            spec.InitialBin( 32 );

            spectra += spec;
        }
    }

    //Note each spectra is implicitly converted from dBm to watts during
    //initialization, so we only need to convert to excess power

    std::cout << "Converting to units of excess power." << std::endl;
    spectra.WattsToExcessPower();

    auto e_spec = spectra.at(20);
    plot( e_spec, "Excess Power Spectra" );

    std::cout << "Weighting spectra by expected axion power." << std::endl;

    spectra.LorentzianWeight();
    spectra.KSVZWeight();

    auto ax_spec = spectra.at(20);
    plot( ax_spec, "Axion Power Spectra" );

    std::cout << "Building grand spectra." << std::endl;
    auto g_spec = spectra.GrandSpectrum();
    plot ( g_spec, "Grand Spectrum" );

    std::cout << "Building limits." << std::endl;
    auto limits = spectra.Limits();
    plot ( limits, "Limits" );

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = end - start;
    auto time_taken = fp_ms.count();

    std::cout<<"Took "<<time_taken<<" ms."<<std::endl;
}
