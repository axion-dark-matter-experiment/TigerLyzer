//Header for this file
#include "plotter.h"
// C System-Headers
//
// C++ System headers
//
// Boost Headers
#include <boost/algorithm/string.hpp>  //split() and is_any_of for parsing .csv files
#include <boost/lexical_cast.hpp>  //lexical cast (unsurprisingly)
// Miscellaneous Headers
#include "/home/bephillips2/gnuplot-iostream/gnuplot-iostream.h"
// Project specific headers
//

void plot ( SingleSpectrum& spec, std::string plot_title, std::string save_file_path ) {

    Gnuplot gp;

    std::vector<double> x_vals;

    double min_freq = spec.center_frequency - spec.frequency_span/2.0;
    double span = spec.frequency_span;
    double spectrum_points = static_cast<double> (spec.size());

    for( uint i = 0 ; i < spec.size() ; i ++ ) {
        double frequency = min_freq + span*i/spectrum_points;
        x_vals.push_back(frequency);
    }

    gp << "set title '"+plot_title+"'\n";
    gp << "set xlabel 'Frequency (MHz)'\n";
    gp << "set ylabel 'Power "+spec.units()+"'\n";

    std::string cmd = "plot '-' using 1:2 with lines title 'Power: "+spec.units()+"'\n";

    if( !save_file_path.empty() ) {
        std::string header = "set terminal png size 1920,1080\n";
        header += "set output '"+save_file_path+"'\n";
        cmd.insert (0,header);
    }
    gp << cmd;
    gp.send( boost::make_tuple( x_vals, spec.sa_power_list ) );
}

void plot ( SingleSpectrum& spec, uint num_plot_points, std::string plot_title, std::string save_file_path ) {

    Gnuplot gp;

    if( num_plot_points > spec.size() ) {
        std::string err_mesg = __FUNCTION__;
        err_mesg += "\nRequested number of points (";
        err_mesg += boost::lexical_cast<std::string>( num_plot_points );
        err_mesg += ") exceeds size of spectrum (";
        err_mesg += boost::lexical_cast<std::string>( spec.size() ) +")";
        throw std::out_of_range(err_mesg);
    }

    uint pivot = static_cast<uint> (spec.size() / num_plot_points);

    double min_freq = spec.center_frequency - spec.frequency_span/2.0;
    double span = spec.frequency_span;
    double spectrum_points = static_cast<double> (spec.size());

    std::vector<uint> pivot_indices;
    for( uint i = 0 ; i < spec.size() ; i ++ ) {
        if ( i%pivot == 0 ) {
            pivot_indices.push_back(i);
        }
    }

    std::vector<double> x_vals;
    for ( const auto& idx : pivot_indices ) {
        double frequency = min_freq + span*idx/spectrum_points;
        x_vals.push_back(frequency);
    }

    std::vector<double> y_vals;
    for ( const auto& idx : pivot_indices ) {
        y_vals.push_back( spec.sa_power_list.at(idx) );
    }

    std::vector<double> delta_y_vals;
    for ( const auto& idx : pivot_indices ) {
        delta_y_vals.push_back( spec.uncertainties.at(idx) );
    }

    gp << "set title '"+plot_title+"'\n";
    gp << "plot '-' using 1:2:3 with yerror title 'Power'\n";
    gp.send( boost::make_tuple(x_vals, y_vals, delta_y_vals) );
}
