#ifndef PLOTTER_H
#define PLOTTER_H

// C System-Headers
//
// C++ System headers
//
// Boost Headers
//
// Miscellaneous Headers
//
// Project specific headers
#include "spectrum.h"

//Friend class of SingleSpectrum
void plot ( SingleSpectrum& spec, std::string plot_title );
//Friend class of SingleSpectrum
void plot ( SingleSpectrum& spec, uint num_plot_points, std::string plot_title );

#endif // PLOTTER_H
