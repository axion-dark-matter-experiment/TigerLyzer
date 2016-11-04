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
#include "singlespectrum.h"

/*!
 * \brief Plot a SingleSpectrum object as a connected line graph
 *
 * Note that uncertainities will -not- be plotted, only power values
 * - further every power value will be plotted. If the spectrum contains a large
 *  number of points plotting will be very slow or fail.
 *
 * Units will be automatically be added to the graph, based on whatever the
 *  units the SingleSpectrum is currently in (tracked automatically).
 *
 * \param spec
 * The SingleSpectrum class to be plotted
 *
 * \param plot_title
 * The title that will appear on the graph
 */
void plot ( SingleSpectrum& spec, std::string plot_title, std::string save_file_path = "" );


/*!
 * \brief Plot a SingleSpectrum object as a point graphic with error bars
 *
 * Error bars are automatically generated using the uncertainities associated
 * with SingleSpectrum.
 *
 * Units will be automatically be added to the graph, based on whatever
 * units the SingleSpectrum is currently in (tracked automatically).
 *
 * \param spec
 * The SingleSpectrum class to be plotted
 *
 * \param num_plot_points
 * The number of elements (and associated uncertainties) that should be plotted.
 * This number should be less than or equal to the size of the SingleSpectrum,
 *  which can be found with the SingleSpectrum.size() function.
 *
 * \param plot_title
 * The title that will appear on the graph
 */
void plot ( SingleSpectrum& spec, uint num_plot_points, std::string plot_title, std::string save_file_path = "");

#endif // PLOTTER_H
