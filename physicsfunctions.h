#ifndef PHYSICSFUNCTIONS_H
#define PHYSICSFUNCTIONS_H

/*! \file
 * \brief Determine the width of an axion in
 * \param Frequency in MHz
 * \return
 */
double axion_width ( double frequency );

/*! \file
 * \brief Don't ask, I don't know. Tends to be ~10^-15
 * \param Frequency in MHz
 * \return
 */
double KSVZ_axion_coupling( double frequency );


double estimate_G_2 ( double freq_mhz );

/*! \file
 * \brief Lorentz Line Shape
 * \param f0, the center frequency
 * \param frequency of interest
 * \param Q, quality factor
 * \return value of lorenzian at omega
 */
double lorentzian (double f0, double omega, double Q );

/*!
 * \brief Get axion-photon conversion power for a Quality Factor of 1
 * Assumes dark mattery density of 0.45 GeV/cm^3
 * \param effective volume of cavity in cm^3
 * \param Magnetic field in Tesla
 * \param Frequency in GHz
 * \return photon conversion power in watts (?)
 */
double max_ksvz_power( double effective_volume, double b_field, double frequency);

//double noise_power=kB*runs[onrun].noise_temperature *(thespectrum.getBinWidth()*1e6);

/*! \file
 * \brief Determine the amount power due to noise
 * \param Noise temperature of the cavity and all attached
 * amplifiers
 * \param The width of the bin in MHz
 * \return Estimated power due to nosie
 */
double power_per_bin( double noise_temperature, double bin_width );

double dbm_to_watts ( double power_dbm );

#endif // PHYSICSFUNCTIONS_H
