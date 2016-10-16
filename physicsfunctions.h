#ifndef PHYSICSFUNCTIONS_H
#define PHYSICSFUNCTIONS_H

/*! \file
 * \brief Determine the width of an axion in
 * \param Frequency in MHz
 * \return
 * Estimated axion width in MHz.
 */
double axion_width ( double frequency );

/*!
 * \brief Estimate a value for \f$ g_{a \gamma\gamma} \f$ using parameters from KSVZ theory.
 * See Ed Daw's theses, page 23, eq.2.25
 *
 * \param frequency
 * frequency (in MHz)
 * \return the coupling term in GeV^-1
 */
double KSVZ_axion_coupling( double frequency );

/*!
 * \brief Lorentz line shape in terms of center frequency and Quality Factor
 *
 * \f[
 *    L(f,f_0, Q) = \frac{\Gamma^2}{ (f - f_0)^2 + \Gamma^2}\\
 *   \text{ for } \Gamma = \frac{f}{2*Q}
 * \f]
 *
 * \param f0
 * center frequency
 *
 * \param omega
 * off-center frequency
 *
 * \param Q
 * Quality Factor
 *
 * \return value of the lorentzian at point omega
 */
double lorentzian (double f0, double omega, double Q );

/*!
 * \brief Compute expected power (in Watts) from an axion to photon conversion
 * Assumes a dark matter halo density of 0.45 GeV/cm^3
 *
 * See Ed Daw's Thesis, pg. 24 eq. 2.28
 * \param effective_volume
 * Form Factor (unitless)
 * \param b_field
 * Magnetic Field (Tesla)
 * \param frequency
 * center frequency of a bin (MHz)
 * \param Q
 * Quality Factor (unitless)
 * \return
 */
double max_ksvz_power(double effective_volume, double b_field, double frequency, double Q);

/*!
 * \brief Power due to noise per bin
 * see Ed Daw's Thesis pg. 72 eq. 4.10
 *
 * \param noise_temperature
 * Noise tempearture of the entire cavity/amplifier system (Kelvin)
 *
 * \param bin_width
 * Width of bins (MHz)
 *
 * \return
 * Power in units of "watts due to noise"
 */
double power_per_bin( double noise_temperature, double bin_width );

/*!
 * \brief Convert from units of dBm to units of watts
 */
double dbm_to_watts ( double power_dbm );

#endif // PHYSICSFUNCTIONS_H
