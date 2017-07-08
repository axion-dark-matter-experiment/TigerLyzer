#ifndef SPECTRUMFILTER_H
#define SPECTRUMFILTER_H

#include "spectrum.h"

/*!
 * \brief Convolve a SingleSpectrum object with a Gaussian kernel.
 *
 * This has the effect of sending a SingleSpectrum object through a
 * low-pass filter.
 *
 *  * Specifically this function computes the following:
 * Given a SingleSpectrum S:
 * \f[
 * S*G \text{ for }\\
 *   G[n] = \frac{1}{\sqrt{ \pi r} } e^{ \frac{2 n}{r}^2 } \text{ for }
 *   n \in [-r,r] \text{ and } r \in \mathbb{N}
 * \f]
 *
 * \param spec
 * SingleSpectrum to be convolved. Note that this object -will be modified- and
 * not copied.
 *
 * \param radius
 * The radius of the Gaussian kernel used in the convolution. The radius implicitly
 * sets the standard deviation of the kernel equal to 5/2*radius.
 */
void GaussianFilter(SingleSpectrum& spec, uint radius);

/*!
 * \brief Subtract background structure from a SingleSpectrum using
 * a high-pass filter.
 *
 * The filter uses the so-called unsharp mask technique- low-frequency
 * structure is identified by convolving the original signal with a gaussian
 * kernel. Then this static structure is subtracted from the original signal.
 *
 * Specifically this function computes the following:
 * Given a SingleSpectrum S:
 * \f[
 * S - S*G \text{ for }\\
 *   G[n] = \frac{1}{\sqrt{ \pi r} } e^{ \frac{2 n}{r}^2 } \text{ for }
 *   n \in [-r,r] \text{ and } r \in \mathbb{N}
 * \f]
 *
 * Edges effects are handled by mirroring edge-points.
 *
 * \param spec
 * The SingleSpectrum that requires background subtraction. Note that the
 * spectrum -will be modified-, this function does not return a new spectrum.
 *
 * \param radius
 * The radius of the Gaussian kernel used in the convolution. The radius implicitly
 * sets the standard deviation of the kernel equal to 5/2*radius.
 */
void UnsharpMask(SingleSpectrum& spec, uint radius, double sigma );

/*!
 * \brief Compute the ideal radius to be used for background subtraction using the UnsharpMask function.
 *
 * We assume each spectra we collect will be a combination of white noise, static structure and (perhaps)
 * and axion signal. We wish to remove the static structure from the rest of the signal using the UnsharpMask
 * function. We need to be careful to remove as much of the static structure as possible without removing any of the actual
 * signal. Since we assume our noise is white noise we know that our signal should follow the relationship that
 * \f$ \frac{\mu}{\sigma} \approx \frac{1}{\sqrt{N}} \f$ where N is the total number of points in the signal.
 * This function will compute what parameter, as used by the UnsharpMask function, will yield a signal
 * that is closest to satisfying the previous statement.
 *
 *
 * \param spec
 * The spectrum that will later have static structure removed.
 *
 * \return
 * The kernel radius (as used by the UnsharpMask function) that will yield optimal background subtraction.
 */
std::pair< uint, double > AutoOptimize( SingleSpectrum& spec, uint max_radius, uint max_sigma );

#endif // SPECTRUMFILTER_H
