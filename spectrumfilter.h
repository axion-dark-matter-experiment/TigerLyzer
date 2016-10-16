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
void UnsharpMask(SingleSpectrum& spec, uint radius);

#endif // SPECTRUMFILTER_H
