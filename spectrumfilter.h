#ifndef SPECTRUMFILTER_H
#define SPECTRUMFILTER_H

#include "spectrum.h"

void GaussianFilter(SingleSpectrum& spec, uint radius);
void UnsharpMask(SingleSpectrum& spec, uint radius);

#endif // SPECTRUMFILTER_H
