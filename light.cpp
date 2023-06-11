#include "light.h"
#include <math.h>

double Light::GetExcitationFrequency(const double ne, const double Te, const HeLvls Term) const{
	const double beta = z * z * Ry / Te;
	int idx = static_cast<int>(Term);
	return ne * 1e-8 * A[idx] / (z*z*z * (beta + Chi[idx]))
			* sqrt(beta / (1 + Elvl[idx] / Te)) * (1 + beta + D[idx]) * exp(- Elvl[idx] / Te); // unit [1/s]
}
