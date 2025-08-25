#ifndef Decomposition
#define Decomposition

#include "utils/global.h"

void fourier_precision_decomposition(mpfr_t * decomposed_flow_tensor, mpfr_t * flow_tensor, int * order, int order_size, int dimensions,
                            int permu_size);

#endif