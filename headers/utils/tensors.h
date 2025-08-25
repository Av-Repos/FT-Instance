#ifndef Tensors
#define Tensors

#include "utils/global.h"

void initializePrecisionTensor(mpfr_t ** tensor, int precision, int dimensions, int permu_size);
void zeroPrecisionTensor(mpfr_t * tensor, int dimensions, int permu_size);
void copyPrecisionTensor(mpfr_t * dest, mpfr_t * source, int dimensions, int permu_size);
void combinePrecisionTensor(mpfr_t * dest, mpfr_t * source_a, mpfr_t * source_b, int dimensions, int permu_size, long double coef_a, long double coef_b);
void customPrecisionStandardizeTensor(mpfr_t * tensor, int dimensions, int permu_size, mpfr_t mean, mpfr_t std);

#endif
