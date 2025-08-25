#ifndef Problem_Instance
#define Problem_Instance

#include "utils/global.h"

std::tuple<int, int> readPrecisionInstance(mpfr_t ** dist_tensor, mpfr_t ** flow_tensor, int precision, char * instance);
void printPrecisionInstance(mpfr_t * dist_tensor, mpfr_t * flow_tensor, int dimensions, int permu_size, char * instance);

#endif
