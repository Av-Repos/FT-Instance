#ifndef Evaluate_MQAP
#define Evaluate_MQAP

#include "utils/global.h"

tuple<mpfr_t *,mpfr_t *> EvaluatePrecisionMQAP_metrics(mpfr_t * dist_tensor, mpfr_t * flow_tensor, int dimensions, int permu_size);

#endif
