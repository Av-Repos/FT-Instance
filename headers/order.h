#ifndef Order
#define Order

#include "utils/global.h"

tuple<int, int**, int*> get_orders(int highest_order, int n);
long double hook_formula(int* order, int rows, int permu_size);

#endif