#ifndef Auxiliar
#define Auxiliar

#include "utils/global.h"

long int multiply(int a, int b);
int index_find(int *arr, int val, int dimensions);
void unravel_index(int ind, int *unraveled_ind, int dimensions, int size);
int unique_unravel_index(int ind, int *unraveled_ind, int dimensions, int size);
int ravel_index(int *unraveled_ind, int dimensions, int size);
int unique(int* vals, int* unique_vals, int * mapping, int size);
void map_values(int * vals, int * mapping, int * mapped_vals, int size);
void concatenate_arrays(int* dest_arr, int* arr1, int size1, int* arr2, int size2);
int rank_array(int * ranked_concat, int * concat, std::unordered_map<int,int> &unique_map, int size);
vector<int> split(const string &s, char delimiter);

#endif
