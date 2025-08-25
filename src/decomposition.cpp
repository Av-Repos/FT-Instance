#include "decomposition.h"
#include "order.h"
#include "utils/auxiliar.h"
#include "utils/global.h"

//Computes the psi multiplicity value (number of times the term d_i1...ik*h_p1...pk appears in the summation) (Precision)
mpfr_t * psi_precision(int * order, int order_size, int * indexes, int * values, int precision, int dimensions, int permu_size){
    std::vector<std::unordered_set<int>> tabloid_1(order_size), tabloid_2(order_size);
    std::unordered_map<int,int> used;
    long int num,den;
    mpfr_t * psi_val;
    mpfr_t aux;
    int * mapping;
    int map_k, ind_k, val_k, cont;

    psi_val = (mpfr_t *)malloc(sizeof(mpfr_t));

    mpfr_init2(*psi_val, precision);
    mpfr_set_d(*psi_val, 0.0, MPFR_RNDD);

    mpfr_init2(aux, precision);

    mapping = (int *)malloc(dimensions*sizeof(int));
    
    //For all possible (j_1...j_k) such that j_a=1,...,|order|...
    //(All possible distributions of k elements in |order| rows)
    for(int j=0;j<pow(order_size,dimensions);j++){

        cont = 1;
        used.clear();

        for(int k = 0; k < order_size; k++){
            tabloid_1[k].clear();
            tabloid_2[k].clear();
        }

        unravel_index(j, mapping, dimensions, order_size);
        
        //Generate a Young tabloid t_1 of shape "order" such that i_a \in t_1(j_a)
        //Generate a Yount tabloid t_2 of shape "order" such that i_a,p_a \in t_2(j_a)
        for(int k=0;k<dimensions;k++){

            map_k = mapping[k];
            ind_k = indexes[k];
            val_k = values[k];

            tabloid_1[map_k].insert(ind_k);
            tabloid_2[map_k].insert(ind_k);
            tabloid_2[map_k].insert(val_k);
            
            //If it is not possible to generate t_2, continue to next iteration
            if(tabloid_2[map_k].size() > (long unsigned int) order[map_k] || 
                (used.find(ind_k) != used.end() && used[ind_k] != map_k) ||
                (used.find(val_k) != used.end() && used[val_k] != map_k)){
                
                cont = 0;
                break;
          
            }

            used[ind_k] = map_k;
            used[val_k] = map_k;
        }

        if(!cont) continue;

        //Compute:
        //
        // ( Number of tabloids that meet the conditions of t_2 ) <---- N_tab
        //                         *
        // ( Number of permutations that meet pi({t_1})=t_2 and pi(i_a)=p_a ) <----- N_perm
        //
        //We avoid unnecesary factorial operations dividing this value by n! and simplifying unnecesary terms
        num = 1;
        for(int k=0;k<order_size;k++) num *= multiply(order[k]-tabloid_2[k].size()+1,order[k]-tabloid_1[k].size());
        den = multiply(permu_size-used.size()+1,permu_size);

        //Update psi_value as psi_value+=N_tab*N_perm
        mpfr_set_ui(aux,(unsigned long int) num, MPFR_RNDD);
        mpfr_div_ui(aux,aux,(unsigned long int) den,MPFR_RNDD);
        mpfr_add(*psi_val,*psi_val,aux,MPFR_RNDD);
    }

    free(mapping);

    return(psi_val);
}

//Computes all possible keys for the cache that stores the psi_values
std::vector<string> get_psi_keys(int dimensions){

    std::vector<string> keys;
    std::vector<string> nums;
    std::vector<string> nums_;
    std::vector<string> nums__;

    //For 1,2,...,dimensions number of unique indexes/values
    for(int i = 1; i <= dimensions; i++){

        //Get all possible individual values in the keys for the current number of unique indexes/values
        nums.clear();
        for(int j=-1;j<i;j++)
            nums.push_back(" " + std::to_string(j));
        
        nums_ = nums;

        //Generate all posible keys for the current number of unique indexes/values
        for(int j = 1; j < i; j++){
            nums__.clear();
            for (string a: nums_) {
                for (int k = 0; k <= i; k++){
                    //Only the -1 value can appear multiple times in the keys
                    if(a.find(nums[k]) == string::npos || nums[k] == " -1") nums__.push_back(a + nums[k]);
                }
            }
            nums_ = nums__;
        }

        keys.insert(keys.end(), nums_.begin(), nums_.end());
    }

    return(keys);
}

//Computes the values in the decomposed flow tensor directly from the values in the remainder flow tensor (Precision)
void fourier_precision_decomposition(mpfr_t * decomposed_flow_tensor, mpfr_t * flow_tensor, int * order, int order_size, int dimensions,
                            int permu_size){

    int * indexes_unique_gen, * values_unique_gen;
    int limit, number, max_number, precision, unique_dimensions_gen;
    mpfr_t * aux_psi;
    std::unordered_map<std::string, mpfr_t*> cache;
    std::vector<string> keys;

    precision = flow_tensor[0]->_mpfr_prec;

    indexes_unique_gen = (int *)malloc(dimensions*sizeof(int));
    values_unique_gen = (int *)malloc(dimensions*sizeof(int));
    
    //Compute hook formula
    long int hook = (long int) (hook_formula(order, order_size, permu_size)+0.5);

    //Compute the cache values to be reutilized later (all possible values of psi_value)

    //Get all possible keys for the cache
    keys = get_psi_keys(dimensions);

    //For each possible key in the cache...
    for(string key: keys){

        std::istringstream iss(key);

        //Get unique index/value arrays for the current key
        unique_dimensions_gen = 0;
        while (iss >> number) {
            values_unique_gen[unique_dimensions_gen] = unique_dimensions_gen;
            indexes_unique_gen[unique_dimensions_gen] = number;
            unique_dimensions_gen++;
        }

        //Replace -1 by actual possible values
        max_number = unique_dimensions_gen;

        for(int i = 0; i < unique_dimensions_gen; i++){
            if(indexes_unique_gen[i] == -1){
                indexes_unique_gen[i] = max_number;
                max_number += 1;
            }
        }

        //If the given arrays are not possible, continue with the next key
        if(max_number > permu_size) continue;

        //Compute psi_value and store it in cache
        aux_psi = (mpfr_t *)malloc(sizeof(mpfr_t));
        mpfr_init2(*aux_psi, precision);
        mpfr_set(*aux_psi, *psi_precision(order,order_size,indexes_unique_gen,values_unique_gen,precision,unique_dimensions_gen,permu_size),MPFR_RNDD);
        cache[key] = aux_psi;
    }

    //For each position in the decomposed flow tensor... (Parallel computation of the for loop)

    #if defined(_OPENMP)
    	omp_set_num_threads(omp_get_max_threads());
    #endif
    
    #pragma omp parallel
    {
        int * indexes, * values, * indexes_unique, * values_unique, * mapping;
        mpfr_t aux_flow;
        std::string key;
        int p, unique_dimensions;

        mpfr_init2(aux_flow, precision);

        indexes = (int *)malloc(dimensions*sizeof(int));
        values = (int *)malloc(dimensions*sizeof(int));
        indexes_unique = (int *)malloc(dimensions*sizeof(int));
        values_unique = (int *)malloc(dimensions*sizeof(int));
        mapping = (int *)malloc(dimensions*sizeof(int));

        limit = pow(permu_size,dimensions);
    
        #pragma omp for
        for(int i=0;i<limit;i++){
            
            mpfr_set_d(decomposed_flow_tensor[i], 0.0, MPFR_RNDD);

            unravel_index(i, indexes, dimensions, permu_size);

            //Get unique indexes to avoid unnecesary computations. Get mapping for the repeated indexes
            unique_dimensions = unique(indexes,indexes_unique,mapping,dimensions);

            //For each position in the remainder flow tensor (avoiding unnecesary computations)...
            for(int p_=0;p_<pow(permu_size,unique_dimensions);p_++){

                //Check if multiple unique indexes have the same associated values. If that is the case, continue to next iteration
                //(psi_value=0)
                if(unique_unravel_index(p_, values_unique, unique_dimensions, permu_size)){

                    //Map unique values according to repeated indexes
                    map_values(values_unique,mapping,values,dimensions);
                    p = ravel_index(values, dimensions, permu_size);

                    //Generate key for the cache storage. The key indicates in which position of unique values appears each unique index
                    //For example: indexes = (1,2) values = (0,2) -> key = "-1 1" (-1 indicates that index does not appear in values)
                    key="";
                    for (int k = 0; k < unique_dimensions; ++k) key += " " + std::to_string(index_find(values_unique, indexes_unique[k], unique_dimensions));

                    //Add remainder_flow_tensor[p]*psi_value to decomposed_flow_tensor[i]
                    mpfr_mul(aux_flow,flow_tensor[p],*(cache[key]),MPFR_RNDD);
                    mpfr_add(decomposed_flow_tensor[i],decomposed_flow_tensor[i],aux_flow,MPFR_RNDD);
                }

            }
            //Multiply decomposed_flow_tensor[i] by the hook value
            mpfr_mul_ui(decomposed_flow_tensor[i],decomposed_flow_tensor[i],(unsigned long int) hook,MPFR_RNDD);
        }

        free(indexes);
        free(indexes_unique);
        free(values);
        free(values_unique);
        free(mapping);
    }

    free(indexes_unique_gen);
    free(values_unique_gen);
    for (auto& pointers : cache) {
        free(pointers.second);
    }
}
