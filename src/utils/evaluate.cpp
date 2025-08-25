#include "utils/instance.h"
#include "utils/evaluate.h"
#include "utils/auxiliar.h"
#include "utils/global.h"

//Computes mean and standard deviation of the objective function of a MQAP over the entire search space (Precision)
tuple<mpfr_t *,mpfr_t *> EvaluatePrecisionMQAP_metrics(mpfr_t * dist_tensor, mpfr_t * flow_tensor, int dimensions, int permu_size){
    int * indexes, * indexes_unique, * values, * values_unique, * mapping, * ranked;
    int p, unique_dimensions, precision;
    mpfr_t * aux_point, * std, * mean;
    mpfr_t mean_2, aux1, aux2;
    long int div;
    std::vector<int> key_vec, key_vec1, key_vec2;
    std::unordered_map<int,int> used_partial;
    std::unordered_map<int,int> used;
    std::unordered_map<std::string,mpfr_t*> distances;
    std::unordered_map<std::string,mpfr_t*> flows;
    std::unordered_map<std::string,long int> divs;
    std::string key,key_r,key_partial;

    precision = flow_tensor[0]->_mpfr_prec;

    mean = (mpfr_t *)malloc(sizeof(mpfr_t));
    std = (mpfr_t *)malloc(sizeof(mpfr_t));

    indexes = (int *)malloc(dimensions*sizeof(int));
    indexes_unique = (int *)malloc(dimensions*sizeof(int));
    values = (int *)malloc(dimensions*sizeof(int));
    values_unique = (int *)malloc(dimensions*sizeof(int));
    mapping = (int *)malloc(dimensions*sizeof(int));
    ranked = (int *)malloc(2*dimensions*sizeof(int));
    
    mpfr_init2(*mean, precision);
    mpfr_set_d(*mean,0.0,MPFR_RNDD);
    mpfr_init2(mean_2, precision);
    mpfr_set_d(mean_2,0.0,MPFR_RNDD);

    mpfr_init2(*std, precision);
    mpfr_init2(aux1, precision);
    mpfr_init2(aux2, precision);

    //Compute mean of the objective function
    //For each position in the distance tensor...
    for(int i = 0; i < pow(permu_size,dimensions); i++){

        key = "";
        used.clear();

        unravel_index(i, indexes, dimensions, permu_size);

        //Get unique indexes to avoid unnecesary computations. Get mapping for the repeated indexes
        unique_dimensions = unique(indexes,indexes_unique,mapping,dimensions);

        //Generate key for the cache storage. The key "normalizes" the indexes to give all indexes with the same characteristics
        //an identical codification. For example: indexes_1 = (0,3,4,0), indexes_2 = (2,1,0,2) ---> key = (0,1,2,0)
        rank_array(ranked,indexes,used,dimensions);

        for(int k = 0; k < dimensions; k++) key += std::to_string(ranked[k]) + " ";

        //If the key is not in the cache storage...
        if(flows.find(key) == flows.end()){

            aux_point = (mpfr_t *)malloc(sizeof(mpfr_t));
            mpfr_init2(*aux_point, precision);
            mpfr_set_d(*aux_point,0.0,MPFR_RNDD);

            flows[key] = aux_point;
       
            //For each position in the flow tensor (avoiding unnecesary computations)...
            for(int p_=0;p_<pow(permu_size,unique_dimensions);p_++){

                //Check if multiple unique indexes have the same associated values. If that is the case, continue to next iteration
                if(unique_unravel_index(p_, values_unique, unique_dimensions, permu_size)){

                    //Map unique values according to repeated indexes
                    map_values(values_unique,mapping,values,dimensions);
                    p = ravel_index(values, dimensions, permu_size);

                    //Add flow_tensor[p] to the cache storage
                    mpfr_add(*(flows[key]),*(flows[key]),flow_tensor[p],MPFR_RNDD);
                
                }

            }
        }

        //Compute div value for the current number of unique indexes in i
        div = multiply(permu_size-unique_dimensions+1,permu_size);

        //Add dist_tensor[i]*(flow_tensor[p]+...)/div to the mean 
        mpfr_mul(aux1,dist_tensor[i],*(flows[key]),MPFR_RNDD);
        mpfr_div_ui(aux1,aux1,(unsigned long int)div,MPFR_RNDD);

        mpfr_add(*mean,*mean,aux1,MPFR_RNDD);
    }

    flows.clear();

    //Compute mean of the squared objective function
    //For each position in the distance/flow tensor...
    for(int i = 0; i < pow(permu_size,dimensions); i++){

        //Generate partial key for the cache storage. The key "normalizes" the indexes to give all indexes with the same 
        //characteristics an identical codification. For example: indexes_1 = (0,3,4,0), indexes_2 = (2,1,0,2) ---> key = (0,1,2,0)
        key_partial = "";
        used_partial.clear();
        
        unravel_index(i, indexes, dimensions, permu_size);

        rank_array(ranked,indexes,used_partial,dimensions);

        for(int k = 0; k < dimensions; k++) key_partial += std::to_string(ranked[k]) + " ";
        
        //For each position in the distance/flow tensor...
        //Note that, due to the symmetric nature of our computations, we ensure that p>=i to avoid unnecessary iterations
        for(p = i; p < pow(permu_size,dimensions); p++){

                //Generate remaining key for the cache storage. The key "normalizes" the values to give all values with the same 
                //characteristics an identical codification. For example: values_1 = (0,3,4,0), values_2 = (2,1,0,2) ---> key = (0,1,2,0)
                //Note that we take into account the elements that appear in the "normalized" indexes. That is,
                //indexes = (0,3,4,0) -> partial_key = (0,1,2,0) , values = (0,0,3,5) -> key = (0,1,2,0|0,0,1,3)
                key = key_partial;
                used = used_partial;
               
                unravel_index(p, values, dimensions, permu_size);

                div = rank_array(ranked,values,used,dimensions);

                for(int k = 0; k < dimensions; k++) key += std::to_string(ranked[k]) + " ";

                //If the key is not in the cache storage...
                if(distances.find(key) == distances.end()){

                    //Initialize cache storage
                    aux_point = (mpfr_t *)malloc(sizeof(mpfr_t));
                    mpfr_init2(*aux_point, precision);
                    mpfr_set_d(*aux_point,0.0,MPFR_RNDD);

                    distances[key] = aux_point;

                    aux_point = (mpfr_t *)malloc(sizeof(mpfr_t));
                    mpfr_init2(*aux_point, precision);
                    mpfr_set_d(*aux_point,0.0,MPFR_RNDD);

                    flows[key] = aux_point;

                    //Compute div value for the current number of unique indexes/values in i+p
                    divs[key] = multiply(permu_size-div+1,permu_size);

                    //If i==p, divide both dist and flow by 2 to avoid duplicated computations
                    if(i==p)  divs[key] *= 4;
                }

                //Add dist_tensor[i]*dist_tensor[p] to the dist cache storage
                mpfr_mul(aux1,dist_tensor[i],dist_tensor[p],MPFR_RNDD);
                mpfr_add(*(distances[key]),*(distances[key]),aux1,MPFR_RNDD);

                //Add flow_tensor[i]*flow_tensor[p] to the flow cache storage
                mpfr_mul(aux1,flow_tensor[i],flow_tensor[p],MPFR_RNDD);
                mpfr_add(*(flows[key]),*(flows[key]),aux1,MPFR_RNDD);

        }
    }

    //For each entry in the cache storage...
    for (const auto& el : distances) {

        //Compute inverse key to account for the symmetric nature of our computations
        //For example: key = (0,0,1,2|1,1,0,0) ----> key_r = (0,0,1,1|1,1,0,2)
        key_r = "";
        used.clear();

        key_vec = split(el.first, ' ');
        key_vec1 = std::vector<int>(key_vec.begin() + dimensions, key_vec.end());
        key_vec2 = std::vector<int>(key_vec.begin(), key_vec.begin() + dimensions);
        key_vec = key_vec1;
        key_vec.insert(key_vec.end(),key_vec2.begin(),key_vec2.end());

        rank_array(ranked,key_vec.data(),used,2*dimensions);

        for(int k = 0; k < 2*dimensions; k++) key_r += std::to_string(ranked[k]) + " ";

        //Add dist_cache[key]+dist_cache[key_r] to account for the symmetric nature of our computations
        mpfr_add(aux1,*(distances[el.first]),*(distances[key_r]),MPFR_RNDD);
        //Add flow_cache[key]+flow_cache[key_r] to account for the symmetric nature of our computations
        mpfr_add(aux2,*(flows[el.first]),*(flows[key_r]),MPFR_RNDD);

        //Add dist_cache*flow_cache/div to the squared mean
        mpfr_mul(aux1,aux1,aux2,MPFR_RNDD);
        mpfr_div_ui(aux1,aux1,(unsigned long int)divs[el.first],MPFR_RNDD);
        mpfr_add(mean_2,mean_2,aux1,MPFR_RNDD);
    }

    //Compute standard deviation of the objective function as sqrt(squared_mean-mean**2)
    mpfr_pow_ui(aux1,*mean,2,MPFR_RNDD);
    mpfr_sub(*std,mean_2,aux1,MPFR_RNDD);

    if(mpfr_cmp_d(*std,0.0) > 0) mpfr_sqrt(*std,*std,MPFR_RNDD);
    else mpfr_set_d(*std,0.0,MPFR_RNDD);

    free(indexes);
    free(indexes_unique);
    free(values);
    free(values_unique);
    free(mapping);
    free(ranked);

    return(make_tuple(mean,std));
}
