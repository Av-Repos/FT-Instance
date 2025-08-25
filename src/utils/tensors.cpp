#include "utils/tensors.h"
#include "utils/global.h"

//Initialize precision tensor
void initializePrecisionTensor(mpfr_t ** tensor, int precision, int dimensions, int permu_size){
    *tensor = (mpfr_t*)malloc(pow(permu_size,dimensions)*sizeof(mpfr_t));
    for(int i = 0; i < pow(permu_size,dimensions); i++){
        mpfr_init2((*tensor)[i], precision); 
    }
}

//Zero precision tensor
void zeroPrecisionTensor(mpfr_t * tensor, int dimensions, int permu_size){
    for(int i = 0; i < pow(permu_size,dimensions); i++) mpfr_set_d(tensor[i], 0.0, MPFR_RNDD);
}

//Copy precision tensor
void copyPrecisionTensor(mpfr_t * dest, mpfr_t * source, int dimensions, int permu_size){
  int i;
  for(i=0;i<pow(permu_size,dimensions);i++){
        mpfr_set(dest[i], source[i], MPFR_RNDD);
  }
}

//Linearly combine two precision tensors with the given coefficients
void combinePrecisionTensor(mpfr_t * dest, mpfr_t * source_a, mpfr_t * source_b, int dimensions, int permu_size, long double coef_a, long double coef_b){
    mpfr_t aux1,aux2;
    int i, precision;

    precision = dest[0]->_mpfr_prec;

    mpfr_init2(aux1, precision); 
    mpfr_init2(aux2, precision);

    for(i=0;i<pow(permu_size,dimensions);i++){
        mpfr_mul_d(aux1,source_a[i],coef_a,MPFR_RNDD);
        mpfr_mul_d(aux2,source_b[i],coef_b,MPFR_RNDD);
        mpfr_add(dest[i],aux1,aux2,MPFR_RNDD);
    }
}

//Standardize precision matrix using the given standard deviation and mean values
void customPrecisionStandardizeTensor(mpfr_t * tensor, int dimensions, int permu_size, mpfr_t mean, mpfr_t std){

    if(mpfr_cmp_d(std, Zero) > 0){
        for (int i=0; i<pow(permu_size,dimensions); i++) {
            mpfr_sub(tensor[i],tensor[i],mean,MPFR_RNDD);
            mpfr_div(tensor[i],tensor[i],std,MPFR_RNDD);
        }
    }
    
}
