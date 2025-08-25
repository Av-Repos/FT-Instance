#include "utils/instance.h"
#include "utils/tensors.h"
#include "utils/auxiliar.h"
#include "utils/global.h"

//Read a QAP instance from file (precision mode)
std::tuple<int, int> readPrecisionInstance(mpfr_t ** dist_tensor, mpfr_t ** flow_tensor, int precision, char * instance){
    FILE * fd;
    int dimensions, permu_size, i;
    long double a;

    fd=fopen(instance,"r");

    if(fd==NULL){
      fprintf(stderr,"The specified file doesn't exist.\n");
      exit(1);
    }

    fscanf(fd,"%d%*[^\n]",&dimensions);
    fscanf(fd,"%d%*[^\n]",&permu_size);

    initializePrecisionTensor(dist_tensor,precision,dimensions,permu_size);

    for(i=0;i<pow(permu_size,dimensions);i++){
        fscanf(fd,"%Lf\t",&a);
        mpfr_set_d((*dist_tensor)[i], (long double) a, MPFR_RNDD);
    }

    initializePrecisionTensor(flow_tensor,precision,dimensions,permu_size);

    for(i=0;i<pow(permu_size,dimensions);i++){
        fscanf(fd,"%Lf\t",&a);
        mpfr_set_d((*flow_tensor)[i], (long double) a, MPFR_RNDD);
    }

    fclose(fd);

    return(std::make_tuple(dimensions,permu_size));
}

//Print a QAP instance to file (precision mode)
void printPrecisionInstance(mpfr_t * dist_tensor, mpfr_t * flow_tensor, int dimensions, int permu_size, char * instance){
    FILE *fd;
    int i,j,k,breaks;

    if(instance!=NULL) fd=fopen(instance,"w");
    else fd = stdout;

    if(fd==NULL){
      fprintf(stderr,"The specified file doesn't exist.\n");
      exit(1);
    }

    fprintf(fd,"%d\n",dimensions);
    fprintf(fd,"%d\n",permu_size);
    for(i=0;i<pow(permu_size,dimensions);i++){
        j = i;
        breaks = 0;
        while(j % permu_size == 0 && j > 0){
            breaks += 1;
            j = std::floor(j/permu_size);
        }
        for(k=0;k<breaks;k++){
            fprintf(fd,"\n");
        }
        mpfr_fprintf(fd,"%.15RNf ",dist_tensor[i]);
    }
    for(k=0;k<dimensions;k++){
        fprintf(fd,"\n");
    }
    for(i=0;i<pow(permu_size,dimensions);i++){
        j = i;
        breaks = 0;
        while(j % permu_size == 0 && j > 0){
            breaks += 1;
            j = std::floor(j/permu_size);
        }
        for(k=0;k<breaks;k++){
            fprintf(fd,"\n");
        }
        mpfr_fprintf(fd,"%.15RNf ",flow_tensor[i]);
    }

    if(instance!=NULL) fclose(fd);
}
