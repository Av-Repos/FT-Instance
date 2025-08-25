#include "decomposition.h"
#include "order.h"
#include "utils/instance.h"
#include "utils/tensors.h"
#include "utils/evaluate.h"
#include "utils/global.h"

//Main function.
int main (int argc, char *argv[]) {

    int i, j, dimensions, permu_size, order_size, MaxOrder, Accumulate, Standardize, Precision;
    mpfr_t * dist_tensor_precision, * flow_tensor_precision, * decomposed_flow_tensor_precision, * combined_flow_tensor_precision, * mean, * std;
    int ** orders, * orders_sizes;
    char * Instance;
    std::string OutputFiles, OutputFile;

    /*
     1.- Instance
     2.- Output files folder
     3.- MaxOrder
     4.- Precision
     5.- Standardize (yes or no)
     6.- Accumulate (yes or no)
     */

    if(argc != 7){
      fprintf(stderr,"Incorrect format.\n");
      exit(1);
    }

    Instance = argv[1];
    OutputFiles = argv[2];
    MaxOrder = atoi(argv[3]);
    Precision = atoi(argv[4]);
    Standardize = atoi(argv[5]);
    Accumulate = atoi(argv[6]);

    auto start = std::chrono::high_resolution_clock::now();

    //Read MQAP instance
    std::tie(dimensions, permu_size) = readPrecisionInstance(&dist_tensor_precision, &flow_tensor_precision, Precision, Instance);

    //Get orders (partitions) to be processed
    std::tie(order_size, orders, orders_sizes) = get_orders(permu_size,std::min(dimensions,MaxOrder));

    initializePrecisionTensor(&decomposed_flow_tensor_precision,Precision,dimensions,permu_size);
    initializePrecisionTensor(&combined_flow_tensor_precision,Precision,dimensions,permu_size);

    zeroPrecisionTensor(combined_flow_tensor_precision,dimensions,permu_size);

    //For each order (partition) to be processed...
    for(i = 0; i < order_size; i++){

        printf("Processing order [ ");
        for(j=0;j<orders_sizes[i];j++) printf("%d ",orders[i][j]);
        printf("]...\n");

        //Find the corresponding decomposed flow tensor from the entries in the remainder flow tensor
        fourier_precision_decomposition(decomposed_flow_tensor_precision,flow_tensor_precision,orders[i], orders_sizes[i],dimensions,permu_size);

        //Add decomposed flow tensor to combined flow tensor
        combinePrecisionTensor(combined_flow_tensor_precision, combined_flow_tensor_precision, decomposed_flow_tensor_precision, dimensions, permu_size, 1.0, 1.0);

        //If the current order has to be printed...
        if(!Accumulate || i+1 == order_size || orders[i+1][0] != orders[i][0]){
            OutputFile = OutputFiles + "/" + to_string(permu_size-orders[i][0]) + "_";
            //"Accumulate" adds the flow tensors corresponding to the same order together (ej., (n-2,2) and (n-2,1,1))
            if(!Accumulate){
                for(j=0; j < orders_sizes[i]; j++) OutputFile = OutputFile + "_" + to_string(orders[i][j]);
                printf("Writing order [ ");
                for(j=0;j<orders_sizes[i];j++) printf("%d ",orders[i][j]);
                printf("] to %s.dat.\n", OutputFile.c_str());
            }else{
                printf("Writing order %d to %s.dat.\n", permu_size-orders[i][0], OutputFile.c_str());
            }
            OutputFile = OutputFile + ".dat";
            //"Standardize" standardizes the entries in the flow tensor according to the standard deviation of the objective function
            if(Standardize){
                std::tie(mean, std) = EvaluatePrecisionMQAP_metrics(dist_tensor_precision, combined_flow_tensor_precision, dimensions, permu_size);
                //Ignore mean
                mpfr_set_d(*mean, 0.0, MPFR_RNDD);
                customPrecisionStandardizeTensor(combined_flow_tensor_precision, dimensions, permu_size, *mean, *std);
                free(mean);
                free(std);
            }
            //Print current decomposed instance
            printPrecisionInstance(dist_tensor_precision, combined_flow_tensor_precision, dimensions, permu_size, OutputFile.data());
            //Reset combined flow tensor
            zeroPrecisionTensor(combined_flow_tensor_precision,dimensions,permu_size);
        }
        
        //Subtract the current decomposed flow tensor from the remainder flow tensor
        combinePrecisionTensor(flow_tensor_precision, flow_tensor_precision, decomposed_flow_tensor_precision, dimensions, permu_size, 1.0, -1.0);
    }
    auto stop = std::chrono::high_resolution_clock::now();

    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

    cout << "Time taken by the decomposition: " << duration.count()/1000000.0 << " seconds" << endl;

    for(i = 0; i < order_size; i++){
        free(orders[i]);
    }

    free(orders);
    free(orders_sizes);

    free(dist_tensor_precision);
    free(flow_tensor_precision);
    free(decomposed_flow_tensor_precision);
    free(combined_flow_tensor_precision);
}
