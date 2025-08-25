#include "order.h"
#include "utils/global.h"

void copy_partition(int* dest, int* src, int size) {
    memcpy(dest, src, size * sizeof(int));
}

//Computes all "orders" (partitions) below max_order (ej., for max_order=2, (n), (n-1,1), (n-2,2) and (n-2,1,1))
tuple<int, int**, int*> get_orders(int permu_size, int max_order) {
    int** partitions = NULL;
    int* sizes = NULL;
    int count = 0;

    partitions = (int**)realloc(partitions, (count + 1) * sizeof(int*));
    sizes = (int*)realloc(sizes, (count + 1) * sizeof(int));
    
    partitions[count] = (int*)malloc(sizeof(int));
    partitions[count][0] = permu_size;
    sizes[count] = 1;
    count++;

    if (permu_size == 1) {
        return make_tuple(count, partitions, sizes);
    } else {
        for (int i = 1; i < permu_size; i++) {
            if (i > max_order) {
                break;
            }

            int sub_count;
            int** partitions_sub;
            int* sizes_sub;
            tie(sub_count, partitions_sub, sizes_sub) = get_orders(i, i);

            for (int j = 0; j < sub_count; j++) {
                if (partitions_sub[j][0] <= permu_size - i) {
                    int size = sizes_sub[j];

                    partitions = (int**)realloc(partitions, (count + 1) * sizeof(int*));
                    sizes = (int*)realloc(sizes, (count + 1) * sizeof(int));
                    
                    partitions[count] = (int*)malloc((size + 1) * sizeof(int));
                    partitions[count][0] = permu_size - i;
                    copy_partition(partitions[count] + 1, partitions_sub[j], size);

                    sizes[count] = size + 1;
                    count++;
                }
            }

            for (int j = 0; j < sub_count; j++) {
                free(partitions_sub[j]);
            }
            free(partitions_sub);
            free(sizes_sub);
        }
    }

    return make_tuple(count, partitions, sizes);
}

//Computes the Hook formula (dimension of the "order" irreducible representation, multiplicity of "order" in (1,1,...,1))
long double hook_formula(int* order, int rows, int permu_size) {

    long double coef = tgamma(permu_size + 1);

    for(int i = 0; i < rows; i++) {

        int row = order[i];
        
        for (int cell = 0; cell < row; cell++) {

            int horizontal = row - cell;

            int vertical = 0;
            for (int j = i + 1; j < rows; j++) {
                if (order[j] <= cell) {
                    break;
                }
                vertical++;
            }
            coef /= (horizontal + vertical);
        }

    }
    return coef;
}
