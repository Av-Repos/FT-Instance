#include "utils/auxiliar.h"
#include "utils/global.h"

long int multiply(int a, int b){
    long int val = 1; 
    for(int i = a; i<=b; i++) val *= i;
    return(val);
}

int index_find(int *arr, int val, int dimensions){
    for(int i = 0; i<dimensions; i++){
        if(arr[i] == val) return(i);
    }
    return(-1);
}

void unravel_index(int ind, int *unraveled_ind, int dimensions, int size) {
    for (int i = dimensions - 1; i >= 0; --i) {
        unraveled_ind[i] = ind % size;
        ind /= size;
    }
}

int unique_unravel_index(int ind, int *unraveled_ind, int dimensions, int size) {
    std::unordered_set<int> unique_set;
    for (int i = dimensions - 1; i >= 0; --i) {
        unraveled_ind[i] = ind % size;
        ind /= size;
        if(unique_set.find(unraveled_ind[i]) != unique_set.end()) return 0;
        unique_set.insert(unraveled_ind[i]);
    }
    return(1);
}

int ravel_index(int *unraveled_ind, int dimensions, int size) {
    int index = 0;
    for (int i = 0; i < dimensions; ++i)
        index = index * size + unraveled_ind[i];
    return index;
}

int unique(int* vals, int* unique_vals, int * mapping, int size) {
    std::unordered_map<int,int> unique_map;

    int i = 0;
    for (int j = 0; j<size; j++) {
        if(unique_map.find(vals[j]) == unique_map.end()){
            unique_map[vals[j]] = i;
            unique_vals[i++] = vals[j];
        }
        mapping[j] = unique_map[vals[j]];
    }

    return i;
}

void map_values(int * vals, int * mapping, int * mapped_vals, int size){
    for(int i = 0; i<size; i++){
        mapped_vals[i] = vals[mapping[i]];
    }
}

void concatenate_arrays(int* dest_arr, int* arr1, int size1, int* arr2, int size2) {

    for (int i = 0; i < size1; i++) {
        dest_arr[i] = arr1[i];
    }
    for (int i = 0; i < size2; i++) {
        dest_arr[size1 + i] = arr2[i];
    }

}

int rank_array(int * ranked_concat, int * concat, std::unordered_map<int,int> &unique_map, int size){
    int rank = unique_map.size();

    for(int i = 0; i<size; i++){
        if(unique_map.find(concat[i]) == unique_map.end()) unique_map[concat[i]] = rank++;
        ranked_concat[i] = unique_map[concat[i]];
    }

    return(rank);
}

vector<int> split(const string &s, char delimiter) {     
    vector<int> tokens;
    string token;     
    istringstream tokenStream(s);     
    while (getline(tokenStream, token, delimiter)) {
        try{
            tokens.push_back(stoi(token));
        }catch(std::exception& e){}
    }     
    return tokens;  
}
