# FT-based instance decomposition for k-adic Assignment Problems

This repository contains the code used in the paper _Fourier Transform-based instance decomposition for k-adic Assignment Problems_.

## k-AP instance format

The accepted **k-AP** instance format, inspired by the one used in QAPLIB, is the following:

> **k** &larr; Number of dimensions in the problem.  
> **n** &larr; Instance size.  
> **D** &larr; Distance tensor of $n^k$ elements. The $k$-th dimension is separated by spaces, while each of the remaining $i$-th dimensions is separated by $k-i$ line breaks.  
> **H** &larr; Flow tensor of $n^k$ elements. The $k$-th dimension is separated by spaces, while each of the remaining $i$-th dimensions is separated by $k-i$ line breaks.  

For reference, the _data_ folder contains examples of valid **2-AP**, **3-AP** and **4-AP** instances.

## Decomposition usage

Given a valid **k-AP** instance, its decomposition can be obtained using the _bash run.sh_ command with the following flags:

  - **i**: Path of the instance to be decomposed.
  - **o**: Highest order components to be considered in the decomposition.
  - **s**: Path in which to save the generated sub-instance folder.
  - **p**: Precision for the MPFR library (256 bits by default).
  - **d**: If set, standardizes the generated sub-instances so that the standard deviation of their objective function is 1 (except for the _(n)_ component).
  - **a**: If set, accumulates generated sub-instances by order.

For example, the following command performs the decomposition of the _./data/2-AP.dat_ instance, storing the standardized sub-instances in the _./data/2-AP-decomposed_ folder:

```
bash run.sh -i ./data/2-AP.dat -o 2 -s ./data/2-AP-decomposed -d
```

Note that the sub-instance associated to the partition (a1,...,am) will be stored in a file called _a1\_...\_am.dat_.
