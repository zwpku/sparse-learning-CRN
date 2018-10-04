# sparse-learning-CRN
Data-based sparsity learning of chemical reaction networks (CRNs)


## DEPENDENCIES

This package relies on the following external libraries.
   
   1.	[libconfig](https://hyperrealm.github.io/libconfig).
   	This library is used to process the configuration file: [./working_dir/sparse_learning.cfg](./working_dir/sparse_learning.cfg).

   2.	[RANLIB.C](http://www.netlib.org/random/) (file: [ranlib.c.tar.gz](http://www.netlib.org/random/ranlib.c.tar.gz)).
       	This library provides random number generators. 
	It is used in the code [./src/ssa.cpp](./src/ssa.cpp) to generate random trajectories of reactions.

To build a parallel code, we also need
   
   3.  	[MPICH](https://www.mpich.org).
	The trajectory data will be distributed and analysed by different processors. 
	Therefore, a parallel code is helpful, when there are multiple trajectories.

## COMPILE & INSTALL

1. Install the above libraries, if necessary.

2. Download the source code.

```
	git clone https://github.com/zwpku/sparse-learning-CRN.git
```

   The code should be avaiable in the directory ./sparse_learning_CRN

3. Enter the directory containing source files. 

```
  	cd ./sparse_learning_CRN/src
```

4. Make sure paths of the directories containing the header and library files are provided.  

5. Compile.

```
        make 
```
