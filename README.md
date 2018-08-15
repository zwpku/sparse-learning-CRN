# sparse-learning-CRN
Data-based sparsity learning of chemical reaction networks 


## DEPENDENCIES

This package relies on the following external libraries.
   
   1.	[libconfig](https://github.com/hyperrealm/libconfig).
   	This library is used to process the configuration file: [./working_dir/sparse_learning.cfg](./working_dir/sparse_learning.cfg).

   2.	[RANLIB.C](http://www.netlib.org/random/ranlib.c.tar.gz).
       	This library contains several random number generators. 
	It is used in the code ssa.cpp to generate random trajectories of reactions.

If one wants to build a parallel code, one also needs
   
   3.  	[MPICH](https://www.mpich.org).
	The trajectory data will be analysed by different processors. 
	Therefore, a parallel code is helpful, when there are multiple trajectories.

## DOWNLOAD

git clone https://github.com/zwpku/sparse-learning-CRN.git

## COMPILE & INSTALL

## LICENSING

