# sparse-learning-CRN
Data-based sparsity learning of chemical reaction networks 


## DEPENDENCIES

This package relies on the following external libraries :
   
   1.	libconfig 

   	Webpage: https://github.com/hyperrealm/libconfig

   	This library is used to process the configuration file :  ./working_dir/sparse_learning.cfg

   2.	RANLIB.C  

   	Webpage: http://www.netlib.org/random/ranlib.c.tar.gz 

       	This library contains several random number generators. 	
       	It is used in the code ssa.cpp to generate random trajectories of reactions.

 In order to build a parallel code, the package also needs 
   
   3.  MPICH

	The trajectory data will be analysed by different processors. Therefore, a parallel code is helpful, when there are multiple trajectories.

## COMPILE & INSTALL

## LICENSING

