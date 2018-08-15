
# 
# If USE_MPI=0, build a sequential code ;  
# Otherwise, build a parallel code using MPI
#
USE_MPI = 0

ifeq ($(USE_MPI), 0) 
  CFLAGS=-g -O2 -DDEBUG 
  CC=g++
else
  CFLAGS=-g -O2 -DMPICH_IGNORE_CXX_SEEK -DDEBUG -DUSE_MPI
  CC=mpicxx
endif

#
# This package relies on the external libraries libconfig and RANLIB.C  
#
# In order to build a parallel code, the package also needs the MPICH library.
# 
# The directories containing the header files and library files of these
# external libraries should be provided.
#
# In particular, please specify 
#
# 	1. the directory containing libconfig.h++. This is required by the library libconfig.
#
# 	2. the directory containing ranlib.h, This is required by the library RANLIB.C.
#
# and, if a parallel version is wanted (USE_MPI=1), 
#
# 	3. the directory containing mpi.h, required by MPICH.
#
LOCAL_DIR=/home/wzhang/local/
INCLUDES=-I$(LOCAL_DIR)/include 

LDFLAGS=-Wl,-rpath,$(LOCAL_DIR)/lib -L$(LOCAL_DIR)/lib -lm -lpthread -lconfig++ -lranlib 

SOURCES=param.cpp config.cpp initialize.cpp utils.cpp ssa.cpp sparse_learning.cpp prepare.cpp

OBJECTS=$(SOURCES:.cpp=.o)

.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

all : ssa sparse_learning prepare 

ssa : ssa.o param.o config.o initialize.o utils.o
	$(CC) $^ -o $@ $(CFLAGS) $(LDFLAGS)

sparse_learning : sparse_learning.o param.o config.o initialize.o utils.o
	$(CC) $^ -o $@ $(CFLAGS) $(LDFLAGS)

prepare : prepare.o param.o config.o initialize.o utils.o
	$(CC) $^ -o $@ $(CFLAGS) $(LDFLAGS)

clean:
	rm -f $(OBJECTS) ssa sparse_learning prepare 
