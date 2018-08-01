
CFLAGS=-g -O2 -DMPICH_IGNORE_CXX_SEEK -DDEBUG

LOCAL_DIR=/home/mi/weizhang/local/
#LOCAL_DIR=/home/mi/weizhang/local/

LDFLAGS=-Wl,-rpath,$(LOCAL_DIR)/lib -L$(LOCAL_DIR)/lib -lm -lpthread -lconfig++ -lranlib 

INCLUDES=-I$(LOCAL_DIR)/include 

CC=mpicxx

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
