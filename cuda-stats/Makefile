NVCC=nvcc 
INCLUDES+=-I/usr/local/NVIDIA_CUDA_SDK/common/inc
CFLAGS=-g -O3 $(INCLUDES) -D_DEBUG
LDFLAGS=$(LIBRARIES)


%.o : %.cpp 
	$(NVCC) -c $^ $(CFLAGS) -o $@ 
	$(NVCC) -M $^ $(CFLAGS)  > $@.dep 
 
%.o : %.c 
	$(NVCC) -c $^ $(CFLAGS) -o $@ 
	$(NVCC) -M $^ $(CFLAGS)  > $@.dep 
 
%.o : %.cu 
	$(NVCC) -c $^ $(CFLAGS) -o $@ 
	$(NVCC) -M $^ $(CFLAGS)  > $@.dep 

include  $(wildcard *.dep) /dev/null 
 
 
OBJECTS = cuda-stats.o
APP     = cuda-stats 
 
$(APP): $(OBJECTS) 
	$(NVCC) $(OBJECTS) $(LDFLAGS) -o $@ 
 
clean:  
	$(RM) $(OBJECTS) *.dep 
  
