all : GPU2DDS
CUDA_PATH = /usr/local/cuda
CC = mpicxx
NVCC = nvcc
CFLAGS = -O3 -w -fopenmp
NVINC = -I$(CUDA_PATH)/include
NVCCFLAGS = -O3 -w 
#NVCCFLAGS = -O3 -gencode arch=compute_37,code=sm_37 -w
ORDER = -DORDER=12
#MACRO = -DELASTIC -DPP  
#MACRO = -DELASTIC -DVECTORIMG -DPP -DPOYNTING 
MACRO = -DELASTIC
LDFLAGS = -L$(CUDA_PATH)/lib64 -lcudart -lm -lcublas -lcufft

GPU2DDS : FD2DGPU.o main.o laplace.o born.o
	$(CC) -o $@ $^ $(NVINC) $(LDFLAGS) $(CFLAGS)

%.o : %.cu 
	$(NVCC)  $(NVINC) $(NVCCFLAGS) $(ORDER) $(MACRO) -o $@ -c $^  

%.o: %.cpp 
	$(CC)  $(NVINC) $(CFLAGS) $(MACRO) -o $@ -c $^

clean:
	rm -f main.o FD2DGPU.o laplace.o born.o GPU2DDS

