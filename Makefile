CXX = g++
NVCC = nvcc

all:
	@echo "Usage: make [target]"
	@echo "Available targets:"
	@echo "  make serial      Compile serial code"
	@echo "  make cuda        Compile cuda code"
	@echo "  make gif         Generate output.gif based on output.txt"
	@echo "  make clean       Remove executable file and output.txt"

serial: fss_serial.cpp main.cpp
	$(CXX) fss_serial.cpp main.cpp -o serial -Wall
# 	./serial -s 1 -o output.txt

cuda: fss_cuda.cu main.cu
	$(NVCC) -O3 fss_cuda.cu main.cu -o fss_cuda
# 	./fss_cuda -s 1 -o output.txt

openmp: fss_openmp.cpp main.cpp
	$(CXX) fss_openmp.cpp main.cpp -o openmp -Wall -fopenmp

gif: output.txt render.py
	python3 render.py output.txt output.gif 0.05

clean:
	rm -f serial
	rm -f fss_cuda
	rm -f output.txt