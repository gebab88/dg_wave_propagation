#verbose with -v flag
#
UNAME_S := $(shell uname -s)

ifeq ($(origin CXX), default)
  ifeq ($(UNAME_S),Darwin)
    CXX = clang++
  else
    CXX = g++
  endif
endif

DEBUGOROPTI		?= -O3 #-g -O0

# Simulation parameters (order, N, time integrator, fluid/domain settings) are
# now read at RUNTIME from config.yaml via yaml-cpp — they are no longer
# compile-time -D macros, so no rebuild is needed to change them.

BASE_CFLAGS = -Wall $(DEBUGOROPTI) -std=c++14 -pthread -Iinclude -DARMA_DONT_USE_WRAPPER

ifeq ($(UNAME_S),Darwin)
  # MacPorts paths. Armadillo headers are under /opt/local/include.
  MACPORTS ?= /opt/local
  INCLUDES = -I$(MACPORTS)/include

  # OpenMP via MacPorts libomp (Apple clang needs -Xpreprocessor -fopenmp and the
  # versioned libomp paths). Used to parallelize the per-cell loop in ClassdXdt.
  OMPFLAGS = -Xpreprocessor -fopenmp -I$(MACPORTS)/include/libomp
  OMPLIBS  = -L$(MACPORTS)/lib/libomp -lomp

  # Armadillo 15 needs at least C++14. ARMA_DONT_USE_WRAPPER -> we link
  # BLAS/LAPACK ourselves via the Accelerate framework below.
  CFLAGS  = $(BASE_CFLAGS) -m64 $(INCLUDES) $(OMPFLAGS)

  # HDF5 isn't used (all hdf5 calls in the code are commented out), so we don't
  # link it. Boost / gnuplot-iostream were dropped too: the solver talks to gnuplot
  # directly (script files + the gnuplot process, see ClassPlot), so neither is needed.
  # yaml-cpp (MacPorts) is linked for runtime config parsing (config.yaml).
  LDFLAGS = -framework Accelerate -pthread $(OMPLIBS) -L$(MACPORTS)/lib -lyaml-cpp
else ifeq ($(UNAME_S),Linux)
  # Linux builds use the distro BLAS/LAPACK packages. In the Docker test image
  # those are OpenBLAS + LAPACK, and g++ supplies OpenMP via -fopenmp.
  CFLAGS  = $(BASE_CFLAGS) -fopenmp
  LDFLAGS = -pthread -fopenmp -lopenblas -llapack -lyaml-cpp
else
  $(error Unsupported OS "$(UNAME_S)"; set CFLAGS and LDFLAGS explicitly)
endif

#LDFLAGS = -L/home/georg/Downloads/armadillo-7.900.1/libarmadillo.so
#LDFLAGS = -larmadillo -lboost_iostreams -lboost_system -lboost_filesystem  -llapack -lopenblas

#LDFLAGS = -larmadillo  -DARMA_DONT_USE_WRAPPER -llapack -lopen-blas

VIDEO = output
PLAYER = vlc
UNAME:=$(shell uname -a)

# Explicit object list so the per-file rules below are actually triggered.
# (The old `obj/Release/*.o` glob matched nothing on a clean tree and broke the build.)
OBJS= obj/Release/main.o obj/Release/InitialAndBoundary.o obj/Release/ClassdXdt.o obj/Release/vandermonde.o obj/Release/ClassPlot.o
BIN= bin/Release/DiscontinousGalerkin

# Rebuild objects when any header changes (the rules below don't track #includes
# otherwise, so e.g. editing include/ClassdXdt.hpp would silently not recompile).
HEADERS= $(wildcard include/*.h include/*.hpp)

# Output directories are gitignored, so make sure they exist on a fresh clone.
$(shell mkdir -p obj/Release bin/Release bin/Debug)

Release:  $(OBJS)
	$(CXX)  -o $(BIN) $(OBJS) $(LDFLAGS)

Debug: $(OBJS)
	$(CXX) -o bin/Debug/DiscontinousGalerkin $(OBJS) $(LDFLAGS)

obj/Release/main.o: src/main.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c $< -o obj/Release/main.o

obj/Release/InitialAndBoundary.o: src/InitialAndBoundary.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c $< -o obj/Release/InitialAndBoundary.o

obj/Release/ClassdXdt.o: src/ClassdXdt.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c $< -o obj/Release/ClassdXdt.o

obj/Release/vandermonde.o: src/vandermonde.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c $< -o obj/Release/vandermonde.o

obj/Release/ClassPlot.o: src/ClassPlot.cpp $(HEADERS)
	$(CXX) $(CFLAGS) -c $< -o obj/Release/ClassPlot.o


video:
	#@echo 'Programm is runnung on System' $(UNAME)
	ffmpeg  -framerate 60 -i %08d.png -c:v libx264 -pix_fmt yuv420p -r 30 $(VIDEO).mp4
	#ffmpeg -framerate 30 -i *.png -c:v libx264 -r 30 $(VIDEO).mp4
	open -a $(PLAYER) $(VIDEO).mp4

cleanOBJ:
	rm -rf $(OBJS)

cleanBIN:
	rm -rf bin/Release/DiscontinousGalerkin

clean:
	rm -rf $(OBJS)
	rm -rf $(BIN)
	rm -rf *.png
	rm -rf *.mp4
	rm -rf *.dat

cleanPictures:
	rm -rf *.png

docker-test-linux-arm64:
	scripts/docker-linux-arm64-test.sh

.PHONY: cleanOBJ cleanBIN cleanAll video cleanPictures Release Debug docker-test-linux-arm64
