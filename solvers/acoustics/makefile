# https://www2.cs.duke.edu/courses/cps108/doc/makefileinfo/sample.html
# compiling and linking libraries: https://tldp.org/HOWTO/Program-Library-HOWTO/static-libraries.html

ifndef OCCA_DIR
ERROR:
	@echo "Error, environment variable [OCCA_DIR] is not set"
endif

#CXXFLAGS = -O0 # no optimization, for debugging
CXXFLAGS = -O2

include ${OCCA_DIR}/scripts/Makefile

# define variables
HDRDIR  = ../../include
OGSDIR  = ../../libs/gatherScatter
LIBSDIR = ../../libs

# set options for this machine
# specify which compilers to use for c, fortran and linking
CC	= mpic++ -std=c++2a
LD	= mpic++ -std=c++2a

# compiler flags to be used (set to compile with debugging on)
CFLAGS = -I. -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -I$(HDRDIR) -I$(OGSDIR) -g -D DHOLMES='"${CURDIR}/../.."' -D DACOUSTICS='"${CURDIR}"'

# if Gaussian Random Fields are included with #define INCLUDE_GRF 1 in acoustics.h, remember to add -lbtwxt to LIBS (see README)
OS := $(shell uname)
ifeq ($(OS),Darwin)
	# link flags to be used (MacOS)
	LDFLAGS	= -Xpreprocessor -fopenmp -lomp -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -g
else
	# link flags to be used (Linux)
	LDFLAGS	= -DOCCA_VERSION_1_0 $(compilerFlags) $(flags) -g -lhdf5 -larmadillo	
endif

# libraries to be linked in
#LIBS	=   -L$(OCCA_DIR)/lib -L$(LIBSDIR) -lbtwxt $(links) # use this if GRFs should be computed (sample initial condition on uniform grid)
LIBS	=   -L$(OCCA_DIR)/lib $(links) # otherwise use this

INCLUDES = acoustics.h

DEPS = $(INCLUDES) \
$(HDRDIR)/mesh.h \
$(HDRDIR)/mesh3D.h

# types of files we are going to construct rules for
.SUFFIXES: .c 

# rule for .c files
.c.o: $(DEPS)
	$(CC) $(CFLAGS) -o $*.o -c $*.c $(paths) 

# list of objects to be compiled
OBJS    = \
./src/acousticsEstimate.o \
./src/acousticsReceiver.o \
./src/acousticsStep.o \
./src/acousticsError.o \
./src/acousticsRun.o \
./src/acousticsSetup.o \
./src/acousticsSetupMain.o \
./src/acousticsSources.o \
./src/acousticsWriters.o \
./src/acousticsIO.o \
./src/acousticsUtils.o \
../../src/meshParallelReaderTet3DCurv.o \
../../src/meshSetupTet3DCurv.o \
../../src/meshGeometricPartition3DCurv.o \
../../src/meshLoadInterpolationCurv.o \
../../src/meshGeometricFactorsTet3DCurv.o \
../../src/meshSurfaceGeometricFactorsTet3DCurv.o \
../../src/meshConnect.o \
../../src/meshConnectBoundary.o \
../../src/meshConnectFaceNodes2D.o \
../../src/meshConnectFaceNodes3D.o \
../../src/meshGeometricFactorsTet3D.o \
../../src/meshGeometricFactorsHex3D.o \
../../src/meshGeometricFactorsTri2D.o \
../../src/meshGeometricFactorsQuad2D.o \
../../src/meshGeometricPartition2D.o \
../../src/meshGeometricPartition3D.o \
../../src/meshHaloExchange.o \
../../src/meshHaloExtract.o \
../../src/meshHaloSetup.o \
../../src/meshLoadReferenceNodesTri2D.o \
../../src/meshLoadReferenceNodesQuad2D.o \
../../src/meshLoadReferenceNodesTet3D.o \
../../src/meshLoadReferenceNodesHex3D.o \
../../src/meshOccaSetup2D.o \
../../src/meshOccaSetup3D.o \
../../src/meshParallelConnectNodes.o \
../../src/meshParallelConnectOpt.o \
../../src/meshParallelPrint2D.o \
../../src/meshParallelReaderTri2D.o \
../../src/meshParallelReaderQuad2D.o \
../../src/meshParallelReaderTet3D.o \
../../src/meshParallelReaderHex3D.o \
../../src/meshPartitionStatistics.o \
../../src/meshPhysicalNodesTri2D.o \
../../src/meshPhysicalNodesQuad2D.o \
../../src/meshPhysicalNodesTet3D.o \
../../src/meshPhysicalNodesHex3D.o \
../../src/meshPlotVTU2D.o \
../../src/meshPlotVTU3D.o \
../../src/meshPrint2D.o \
../../src/meshPrint3D.o \
../../src/meshSetupTri2D.o \
../../src/meshSetupQuad2D.o \
../../src/meshSetupTet3D.o \
../../src/meshSetupHex3D.o \
../../src/meshSurfaceGeometricFactorsTri2D.o \
../../src/meshSurfaceGeometricFactorsQuad2D.o \
../../src/meshSurfaceGeometricFactorsTet3D.o \
../../src/meshSurfaceGeometricFactorsHex3D.o \
../../src/meshVTU2D.o \
../../src/meshVTU3D.o \
../../src/mysort.o \
../../src/parallelSort.o \
../../src/matrix.o \
../../src/setupAide.o \
../../src/trace.o \
../../src/readArray.o \
../../src/occaDeviceConfig.o\
../../src/occaHostMallocPinned.o \
../../src/timer.o

OBJS_MAIN    = \
./src/acousticsMain.o

OBJS_TESTS    = \
./tests/acousticsTestMain.o \
./tests/utils/comparePressureFiles.o \
./tests/utils/compareFiles.o \
./tests/utils/openFile.o

acousticsMain:$(OBJS) $(OBJS_MAIN)
	$(LD)  $(LDFLAGS)  -o acousticsMain $(OBJS) $(OBJS_MAIN) $(paths) $(LIBS)

tests:$(OBJS) $(OBJS_TESTS) 
	$(LD)  $(LDFLAGS)  -o acousticsTestMain $(OBJS) $(OBJS_TESTS) $(paths) $(LIBS)

clean :
	rm -r $(OBJS) $(OBJS_MAIN) acousticsMain
	rm -r $(OBJS) $(OBJS_TESTS) acousticsTestMain
