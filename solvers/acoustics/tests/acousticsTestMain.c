// https://github.com/catchorg/Catch2/blob/master/docs/tutorial.md#writing-tests
// http://www.electronvector.com/blog/using-catch-to-write-bdd-style-unit-tests-for-c

//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#define CATCH_CONFIG_RUNNER
#include<acoustics.h>
#include<acousticsTests.h>
#include <catch2/catch.hpp>
#include <limits.h>
#include<stdio.h> 
#include<string.h> 
#include <mpi.h>

int main( int argc, char* argv[] ) {
    MPI_Init(&argc, &argv);
    int result = Catch::Session().run( argc, argv );
    MPI_Finalize();
    return result;
}

TEST_CASE( "Simulation in cube, freq. indep. boundaries", "[factorial]" ) {
    string filepathGen = "tests/data/generated/cube_500hz_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/cube_500hz_RecvPoint_00_REFERENCE.txt";

    setupAide newOptions("tests/setups/setup_cube");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    FILE *iFPGen = NULL;
    FILE *iFPRef = NULL;

    REQUIRE( openFile(filepathGen, &iFPGen) == 0 );
    REQUIRE( openFile(filepathRef, &iFPRef) == 0 );
        
    int line, col;        
    REQUIRE( compareFiles(iFPGen, iFPRef, &line, &col) == 0 );

    // cleanup
    if (remove((char*)filepathGen.c_str()) == 0) {
        printf("Cleanup successfully"); 
    }      
}