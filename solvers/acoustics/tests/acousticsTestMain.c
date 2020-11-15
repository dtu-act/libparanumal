// https://github.com/catchorg/Catch2/blob/master/docs/tutorial.md#writing-tests
// http://www.electronvector.com/blog/using-catch-to-write-bdd-style-unit-tests-for-c

//#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file

#define CATCH_CONFIG_RUNNER
#include<acoustics.h>
#include<acousticsTests.h>
//#include <catch2/catch.hpp>
#include "catch.hpp"
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

void compareAndCleanup(string filepathGen, string filepathRef) {
    FILE *iFPGen = NULL;
    FILE *iFPRef = NULL;

    REQUIRE( openFile(filepathGen, &iFPGen) == 0 );
    REQUIRE( openFile(filepathRef, &iFPRef) == 0 );
        
    int line;        
    REQUIRE( comparePressureFiles(iFPGen, iFPRef, &line, 10e-8) == 0 );

    fclose(iFPGen);
    fclose(iFPRef);

    // cleanup
    if (remove((char*)filepathGen.c_str()) == 0) {
        printf("Cleanup successfully\n"); 
    }
}

// STUDIO
TEST_CASE( "Studio: Freq. indep. boundaries", "[]" ) {
    string filepathGen = "tests/data/generated/studio_250hz_p4_5ppw_freq_indep_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_freq_indep_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_studio_250hz_freq_indep");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Studio: Freq. dep. LR boundaries", "[]" ) {
    string filepathGen = "tests/data/generated/studio_250hz_p4_5ppw_freq_dep_lr_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_freq_dep_lr_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_studio_250hz_freq_dep_lr");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Studio: Perf. refl. boundaries", "[]" ) {
    string filepathGen = "tests/data/generated/studio_250hz_p4_5ppw_perf_refl_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_perf_refl_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_studio_250hz_perf_refl");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Studio: Combined boundaries (freq. indep. + LR + perf. refl)", "[]" ) {
    string filepathGen = "tests/data/generated/studio_250hz_p4_5ppw_combined_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_combined_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_studio_250hz_combined");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

// CURVILINEAR
TEST_CASE( "Cylinder: Freq. indep. boundaries", "[]" ) {
    string filepathGen = "tests/data/generated/cylinder_250hz_p8_5ppw_freq_indep_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/cylinder_250hz_p8_5ppw_freq_indep_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_cylinder_250hz_freq_indep");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Cylinder: Freq. dep. LR boundaries", "[]" ) {
    string filepathGen = "tests/data/generated/cylinder_250hz_p8_5ppw_freq_dep_lr_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/cylinder_250hz_p8_5ppw_freq_dep_lr_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_cylinder_250hz_freq_dep_lr");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Cylinder: Combined boundaries (freq. indep. + LR + perf. refl)", "[]" ) {
    string filepathGen = "tests/data/generated/cylinder_250hz_p8_5ppw_combined_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/cylinder_250hz_p8_5ppw_combined_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_cylinder_250hz_combined");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Cylinder: Perf. refl. boundaries", "[]" ) {
    string filepathGen = "tests/data/generated/cylinder_250hz_p8_5ppw_perf_refl_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/cylinder_250hz_p8_5ppw_perf_refl_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_cylinder_250hz_perf_refl");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

// TEST_CASE( "Freq. indep. all boundaries", "[cube]" ) {
//     string filepathGen = "tests/data/generated/cube_500hz_p4_5ppw_freq_indep_RecvPoint_00.txt";
//     string filepathRef = "tests/data/ref/cube_500hz_p4_5ppw_freq_indep_RecvPoint_00_REF.txt";

//     setupAide newOptions("tests/setups/setup_cube_500hz_freq_indep");

//     REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
//     compareAndCleanup(filepathGen, filepathRef);
// }