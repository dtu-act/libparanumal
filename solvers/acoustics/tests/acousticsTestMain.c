/*

The MIT License (MIT)

Copyright (c) 2017 Tim Warburton, Noel Chalmers, Jesse Chan, Ali Karakus

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

*/

// https://github.com/catchorg/Catch2/blob/master/docs/tutorial.md#writing-tests

// TO SPECIFY WHICH TEST CASE TO RUN: https://github.com/catchorg/Catch2/blob/v2.x/docs/command-line.md#specifying-which-tests-to-run
// e.g.: mpirun ./acousticsTestMain  [studio][freq_indep]

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
    if (remove((char*)filepathGen.c_str()) != 0) {
        printf("Cleanup failed!\n"); 
    }
}

// STUDIO
TEST_CASE( "Studio with freq. indep. boundaries", "[studio][freq_indep]" ) {
    string filepathGen = "tests/data/generated/studio_250hz_p4_5ppw_freq_indep_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_freq_indep_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_studio_250hz_freq_indep");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Studio with freq. dep. LR boundaries", "[studio][freq_dep_lr]" ) {
    string filepathGen = "tests/data/generated/studio_250hz_p4_5ppw_freq_dep_lr_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_freq_dep_lr_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_studio_250hz_freq_dep_lr");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Studio with perf. refl. boundaries", "[studio][perf_refl]" ) {
    string filepathGen = "tests/data/generated/studio_250hz_p4_5ppw_perf_refl_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_perf_refl_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_studio_250hz_perf_refl");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Studio with combined boundaries (freq. indep. + LR + perf. refl)", "[studio][combined]" ) {
    string filepathGen = "tests/data/generated/studio_250hz_p4_5ppw_combined_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_combined_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_studio_250hz_combined");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

// CURVILINEAR
TEST_CASE( "Cylinder with freq. indep. boundaries", "[cylinder][freq_indep]" ) {
    string filepathGen = "tests/data/generated/cylinder_250hz_p8_5ppw_freq_indep_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/cylinder_250hz_p8_5ppw_freq_indep_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_cylinder_250hz_freq_indep");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Cylinder with freq. dep. LR boundaries", "[cylinder][freq_dep_lr]" ) {
    string filepathGen = "tests/data/generated/cylinder_250hz_p8_5ppw_freq_dep_lr_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/cylinder_250hz_p8_5ppw_freq_dep_lr_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_cylinder_250hz_freq_dep_lr");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Cylinder with perf. refl. boundaries", "[cylinder][perf_refl]" ) {
    string filepathGen = "tests/data/generated/cylinder_250hz_p8_5ppw_perf_refl_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/cylinder_250hz_p8_5ppw_perf_refl_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_cylinder_250hz_perf_refl");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

TEST_CASE( "Cylinder with combined boundaries (freq. indep. + LR + perf. refl)", "[cylinder][combined]" ) {
    string filepathGen = "tests/data/generated/cylinder_250hz_p8_5ppw_combined_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/cylinder_250hz_p8_5ppw_combined_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_cylinder_250hz_combined");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}

// not run by default - testing purposes
TEST_CASE( "Cube with freq. indep. (quick)", "[!hide][cube]" ) {
    string filepathGen = "tests/data/generated/cube_500hz_p4_5ppw_freq_indep_RecvPoint_00.txt";
    string filepathRef = "tests/data/ref/cube_500hz_p4_5ppw_freq_indep_RecvPoint_00_REF.txt";

    setupAide newOptions("tests/setups/setup_cube_500hz_freq_indep");

    REQUIRE( acousticsSetupMain(newOptions) == 0 );
    
    compareAndCleanup(filepathGen, filepathRef);
}