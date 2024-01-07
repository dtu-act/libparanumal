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
#include <acoustics.h>
#include <acousticsTests.h>
#include "catch.hpp" // <catch2/catch.hpp>
#include <mpi.h>
#include <tinywav.h>

#include <filesystem>
#include <limits.h>
#include <stdio.h> 
#include <string.h> 
#include <vector>


namespace fs = std::filesystem;

#define BLOCK_SIZE 512

int main( int argc, char* argv[] ) {
    MPI_Init(&argc, &argv);
    int result = Catch::Session().run( argc, argv );
    MPI_Finalize();
    return result;
}

void compareAndCleanup(string filepathGen, string filepathRef, string dirGenCleanup) {

    TinyWav wavefileGen;
    tinywav_open_read(&wavefileGen, 
        filepathGen.c_str(),
        TW_INLINE
    );

    float samplesGen[BLOCK_SIZE]{0};
    tinywav_read_f(&wavefileGen, samplesGen, BLOCK_SIZE);

    TinyWav wavefileRef;
    tinywav_open_read(&wavefileRef, 
        filepathRef.c_str(),
        TW_INLINE
    );

    float samplesRef[BLOCK_SIZE]{0};
    tinywav_read_f(&wavefileRef, samplesRef, BLOCK_SIZE);

    std::vector<float> samplesGen_vec(samplesGen, samplesGen + BLOCK_SIZE);
    std::vector<float> samplesRef_vec(samplesRef, samplesRef + BLOCK_SIZE);   

    REQUIRE( samplesGen_vec == samplesRef_vec);
    
    tinywav_close_read(&wavefileGen);
    tinywav_close_read(&wavefileRef);

    // cleanup
    std::error_code ec;
    fs::remove_all((char*)dirGenCleanup.c_str(), ec);
    if (ec) {
        printf("Cleanup failed!\n"); 
    }
}

// STUDIO
TEST_CASE( "Studio with freq. indep. boundaries", "[studio][freq_indep]" ) {
    string dirGen = "tests/data/generated/studio_250hz_p4_5ppw_freq_indep/";
    string filepathGen = dirGen + "00_x0=['0.70', '1.00', '1.50']_r0=['1.00', '3.00', '1.70'].wav";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_freq_indep_00_REF.wav";

    setupAide newOptions("tests/setups/setup_studio_250hz_freq_indep");

    REQUIRE( acousticsSetupMain(newOptions) == EXIT_SUCCESS );
    
    compareAndCleanup(filepathGen, filepathRef, dirGen);
}

TEST_CASE( "Studio with freq. dep. LR boundaries", "[studio][freq_dep_lr]" ) {
    string dirGen = "tests/data/generated/studio_250hz_p4_5ppw_freq_dep_lr/";
    string filepathGen = dirGen + "00_x0=['0.70', '1.00', '1.50']_r0=['1.00', '3.00', '1.70'].wav";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_freq_dep_lr_00_REF.wav";

    setupAide newOptions("tests/setups/setup_studio_250hz_freq_dep_lr");

    REQUIRE( acousticsSetupMain(newOptions) == EXIT_SUCCESS );
    
    compareAndCleanup(filepathGen, filepathRef, dirGen);
}

TEST_CASE( "Studio with perf. refl. boundaries", "[studio][perf_refl]" ) {
    string dirGen = "tests/data/generated/studio_250hz_p4_5ppw_perf_refl/";
    string filepathGen = dirGen + "00_x0=['0.70', '1.00', '1.50']_r0=['1.00', '3.00', '1.70'].wav";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_perf_refl_00_REF.wav";

    setupAide newOptions("tests/setups/setup_studio_250hz_perf_refl");

    REQUIRE( acousticsSetupMain(newOptions) == EXIT_SUCCESS );
    
    compareAndCleanup(filepathGen, filepathRef, dirGen);
}

TEST_CASE( "Studio with combined boundaries (freq. indep. + LR + perf. refl)", "[studio][combined]" ) {
    string dirGen = "tests/data/generated/studio_250hz_p4_5ppw_combined/";
    string filepathGen = dirGen + "00_x0=['0.70', '1.00', '1.50']_r0=['1.00', '3.00', '1.70'].wav";
    string filepathRef = "tests/data/ref/studio_250hz_p4_5ppw_combined_00_REF.wav";

    setupAide newOptions("tests/setups/setup_studio_250hz_combined");

    REQUIRE( acousticsSetupMain(newOptions) == EXIT_SUCCESS );
    
    compareAndCleanup(filepathGen, filepathRef, dirGen);
}

// Cylinder
TEST_CASE( "Cylinder with perf. refl. boundaries", "[cylinder][perf_refl]" ) {
    string dirGen = "tests/data/generated/cylinder_250hz_p8_5ppw_perf_refl/";
    string filepathGen = dirGen + "00_x0=['0.00', '0.00', '0.50']_r0=['0.75', '0.00', '0.75'].wav";
    string filepathRef = "tests/data/ref/cylinder_250hz_p8_5ppw_perf_refl_00_REF.wav";

    setupAide newOptions("tests/setups/setup_cylinder_250hz_perf_refl");

    REQUIRE( acousticsSetupMain(newOptions) == EXIT_SUCCESS );
    
    compareAndCleanup(filepathGen, filepathRef, dirGen);
}

// CUBE (use this for testing purposes using argument [cube])
TEST_CASE( "Cube with freq. indep. (quick)", "[cube]" ) {
    string dirGen = "tests/data/generated/cube_500hz_p4_5ppw_freq_indep/";
    string filepathGen = dirGen + "00_x0=['0.50', '0.50', '0.50']_r0=['0.10', '0.10', '0.10'].wav";
    string filepathRef = "tests/data/ref/cube_500hz_p4_5ppw_freq_indep_00_REF.wav";

    setupAide newOptions("tests/setups/setup_cube_500hz_freq_indep");

    REQUIRE( acousticsSetupMain(newOptions) == EXIT_SUCCESS );
    
    compareAndCleanup(filepathGen, filepathRef, dirGen);
}