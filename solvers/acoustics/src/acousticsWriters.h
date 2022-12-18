#ifndef ACOUSTICS_WRITERS_H
#define ACOUSTICS_WRITERS_H

#include "acoustics.h"
#include <fstream>
#include <vector>
#include "acousticsUtils.h"
#include <highfive/H5Easy.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>

class IAcousticWriter
{
public:
    virtual void write(acoustics_t *acoustics, uint iter) = 0;
};

class AcousticXdmfWriter : public IAcousticWriter
{
public:
    AcousticXdmfWriter(acoustics_t *acoustics, std::vector<dfloat> timeVector);
    void write(acoustics_t *acoustics, uint iter);

private:
    std::string _filepathH5;
    void writeXdmfHeader(acoustics_t *acoustics, size_t Nelem,
                         std::string fileNameXdmf, std::string filenameH5, std::string dataTag, std::vector<dfloat> timeVector);

    void writeMeshH5(std::string fileName, std::string fileTag0, std::string fileTag1,
                     std::vector<dfloat> &x1d, std::vector<dfloat> &y1d, std::vector<dfloat> &z1d, uint fileAttr);
};

// https://github.com/BlueBrain/HighFive/blob/master/src/examples/create_extensible_dataset.cpp
class AcousticH5CompactWriter : public IAcousticWriter
{
public:
    AcousticH5CompactWriter(acoustics_t *acoustics, uint nPressurePoints, std::vector<dfloat> timeVector,
        bool writeConnTable=false);
    void write(acoustics_t *acoustics, uint iter);

private:
    HighFive::DataSet _presssureDataset; // constructor is called - dummy
    void writeMeshH5(std::string fileName, std::string fileTag0,
                     std::vector<dfloat> &x1d, std::vector<dfloat> &y1d, std::vector<dfloat> &z1d, uint fileAttr);
    void writeMeshH5(std::string fileName, std::string fileTag0, std::string fileTag1,
                     std::vector<dfloat> &x1d, std::vector<dfloat> &y1d, std::vector<dfloat> &z1d, 
                     std::vector<std::vector<uint>> &conn, uint fileAttr);
};

#endif