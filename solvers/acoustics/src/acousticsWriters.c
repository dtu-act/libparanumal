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

#include "acousticsWriters.h"
#include "acoustics.h"
#include <fstream>
#include <vector>
#include <highfive/H5Easy.hpp>
#include <highfive/H5File.hpp>
#include "acousticsUtils.h"

using namespace HighFive;

AcousticXdmfWriter::AcousticXdmfWriter(acoustics_t *acoustics, std::vector<dfloat> timeVector)
{
  const std::string FILE_TAG = "/data";
  const std::string FILE_TAG0 = "/data0";
  const std::string FILE_TAG1 = "/data1";
  const std::string DATA_UTAG = "/udata";
  const std::string FILE_UTAG0 = "/udata0";
  const std::string FILE_UTAG1 = "/udata1";
  const std::string FILE_UTAG2 = "/udata2";

  string filenameH5 = acoustics->simulationID + ".h5";
  _filepathH5 = acoustics->outDir + "/" + filenameH5;

  // write mesh and Xdmf only once
  std::string filepathXdmf = acoustics->outDir + "/" + acoustics->simulationID + ".xdmf";
  std::string filepathRectilinearXdmf = acoustics->outDir + "/" + acoustics->simulationID + "_rectilinear.xdmf";

  auto conn_not_used = std::vector<std::vector<uint>>();
  auto x1d = std::vector<dfloat>();
  auto y1d = std::vector<dfloat>();
  auto z1d = std::vector<dfloat>();
  auto p1d = std::vector<dfloat>();
  extractUniquePoints(acoustics->mesh, acoustics, conn_not_used, x1d, y1d, z1d, p1d);

  writeXdmfHeader(acoustics, x1d.size(), filepathXdmf, filenameH5, FILE_TAG, timeVector);
  writeMeshH5(_filepathH5, FILE_TAG0, FILE_TAG1, x1d, y1d, z1d, File::Overwrite);

  H5Easy::File file(_filepathH5, File::OpenOrCreate);

  if (acoustics->x1dRectilinear.size() > 0)
  {
    auto x1d_u = acoustics->x1dRectilinear;
    auto y1d_u = acoustics->y1dRectilinear;
    auto z1d_u = acoustics->z1dRectilinear;
    writeXdmfHeader(acoustics, acoustics->x1dRectilinear.size(), filepathRectilinearXdmf, filenameH5, DATA_UTAG, {0.0});
    writeMeshH5(_filepathH5, FILE_UTAG0, FILE_UTAG1, x1d_u, y1d_u, z1d_u, File::ReadWrite);
    H5Easy::dump(file, FILE_UTAG2, acoustics->pRectilinearMesh);

    auto umesh_shape = std::vector<dfloat>({acoustics->rectilinearMeshShape[0], acoustics->rectilinearMeshShape[1], acoustics->rectilinearMeshShape[2]});
    H5Easy::dumpAttribute(file, FILE_UTAG0, "umesh_shape", umesh_shape);
  }
  if (acoustics->sourceType == GaussianFunction)
  {
    // write source position
    auto src_pos = std::vector<dfloat>({acoustics->sourcePosition[0], acoustics->sourcePosition[1], acoustics->sourcePosition[2]});
    H5Easy::dump(file, "/source_position", src_pos);
  }
}

void AcousticXdmfWriter::write(acoustics_t *acoustics, uint iter)
{
  const std::string FILE_TAG = "/data";

  auto conn = std::vector<std::vector<uint>>();
  auto x1d = std::vector<dfloat>();
  auto y1d = std::vector<dfloat>();
  auto z1d = std::vector<dfloat>();
  auto p1d = std::vector<dfloat>();

  extractUniquePoints(acoustics->mesh, acoustics, conn, x1d, y1d, z1d, p1d);

  std::string fileTag = FILE_TAG + std::to_string(iter + 2);
  H5Easy::File file(_filepathH5, File::OpenOrCreate);
  H5Easy::dump(file, fileTag, p1d);
}

void AcousticXdmfWriter::writeXdmfHeader(acoustics_t *acoustics, size_t Nelem,
                                         string fileNameXdmf, string filenameH5, string dataTag, std::vector<dfloat> timeVector)
{
  std::ofstream ofs(fileNameXdmf.c_str(), std::ofstream::out);
  // https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/data_into_visit/XdmfFormat.html#an-example-of-a-point-mesh

  // HEADER
  ofs << "<?xml version=\"1.0\" ?>" << std::endl;
  ofs << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
  ofs << "<Xdmf Version=\"3.0\">" << std::endl;

  // Begin
  ofs << "  <Domain>" << std::endl;
  ofs << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;

  for (int i = 0; i < timeVector.size(); i++)
  {
    std::string fileTag = "/" + dataTag + std::to_string(i + 2);

    ofs << "      <Grid>" << std::endl;
    ofs << "        <include xpointer=\"xpointer(//Grid[@Name=&quot;mesh&quot;]/*[self::Topology or self::Geometry])\" />" << std::endl;
    ofs << "        <Time Value=\"" << timeVector[i] << "\" />" << std::endl;
    ofs << "        <Attribute Name=\"p\" AttributeType=\"Scalar\" Center=\"Node\">" << std::endl;
    ofs << "          <DataItem DataType=\"Float\" Dimensions=\"" << Nelem << "\" Format=\"HDF\" Precision=\"8\">" << std::endl;
    ofs << "            " << filenameH5 << ":" << fileTag << std::endl;
    ofs << "          </DataItem>" << std::endl;
    ofs << "        </Attribute>" << std::endl;
    ofs << "      </Grid>" << std::endl;
  }
  ofs << "    </Grid>" << std::endl;

  ofs << "    <Grid Name=\"mesh\" GridType=\"Uniform\">" << std::endl;
  ofs << "      <Geometry GeometryType=\"XYZ\">" << std::endl;
  ofs << "        <DataItem DataType=\"Float\" Dimensions=\"" << Nelem << " " << 3 << "\" Format=\"HDF\" Precision=\"8\">" << std::endl;
  ofs << "          " << filenameH5 << ":/" + dataTag + "0" << endl;
  ofs << "        </DataItem>" << std::endl;
  ofs << "      </Geometry>" << std::endl;
  ofs << "      <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << Nelem << "\">" << std::endl;
  ofs << "        <DataItem DataType=\"Int\" Dimensions=\"" << Nelem << " " << 1 << "\" Format=\"HDF\" Precision=\"8\">" << std::endl;
  ofs << "          " << filenameH5 << ":/" + dataTag + "1" << endl;
  ofs << "        </DataItem>" << std::endl;
  ofs << "      </Topology>" << std::endl;
  ofs << "    </Grid>" << std::endl;
  ofs << "  </Domain>" << std::endl;
  ofs << "</Xdmf>" << std::endl;
}

void AcousticXdmfWriter::writeMeshH5(string fileName, string fileTag0, string fileTag1,
                                     std::vector<dfloat> &x1d, std::vector<dfloat> &y1d, std::vector<dfloat> &z1d, uint fileAttr)
{
  File file(fileName, fileAttr);

  size_t Nelem = x1d.size();

  std::vector<dfloat> rowVec(3);
  std::vector<std::vector<dfloat>> meshData(Nelem, rowVec);
  std::vector<int> vertexData(Nelem);

  // compute plot node coordinates on the fly
  for (dlong i = 0; i < x1d.size(); i++)
  {
    meshData[i][0] = x1d[i];
    meshData[i][1] = y1d[i];
    meshData[i][2] = z1d[i];
    vertexData[i] = i;
  }

  DataSet dataset = file.createDataSet<dfloat>(fileTag0, DataSpace::From(meshData));
  dataset.write(meshData);

  dataset = file.createDataSet<int>(fileTag1, DataSpace::From(vertexData));
  dataset.write(vertexData);
}

// https://github.com/BlueBrain/HighFive/blob/master/src/examples/create_extensible_dataset.cpp
AcousticH5CompactWriter::AcousticH5CompactWriter(acoustics_t *acoustics, uint nPressurePoints, std::vector<dfloat> timeVector, bool writeConnTable)
{
  const string TAG_MESH = "/mesh";
  const string TAG_PRESSURES = "/pressures";
  const string TAG_UMESH = "/umesh";
  const string TAG_UPRESSURES = "/upressures";
  const string TAG_CONN = "/conn";

  string filenameH5 = acoustics->simulationID + ".h5";
  string filepathH5 = acoustics->outDir + "/" + filenameH5;

  // write mesh only once
  auto conn = std::vector<std::vector<uint>>();
  auto x1d = std::vector<dfloat>();
  auto y1d = std::vector<dfloat>();
  auto z1d = std::vector<dfloat>();
  auto p1d = std::vector<dfloat>();
  extractUniquePoints(acoustics->mesh, acoustics, conn, x1d, y1d, z1d, p1d);

  if (writeConnTable) {
    writeMeshH5(filepathH5, TAG_MESH, TAG_CONN, x1d, y1d, z1d, conn, File::Overwrite);
  } 
  else {
    writeMeshH5(filepathH5, TAG_MESH, x1d, y1d, z1d, File::Overwrite);
  }
  
  H5Easy::File file(filepathH5, File::OpenOrCreate);

  if (acoustics->x1dRectilinear.size() > 0)
  {    
    auto x1d_u = acoustics->x1dRectilinear;
    auto y1d_u = acoustics->y1dRectilinear;
    auto z1d_u = acoustics->z1dRectilinear;
    writeMeshH5(filepathH5, TAG_UMESH, x1d_u, y1d_u, z1d_u, File::ReadWrite);
    H5Easy::dump(file, TAG_UPRESSURES, acoustics->pRectilinearMesh);

    auto umesh_shape = std::vector<dfloat>({acoustics->rectilinearMeshShape[0], acoustics->rectilinearMeshShape[1], acoustics->rectilinearMeshShape[2]});
    H5Easy::dumpAttribute(file, TAG_UMESH, "umesh_shape", umesh_shape);
  }
  if (acoustics->sourceType == GaussianFunction)
  {
    // write source position
    auto src_pos = std::vector<dfloat>({acoustics->sourcePosition[0], acoustics->sourcePosition[1], acoustics->sourcePosition[2]});
    H5Easy::dump(file, "/source_position", src_pos);
  }

  // Create a dataspace with initial shape and max shape
  DataSpace dataspace = DataSpace({timeVector.size()}, nPressurePoints);

  DataSetCreateProps props;
  props.add(Chunking(std::vector<hsize_t>{1, nPressurePoints}));
  _presssureDataset = file.createDataSet(TAG_PRESSURES, dataspace, create_datatype<dfloat>(), props);
  
  H5Easy::dumpAttribute(file, TAG_PRESSURES, "time_steps", timeVector);
}

void AcousticH5CompactWriter::write(acoustics_t *acoustics, uint iter)
{
  auto conn = std::vector<std::vector<uint>>();
  auto x1d = std::vector<dfloat>();
  auto y1d = std::vector<dfloat>();
  auto z1d = std::vector<dfloat>();
  auto p1d = std::vector<dfloat>();

  extractUniquePoints(acoustics->mesh, acoustics, conn, x1d, y1d, z1d, p1d);

  // Create a dataspace with shape [pressures, timesteps]
  _presssureDataset.select({iter, 0}, {1, p1d.size()}).write(p1d);
}

void AcousticH5CompactWriter::writeMeshH5(string fileName, string fileTag,
                                          std::vector<dfloat> &x1d, 
                                          std::vector<dfloat> &y1d, 
                                          std::vector<dfloat> &z1d, 
                                          uint fileAttr)
{
  File file(fileName, fileAttr);

  size_t Nelem = x1d.size();

  std::vector<dfloat> rowVec(3);
  std::vector<std::vector<dfloat>> meshData(Nelem, rowVec);

  // compute plot node coordinates on the fly
  for (dlong i = 0; i < x1d.size(); i++)
  {
    meshData[i][0] = x1d[i];
    meshData[i][1] = y1d[i];
    meshData[i][2] = z1d[i];
  }

  DataSet dataset = file.createDataSet<dfloat>(fileTag, DataSpace::From(meshData));
  dataset.write(meshData);
}

void AcousticH5CompactWriter::writeMeshH5(string fileName, string fileTag, string connTag,
                                          std::vector<dfloat> &x1d, 
                                          std::vector<dfloat> &y1d, 
                                          std::vector<dfloat> &z1d,
                                          std::vector<std::vector<uint>> &connData, 
                                          uint fileAttr)
{
  File file(fileName, fileAttr);

  size_t Nelem = x1d.size();

  std::vector<dfloat> rowMesh(3);
  std::vector<std::vector<dfloat>> meshData(Nelem, rowMesh);  

  // compute plot node coordinates on the fly
  for (dlong i = 0; i < x1d.size(); i++)
  {
    meshData[i][0] = x1d[i];
    meshData[i][1] = y1d[i];
    meshData[i][2] = z1d[i];
  }

  DataSet dataset = file.createDataSet<dfloat>(fileTag, DataSpace::From(meshData));
  dataset.write(meshData);

  dataset = file.createDataSet<int>(connTag, DataSpace::From(connData));
  dataset.write(connData);
}