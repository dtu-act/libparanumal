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

#include "acoustics.h"
#include <fstream>
#include <vector>
#include <highfive/H5Easy.hpp>
#include <highfive/H5File.hpp>

using namespace HighFive;

void acousticsWriteXdmfHeader(acoustics_t *acoustics, size_t Nelem, string fileNameXdmf, string filenameH5, 
  string dataTag, std::vector<dfloat> timeVector);
void acousticsWriteMeshH5(string fileName, string fileTag0, string fileTag1, 
  std::vector<dfloat> &x1d, std::vector<dfloat> &y1d, std::vector<dfloat> &z1d, uint fileAttr);

void acousticsWriteXdmf(acoustics_t *acoustics, std::vector<dfloat> timeVector, int iter) {
  string filenameH5 = acoustics->simulationID + ".h5";
  std::string filepathH5 = acoustics->outDir + "/" + filenameH5;  

  auto conn = std::vector<uint>();
  auto x1d = std::vector<dfloat>();
  auto y1d = std::vector<dfloat>();
  auto z1d = std::vector<dfloat>();
  auto p1d = std::vector<dfloat>();

  extractUniquePoints(acoustics->mesh, acoustics, conn, x1d, y1d, z1d, p1d);

  if (iter == 0) {
    // write mesh and Xdmf once
    std::string filepathXdmf = acoustics->outDir + "/" + acoustics->simulationID + ".xdmf";
    std::string filepathUniformXdmf = acoustics->outDir + "/" + acoustics->simulationID + "_uniform.xdmf";
    std::string dataTag = "/data";
    std::string fileTag0 = "/data0";
    std::string fileTag1 = "/data1";

    acousticsWriteXdmfHeader(acoustics, x1d.size(), filepathXdmf, filenameH5, dataTag, timeVector);
    acousticsWriteMeshH5(filepathH5, fileTag0, fileTag1, x1d, y1d, z1d, File::Overwrite);

    H5Easy::File file(filepathH5, File::OpenOrCreate);
    
    if (acoustics->x1d_uniform.size() > 0) {
      dataTag = "/udata";
      fileTag0 = "/udata0";
      fileTag1 = "/udata1";
      std::string fileTag2 = "/udata2";

      auto x1d_u = acoustics->x1d_uniform;
      auto y1d_u = acoustics->y1d_uniform;
      auto z1d_u = acoustics->z1d_uniform;
      acousticsWriteXdmfHeader(acoustics, acoustics->z1d_uniform.size(), filepathUniformXdmf, filenameH5, dataTag, {0.0});
      acousticsWriteMeshH5(filepathH5, fileTag0, fileTag1, x1d_u, y1d_u, z1d_u, File::ReadWrite);
      
      H5Easy::dump(file, "ic_uniform_shape", std::vector<dfloat>({acoustics->ic_uniform_shape[0],
                                                                  acoustics->ic_uniform_shape[1],
                                                                  acoustics->ic_uniform_shape[2]}));
      H5Easy::dump(file, fileTag2, acoustics->ic_uniform);
    }

    if (acoustics->sourceType == GaussianFunction) {
      // write source position      
      DataSet dataset = file.createDataSet<dfloat>("/source_position",  DataSpace::From(acoustics->sourcePosition));
      dataset.write(acoustics->sourcePosition);
    }
  }

  std::string fileTag = "/data" + std::to_string(iter + 2);
  H5Easy::File file(filepathH5, File::OpenOrCreate);
  H5Easy::dump(file, fileTag, p1d);
}

void acousticsWriteXdmfHeader(acoustics_t *acoustics, size_t Nelem, 
  string fileNameXdmf, string filenameH5, string dataTag, std::vector<dfloat> timeVector)
{
  std::ofstream ofs(fileNameXdmf.c_str(), std::ofstream::out);
  // https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/data_into_visit/XdmfFormat.html#an-example-of-a-point-mesh
  
  // HEADER
  ofs << "<?xml version=\"1.0\" ?>" << std::endl;
  ofs << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
  ofs << "<Xdmf Version=\"3.0\">"  << std::endl;

  // Begin
  ofs << "  <Domain>" << std::endl;
  ofs << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;

  for (int i=0; i < timeVector.size(); i++) {
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

void acousticsWriteMeshH5(string fileName, string fileTag0, string fileTag1, 
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

  DataSet dataset = file.createDataSet<dfloat>(fileTag0,  DataSpace::From(meshData));
  dataset.write(meshData);

  dataset = file.createDataSet<int>(fileTag1,  DataSpace::From(vertexData));
  dataset.write(vertexData);
}

// interpolate data to plot nodes and save to file (one per process
void acousticsWriteVTU(acoustics_t *acoustics, bool writeVelocity /*= false */)
{
  mesh_t *mesh = acoustics->mesh;
  
  char fname[BUFSIZ];  
  sprintf(fname, "%s/pressure_field_%04d_%04d.vtu", (char*)acoustics->outDir.c_str(), mesh->rank, acoustics->frame++);
  std::ofstream ofs(fname, std::ofstream::out | std::ofstream::binary);

  ofs << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"BigEndian\">" << std::endl;
  ofs << "  <UnstructuredGrid>" << std::endl;
  ofs << "    <Piece NumberOfPoints=\"" << mesh->Nelements * mesh->plotNp
      << "\" NumberOfCells=\"" << mesh->Nelements * mesh->plotNelements << "\">" << std::endl;
  
  // write out nodes
  ofs << "      <Points>" << std::endl;
  ofs << "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;

  // compute plot node coordinates on the fly
  for (dlong e = 0; e < mesh->Nelements; ++e)
  {
    for (int n = 0; n < mesh->plotNp; ++n)
    {
      dfloat plotxn = 0, plotyn = 0, plotzn = 0;

      for (int m = 0; m < mesh->Np; ++m)
      {
        plotxn += mesh->plotInterp[n * mesh->Np + m] * mesh->x[m + e * mesh->Np];
        plotyn += mesh->plotInterp[n * mesh->Np + m] * mesh->y[m + e * mesh->Np];
        plotzn += mesh->plotInterp[n * mesh->Np + m] * mesh->z[m + e * mesh->Np];
      }

      ofs << "       ";
      ofs << plotxn << " " << plotyn << " " << plotzn;
    }
  }
  ofs << "        </DataArray>" << std::endl;
  ofs << "      </Points>" << std::endl;

  // write out pressure
  ofs << "      <PointData Scalars=\"scalars\">" << std::endl;
  ofs << "        <DataArray type=\"Float64\" Name=\"Pressure\" Format=\"ascii\">" << std::endl;
  for (dlong e = 0; e < mesh->Nelements; ++e)
  {
    for (int n = 0; n < mesh->plotNp; ++n)
    {
      dfloat plotpn = 0;
      for (int m = 0; m < mesh->Np; ++m)
      {
        dfloat pm = acoustics->q[e * mesh->Np * mesh->Nfields + m];
        plotpn += mesh->plotInterp[n * mesh->Np + m] * pm;
      }

      ofs << "       ";
      ofs << plotpn << std::endl;
    }
  }
  ofs << "       </DataArray>" << std::endl;

  if (writeVelocity)
  {
    // write out velocity
    ofs << "        <DataArray type=\"Float64\" Name=\"Velocity\" NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;
    for (dlong e = 0; e < mesh->Nelements; ++e)
    {
      for (int n = 0; n < mesh->plotNp; ++n)
      {
        dfloat plotun = 0, plotvn = 0, plotwn = 0;
        for (int m = 0; m < mesh->Np; ++m)
        {
          dfloat rm = acoustics->q[e * mesh->Np * mesh->Nfields + m];
          dfloat um = acoustics->q[e * mesh->Np * mesh->Nfields + m + mesh->Np];
          dfloat vm = acoustics->q[e * mesh->Np * mesh->Nfields + m + mesh->Np * 2];
          //
          plotun += mesh->plotInterp[n * mesh->Np + m] * um;
          plotvn += mesh->plotInterp[n * mesh->Np + m] * vm;

          if (acoustics->dim == 3)
          {
            dfloat wm = acoustics->q[e * mesh->Np * mesh->Nfields + m + mesh->Np * 3];

            plotwn += mesh->plotInterp[n * mesh->Np + m] * wm;
          }
        }

        ofs << "       ";
        ofs << plotun << " " << plotvn << " " << plotwn;
      }
    }
    ofs << "       </DataArray>" << std::endl;
  }

  ofs << "     </PointData>" << std::endl;


  ofs << "    <Cells>" << std::endl;
  ofs << "      <DataArray type=\"Int32\" Name=\"connectivity\" Format=\"ascii\">" << std::endl;

  for (dlong e = 0; e < mesh->Nelements; ++e)
  {
    for (int n = 0; n < mesh->plotNelements; ++n)
    {
      ofs << "       ";
      for (int m = 0; m < mesh->plotNverts; ++m)
      {
        ofs << e * mesh->plotNp + mesh->plotEToV[n * mesh->plotNverts + m] << " ";
      }
      ofs << std::endl;
    }
  }
  ofs << "        </DataArray>" << std::endl;

  ofs << "        <DataArray type=\"Int32\" Name=\"offsets\" Format=\"ascii\">" << std::endl;
  dlong cnt = 0;
  for (dlong e = 0; e < mesh->Nelements; ++e)
  {
    for (int n = 0; n < mesh->plotNelements; ++n)
    {
      cnt += mesh->plotNverts;
      ofs << "       ";
      ofs << cnt;
    }
  }
  ofs << "       </DataArray>" << std::endl;

  ofs << "       <DataArray type=\"Int32\" Name=\"types\" Format=\"ascii\">" << std::endl;
  for (dlong e = 0; e < mesh->Nelements; ++e)
  {
    for (int n = 0; n < mesh->plotNelements; ++n)
    {
      if (acoustics->dim == 2)
        ofs << "5" << std::endl;
      else
        ofs << "10" << std::endl;
    }
  }
  ofs << "        </DataArray>" << std::endl;
  ofs << "      </Cells>" << std::endl;


  ofs << "    </Piece>" << std::endl;
  ofs << "  </UnstructuredGrid>" << std::endl;
  ofs << "</VTKFile>" << std::endl;

  ofs.close();
}