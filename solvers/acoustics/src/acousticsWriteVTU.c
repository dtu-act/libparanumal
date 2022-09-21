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
#include <highfive/H5Easy.hpp>
#include <highfive/H5File.hpp>

using namespace HighFive;

void acousticsWriteXdmfHeader(acoustics_t *acoustics, string fileNameXdmf, string filenameH5, std::vector<dfloat> timeVector);
void acousticsWriteMeshH5(acoustics_t *acoustics, string fileName, string fileTag0, string fileTag1);
void acousticsWriteH5(acoustics_t *acoustics, std::string fileName, std::string fileTag);

void acousticsWriteXdmf(acoustics_t *acoustics, std::vector<dfloat> timeVector, int iter) {
  string filenameH5 = acoustics->simulationID + ".h5";
  std::string filepathH5 = acoustics->outDir + "/" + filenameH5;  

  if (iter == 0) {
    // write mesh and Xdmf once
    std::string filepathXdmf = acoustics->outDir + "/" + acoustics->simulationID + ".xdmf";    
    std::string fileTag0 = "/data0";
    std::string fileTag1 = "/data1";
    
    acousticsWriteXdmfHeader(acoustics, filepathXdmf, filenameH5, timeVector);
    acousticsWriteMeshH5(acoustics, filepathH5, fileTag0, fileTag1);    
  }

  std::string fileTag = "/data" + std::to_string(iter + 2);
  acousticsWriteH5(acoustics, filepathH5, fileTag);
}

void acousticsWriteXdmfHeader(acoustics_t *acoustics, string fileNameXdmf, string filenameH5, std::vector<dfloat> timeVector)
{
  mesh_t *mesh = acoustics->mesh;  

  std::ofstream ofs(fileNameXdmf.c_str(), std::ofstream::out);

  size_t N_elem = mesh->Nelements * mesh->plotNp;
  // https://visit-sphinx-github-user-manual.readthedocs.io/en/develop/data_into_visit/XdmfFormat.html#an-example-of-a-point-mesh
  
  // HEADER
  ofs << "<?xml version=\"1.0\" ?>" << std::endl;
  ofs << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>" << std::endl;
  ofs << "<Xdmf Version=\"2.0\">"  << std::endl;

  // Begin
  ofs << "  <Domain>" << std::endl;
  ofs << "    <Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">" << std::endl;

  for (int i=0; i < timeVector.size(); i++) {
    std::string fileTag = "/data" + std::to_string(i + 2);

    ofs << "      <Grid>" << std::endl;
    ofs << "        <include xpointer=\"xpointer(//Grid[@Name=&quot;mesh&quot;]/*[self::Topology or self::Geometry])\" />" << std::endl;
    ofs << "        <Time Value=\"" << timeVector[i] << "\" />" << std::endl;
		ofs << "        <Attribute Name=\"p\" AttributeType=\"Scalar\" Center=\"Node\">" << std::endl;
		ofs << "          <DataItem DataType=\"Float\" Dimensions=\"" << N_elem << "\" Format=\"HDF\" Precision=\"8\">" << std::endl;
		ofs << "            " << filenameH5 << ":" << fileTag << std::endl;
		ofs << "          </DataItem>" << std::endl;
    ofs << "        </Attribute>" << std::endl;
    ofs << "      </Grid>" << std::endl;
  }
  ofs << "    </Grid>" << std::endl;

  ofs << "    <Grid Name=\"mesh\" GridType=\"Uniform\">" << std::endl;
  ofs << "      <Geometry GeometryType=\"XYZ\">" << std::endl;
  ofs << "        <DataItem DataType=\"Float\" Dimensions=\"" << N_elem << " " << 3 << "\" Format=\"HDF\" Precision=\"8\">" << std::endl;
  ofs << "          " << filenameH5 << ":/data0" << endl;
  ofs << "        </DataItem>" << std::endl;
  ofs << "      </Geometry>" << std::endl;
  ofs << "      <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << N_elem << "\">" << std::endl;
  ofs << "        <DataItem DataType=\"Int\" Dimensions=\"" << N_elem << " " << 1 << "\" Format=\"HDF\" Precision=\"8\">" << std::endl;
  ofs << "          " << filenameH5 << ":/data1" << endl;
  ofs << "        </DataItem>" << std::endl;
  ofs << "      </Topology>" << std::endl;
  ofs << "    </Grid>" << std::endl;
  ofs << "  </Domain>" << std::endl;
  ofs << "</Xdmf>" << std::endl;
}

void acousticsWriteMeshH5(acoustics_t *acoustics, string fileName, string fileTag0, string fileTag1)
{  
  File file(fileName, File::Overwrite);
  
  auto mesh = acoustics->mesh;

  size_t Nelem = mesh->Nelements * mesh->plotNp;
  
  std::vector<dfloat> rowVec(3);
  std::vector<std::vector<dfloat>> meshData(Nelem, rowVec);
  std::vector<int> vertexData(Nelem);

  int i = 0;
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

      meshData[i][0] = plotxn;
      meshData[i][1] = plotyn;
      meshData[i][2] = plotzn;
      vertexData[i] = i;

      i++; 
    }
  }

  DataSet dataset = file.createDataSet<dfloat>(fileTag0,  DataSpace::From(meshData));
  dataset.write(meshData);

  dataset = file.createDataSet<int>(fileTag1,  DataSpace::From(vertexData));
  dataset.write(vertexData);
}

void acousticsWriteH5(acoustics_t *acoustics, std::string fileName, std::string fileTag)
{
  H5Easy::File file(fileName, File::OpenOrCreate);
  mesh_t *mesh = acoustics->mesh;
  std::vector<dfloat> pressures;

  // write out pressure
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
      pressures.push_back(plotpn);
    }    
  }
                             
  H5Easy::dump(file, fileTag, pressures);

  // if (writeVelocity) {

  //   std::vector vel_x = std::vector();
  //   std::vector vel_y = std::vector();
  //   std::vector vel_z = std::vector();
    
  //   for (dlong e = 0; e < mesh->Nelements; ++e)
  //   {
  //     for (int n = 0; n < mesh->plotNp; ++n)
  //     {
  //       dfloat plotun = 0, plotvn = 0, plotwn = 0;
  //       for (int m = 0; m < mesh->Np; ++m)
  //       {
  //         dfloat rm = acoustics->q[e * mesh->Np * mesh->Nfields + m];
  //         dfloat um = acoustics->q[e * mesh->Np * mesh->Nfields + m + mesh->Np];
  //         dfloat vm = acoustics->q[e * mesh->Np * mesh->Nfields + m + mesh->Np * 2];
  //         //
  //         plotun += mesh->plotInterp[n * mesh->Np + m] * um;
  //         plotvn += mesh->plotInterp[n * mesh->Np + m] * vm;

  //         if (acoustics->dim == 3)
  //         {
  //           dfloat wm = acoustics->q[e * mesh->Np * mesh->Nfields + m + mesh->Np * 3];

  //           plotwn += mesh->plotInterp[n * mesh->Np + m] * wm;
  //         }
  //       }

  //       vel_x.push_back(plotun);
  //       vel_y.push_back(plotvn);
  //       vel_z.push_back(plotwn);
  //     }
  //   }

  //   H5Easy::dump(file, fileTag, plotun);
  //   H5Easy::dump(file, fileTag, plotvn);
  //   H5Easy::dump(file, fileTag, plotwn);
  // }
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