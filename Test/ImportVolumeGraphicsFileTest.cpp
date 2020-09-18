/* ============================================================================
 * Copyright (c) 2020 BlueQuartz Software, LLC
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice, this
 * list of conditions and the following disclaimer in the documentation and/or
 * other materials provided with the distribution.
 *
 * Neither the names of any of the BlueQuartz Software contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include <cstdio>
#include <string>

#include <QtCore/QDebug>
#include <QtCore/QDir>
#include <QtCore/QJsonDocument>

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/SIMPLArray.hpp"
#include "SIMPLib/DataArrays/DataArray.hpp"
#include "SIMPLib/DataContainers/AttributeMatrix.h"
#include "SIMPLib/DataContainers/DataContainer.h"
#include "SIMPLib/DataContainers/DataContainerArray.h"

#include "DREAM3DReview/DREAM3DReviewFilters/ImportVolumeGraphicsFile.h"
#include "UnitTestSupport.hpp"

#include "DREAM3DReviewTestFileLocations.h"

class ImportVolumeGraphicsFileTest
{
public:
  ImportVolumeGraphicsFileTest() = default;
  ~ImportVolumeGraphicsFileTest() = default;

  const QString k_VgiFile = UnitTest::PluginSourceDir + "/Test/TestFiles/VolumeGraphicsTest.vgi";

  const QString k_VolFile = UnitTest::TestTempDir + "/VolumeGraphicsTest.vol";
  const SizeVec3Type k_Dimensions = {50, 20, 80};

  // -----------------------------------------------------------------------------
  void RemoveTestFiles()
  {
#if REMOVE_TEST_FILES
    std::string output = QString(UnitTest::TestTempDir + "/VolumeGraphicsTest.vgi").toStdString();
    QFile::remove(QString::fromStdString(output));
    QFile::remove(k_VolFile);
#endif
  }

  // -----------------------------------------------------------------------------
  void PrepareFiles()
  {
    // Write out the data file
    std::string volFile = k_VolFile.toStdString();
    FILE* f = fopen(volFile.c_str(), "wb");
    size_t count = k_Dimensions[0] * k_Dimensions[1] * k_Dimensions[2];
    std::vector<float> data(count, 2.0F);
    if(fwrite(data.data(), sizeof(float), count, f) != count)
    {
      DREAM3D_REQUIRE_EQUAL(1, 0)
    }
    fclose(f);

    // Copy the vgi file
    try
    {
      std::string output = QString(UnitTest::TestTempDir + "/VolumeGraphicsTest.vgi").toStdString();
      fs::copy(k_VgiFile.toStdString(), output); // copy vgi file
    } catch(std::exception& e)
    {
    }
  }

  // -----------------------------------------------------------------------------
  void TestFilter()
  {
    ImportVolumeGraphicsFile::Pointer filter = ImportVolumeGraphicsFile::New();
    {
      DataContainerArray::Pointer dca = DataContainerArray::New();
      filter->setDataContainerArray(dca);
      filter->preflight();
      int32_t err = filter->getErrorCode();
      DREAM3D_REQUIRED(err, ==, -91500)
    }
    {
      DataContainerArray::Pointer dca = DataContainerArray::New();
      filter->setDataContainerArray(dca);
      filter->setVGHeaderFile(QString(UnitTest::TestTempDir));
      filter->preflight();
      int32_t err = filter->getErrorCode();
      DREAM3D_REQUIRED(err, ==, -91507)
    }
    {
      DataContainerArray::Pointer dca = DataContainerArray::New();
      filter->setDataContainerArray(dca);
      filter->setVGHeaderFile(QString(UnitTest::TestTempDir + "/VolumeGraphicsTest.vgi"));
      filter->preflight();
      int32_t err = filter->getErrorCode();
      DREAM3D_REQUIRED(err, ==, 0)
    }
    {
      DataContainerArray::Pointer dca = DataContainerArray::New();
      filter->setDataContainerArray(dca);
      filter->setVGHeaderFile(QString(UnitTest::TestTempDir + "/VolumeGraphicsTest.vgi"));
      filter->execute();
      int32_t err = filter->getErrorCode();
      DREAM3D_REQUIRED(err, ==, 0)

      DataContainer::Pointer dc = dca->getDataContainer("VolumeGraphics");
      ImageGeom::Pointer imageGeom = dc->getGeometryAs<ImageGeom>();
      SizeVec3Type dims = imageGeom->getDimensions();
      DREAM3D_REQUIRED(dims[0], ==, 50)
      DREAM3D_REQUIRED(dims[1], ==, 20)
      DREAM3D_REQUIRED(dims[2], ==, 80)

      AttributeMatrix::Pointer am = dc->getAttributeMatrix("CT Data");
      DREAM3D_REQUIRE_VALID_POINTER(am.get())

      FloatArrayType& data = *(am->getAttributeArrayAs<FloatArrayType>("Density"));

      size_t numTuples = data.getNumberOfTuples();
      for(size_t i = 0; i < numTuples; i++)
      {
        DREAM3D_REQUIRED(data[i], ==, 2.0F)
      }
    }
  }

  // -----------------------------------------------------------------------------
  void operator()()
  {
    std::cout << "###### ImportVolumeGraphicsFileTest ######" << std::endl;
    int err = EXIT_SUCCESS;

    DREAM3D_REGISTER_TEST(PrepareFiles())
    DREAM3D_REGISTER_TEST(TestFilter())
    DREAM3D_REGISTER_TEST(RemoveTestFiles())
  }

public:
  ImportVolumeGraphicsFileTest(const ImportVolumeGraphicsFileTest&) = delete;            // Copy Constructor Not Implemented
  ImportVolumeGraphicsFileTest(ImportVolumeGraphicsFileTest&&) = delete;                 // Move Constructor Not Implemented
  ImportVolumeGraphicsFileTest& operator=(const ImportVolumeGraphicsFileTest&) = delete; // Copy Assignment Not Implemented
  ImportVolumeGraphicsFileTest& operator=(ImportVolumeGraphicsFileTest&&) = delete;      // Move Assignment Not Implemented
};
