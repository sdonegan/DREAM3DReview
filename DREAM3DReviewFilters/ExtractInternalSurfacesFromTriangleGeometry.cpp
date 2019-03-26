/* ============================================================================
* Copyright (c) 2009-2016 BlueQuartz Software, LLC
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
* Neither the name of BlueQuartz Software, the US Air Force, nor the names of its
* contributors may be used to endorse or promote products derived from this software
* without specific prior written permission.
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
* The code contained herein was partially funded by the followig contracts:
*    United States Air Force Prime Contract FA8650-07-D-5800
*    United States Air Force Prime Contract FA8650-10-D-5210
*    United States Prime Contract Navy N00173-07-C-2068
*
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#include "ExtractInternalSurfacesFromTriangleGeometry.h"

#include <array>
#include <cassert>
#include <unordered_map>

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/DataContainerSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Geometry/TriangleGeom.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ExtractInternalSurfacesFromTriangleGeometry::ExtractInternalSurfacesFromTriangleGeometry()
: m_TriangleDataContainerName(SIMPL::Defaults::TriangleDataContainerName)
, m_NodeTypesArrayPath(SIMPL::Defaults::TriangleDataContainerName, SIMPL::Defaults::VertexAttributeMatrixName, SIMPL::VertexData::SurfaceMeshNodeType)
, m_InternalTrianglesName("InternalTrianglesDataContainer")
{
  initialize();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ExtractInternalSurfacesFromTriangleGeometry::~ExtractInternalSurfacesFromTriangleGeometry()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractInternalSurfacesFromTriangleGeometry::initialize()
{
  setErrorCondition(0);
  setWarningCondition(0);
  setCancel(false);
  m_AttrMatList.clear();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractInternalSurfacesFromTriangleGeometry::setupFilterParameters()
{
  FilterParameterVector parameters;
  DataContainerSelectionFilterParameter::RequirementType dcsReq;
  IGeometry::Types geomTypes = {IGeometry::Type::Triangle};
  parameters.push_back(SIMPL_NEW_DC_SELECTION_FP("Triangle Data Container", TriangleDataContainerName, FilterParameter::RequiredArray, ExtractInternalSurfacesFromTriangleGeometry, dcsReq));
  parameters.push_back(SeparatorFilterParameter::New("Vertex Data", FilterParameter::RequiredArray));
  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int8, 1, AttributeMatrix::Type::Vertex, IGeometry::Type::Triangle);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Node Types", NodeTypesArrayPath, FilterParameter::RequiredArray, ExtractInternalSurfacesFromTriangleGeometry, dasReq));
  parameters.push_back(SIMPL_NEW_STRING_FP("Internal Triangles Data Container", InternalTrianglesName, FilterParameter::CreatedArray, ExtractInternalSurfacesFromTriangleGeometry));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractInternalSurfacesFromTriangleGeometry::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setTriangleDataContainerName(reader->readString("TriangleDataContainerName", getTriangleDataContainerName()));
  setNodeTypesArrayPath(reader->readDataArrayPath("NodeTypesArrayPath", getNodeTypesArrayPath()));
  setInternalTrianglesName(reader->readString("InternalTrianglesName", getInternalTrianglesName()));
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractInternalSurfacesFromTriangleGeometry::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);

  QVector<IDataArray::Pointer> arrays;

  TriangleGeom::Pointer tris = getDataContainerArray()->getPrereqGeometryFromDataContainer<TriangleGeom, AbstractFilter>(this, getTriangleDataContainerName());

  if(getErrorCondition() < 0)
  {
    return;
  }

  arrays.push_back(tris->getVertices());

  QVector<size_t> cDims(1, 1);

  m_NodeTypesPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int8_t>, AbstractFilter>(this, getNodeTypesArrayPath(), cDims);
  if(m_NodeTypesPtr.lock())
  {
    m_NodeTypes = m_NodeTypesPtr.lock()->getPointer(0);
  }
  if(getErrorCondition() >= 0)
  {
    arrays.push_back(m_NodeTypesPtr.lock());
  }

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, arrays);

  DataContainer::Pointer dc = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getInternalTrianglesName());

  if(getErrorCondition() < 0)
  {
    return;
  }

  SharedVertexList::Pointer vertices = TriangleGeom::CreateSharedVertexList(0, !getInPreflight());
  TriangleGeom::Pointer internalTris = TriangleGeom::CreateGeometry(0, vertices, SIMPL::Geometry::TriangleGeometry, !getInPreflight());
  dc->setGeometry(internalTris);

  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getTriangleDataContainerName());
  m_AttrMatList = m->getAttributeMatrixNames();
  QVector<size_t> tDims(1, 0);
  QList<QString> tempDataArrayList;
  DataArrayPath tempPath;
  AttributeMatrix::Type tempAttrMatType = AttributeMatrix::Type::Vertex;

  for(auto&& attr_mat : m_AttrMatList)
  {
    AttributeMatrix::Pointer tmpAttrMat = m->getPrereqAttributeMatrix(this, attr_mat, -301);
    if(getErrorCondition() >= 0)
    {
      tempAttrMatType = tmpAttrMat->getType();
      if(tempAttrMatType == AttributeMatrix::Type::Vertex || tempAttrMatType == AttributeMatrix::Type::Face)
      {
        dc->createNonPrereqAttributeMatrix(this, tmpAttrMat->getName(), tDims, tempAttrMatType);
        tempDataArrayList = tmpAttrMat->getAttributeArrayNames();
        for(auto&& data_array : tempDataArrayList)
        {
          tempPath.update(getInternalTrianglesName(), tmpAttrMat->getName(), data_array);
          IDataArray::Pointer tmpDataArray = tmpAttrMat->getPrereqIDataArray<IDataArray, AbstractFilter>(this, data_array, -90002);
          if(getErrorCondition() >= 0)
          {
            QVector<size_t> cDims = tmpDataArray->getComponentDimensions();
            TemplateHelpers::CreateNonPrereqArrayFromArrayType()(this, tempPath, cDims, tmpDataArray);
          }
        }
      }
      else
      {
        AttributeMatrix::Pointer attrMat = tmpAttrMat->deepCopy(getInPreflight());
        dc->addAttributeMatrix(attr_mat, attrMat);
      }
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractInternalSurfacesFromTriangleGeometry::preflight()
{
  // These are the REQUIRED lines of CODE to make sure the filter behaves correctly
  setInPreflight(true);              // Set the fact that we are preflighting.
  emit preflightAboutToExecute();    // Emit this signal so that other widgets can do one file update
  emit updateFilterParameters(this); // Emit this signal to have the widgets push their values down to the filter
  dataCheck();                       // Run our DataCheck to make sure everthing is setup correctly
  emit preflightExecuted();          // We are done preflighting this filter
  setInPreflight(false);             // Inform the system this filter is NOT in preflight mode anymore.
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <class T> inline void hashCombine(size_t& seed, const T& obj)
{
  std::hash<T> hasher;
  seed ^= hasher(obj) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void copyData(IDataArray::Pointer inDataPtr, IDataArray::Pointer outDataPtr, std::unordered_map<int64_t, int64_t>& elementMap)
{
  typename DataArray<T>::Pointer inputDataPtr = std::dynamic_pointer_cast<DataArray<T>>(inDataPtr);
  T* inputData = static_cast<T*>(inputDataPtr->getPointer(0));
  typename DataArray<T>::Pointer croppedDataPtr = std::dynamic_pointer_cast<DataArray<T>>(outDataPtr);
  T* outputData = static_cast<T*>(croppedDataPtr->getPointer(0));

  size_t nTuples = outDataPtr->getNumberOfTuples();
  size_t nComps = inDataPtr->getNumberOfComponents();
  size_t tmpIndex = 0;
  size_t ptrIndex = 0;

  for(size_t i = 0; i < nTuples; i++)
  {
    for(size_t d = 0; d < nComps; d++)
    {
      tmpIndex = nComps * i + d;
      ptrIndex = nComps * elementMap[i] + d;
      outputData[tmpIndex] = inputData[ptrIndex];
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ExtractInternalSurfacesFromTriangleGeometry::execute()
{
  initialize();
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  TriangleGeom::Pointer tris = getDataContainerArray()->getDataContainer(m_TriangleDataContainerName)->getGeometryAs<TriangleGeom>();
  float* vertices = tris->getVertexPointer(0);
  int64_t* triangles = tris->getTriPointer(0);
  int64_t numVerts = tris->getNumberOfVertices();
  int64_t numTris = tris->getNumberOfTris();

  typedef std::array<float, 3> Vertex;

  struct ArrayHasher
  {
    size_t operator()(const Vertex& vert) const
    {
      size_t hash = std::hash<float>()(vert[0]);
      hashCombine(hash, vert[1]);
      hashCombine(hash, vert[2]);
      return hash;
    }
  };

  typedef std::array<int64_t, 3> Triangle;
  typedef std::unordered_map<Vertex, int64_t, ArrayHasher> VertexMap;
  VertexMap vertexMap;
  std::unordered_map<int64_t, int64_t> internalVertexMap;
  std::unordered_map<int64_t, int64_t> internalTriMap;
  std::vector<Vertex> tmpVerts;
  tmpVerts.reserve(numVerts);
  std::vector<Triangle> tmpTris;
  tmpTris.reserve(numTris);
  std::vector<int8_t> tmpNodeTypes;
  tmpNodeTypes.reserve(numVerts);
  int64_t vertCounter = 0;
  int64_t triCounter = 0;
  int64_t tmpVert0 = 0;
  int64_t tmpVert1 = 0;
  int64_t tmpVert2 = 0;

  int64_t progIncrement = numTris / 100;
  int64_t prog = 1;
  int64_t progressInt = 0;
  int64_t counter = 0;

  for(int64_t i = 0; i < numTris; i++)
  {
    if(getCancel())
    {
      return;
    }
    if((m_NodeTypes[triangles[3 * i + 0]] == 2 || m_NodeTypes[triangles[3 * i + 0]] == 3 || m_NodeTypes[triangles[3 * i + 0]] == 4) &&
       (m_NodeTypes[triangles[3 * i + 1]] == 2 || m_NodeTypes[triangles[3 * i + 1]] == 3 || m_NodeTypes[triangles[3 * i + 1]] == 4) &&
       (m_NodeTypes[triangles[3 * i + 2]] == 2 || m_NodeTypes[triangles[3 * i + 2]] == 3 || m_NodeTypes[triangles[3 * i + 2]] == 4))
    {
      Vertex tmpCoords1 = {{vertices[3 * triangles[3 * i + 0] + 0], vertices[3 * triangles[3 * i + 0] + 1], vertices[3 * triangles[3 * i + 0] + 2]}};
      Vertex tmpCoords2 = {{vertices[3 * triangles[3 * i + 1] + 0], vertices[3 * triangles[3 * i + 1] + 1], vertices[3 * triangles[3 * i + 1] + 2]}};
      Vertex tmpCoords3 = {{vertices[3 * triangles[3 * i + 2] + 0], vertices[3 * triangles[3 * i + 2] + 1], vertices[3 * triangles[3 * i + 2] + 2]}};

      auto iter = vertexMap.find(tmpCoords1);
      if(iter == vertexMap.end())
      {
        tmpVerts.push_back(tmpCoords1);
        tmpNodeTypes.push_back(m_NodeTypes[triangles[3 * i + 0]]);
        vertexMap[tmpCoords1] = vertCounter;
        tmpVert0 = vertCounter;
        internalVertexMap[tmpVert0] = triangles[3 * i + 0];
        vertCounter++;
      }
      else
      {
        tmpVert0 = vertexMap[tmpCoords1];
      }

      iter = vertexMap.find(tmpCoords2);
      if(iter == vertexMap.end())
      {
        tmpVerts.push_back(tmpCoords2);
        tmpNodeTypes.push_back(m_NodeTypes[triangles[3 * i + 1]]);
        vertexMap[tmpCoords2] = vertCounter;
        tmpVert1 = vertCounter;
        internalVertexMap[tmpVert1] = triangles[3 * i + 1];
        vertCounter++;
      }
      else
      {
        tmpVert1 = vertexMap[tmpCoords2];
      }

      iter = vertexMap.find(tmpCoords3);
      if(iter == vertexMap.end())
      {
        tmpVerts.push_back(tmpCoords3);
        tmpNodeTypes.push_back(m_NodeTypes[triangles[3 * i + 2]]);
        vertexMap[tmpCoords3] = vertCounter;
        tmpVert2 = vertCounter;
        internalVertexMap[tmpVert2] = triangles[3 * i + 2];
        vertCounter++;
      }
      else
      {
        tmpVert2 = vertexMap[tmpCoords3];
      }

      Triangle tmpTri = { { tmpVert0, tmpVert1, tmpVert2} };
      tmpTris.push_back(tmpTri);
      internalTriMap[triCounter] = i;
      triCounter++;
    }

    if(counter > prog)
    {
      progressInt = static_cast<int64_t>((static_cast<float>(counter) / numTris) * 100.0f);
      QString ss = QObject::tr("Checking Triangle %1 of %2 || %3% Completed").arg(counter).arg(numTris).arg(progressInt);
      notifyStatusMessage(getMessagePrefix(), ss);
      prog = prog + progIncrement;
    }
    counter++;
  }

  QString ss = QObject::tr("Finished Checking Triangles || Updating Array Information...");
  notifyStatusMessage(getMessagePrefix(), ss);

  tmpVerts.shrink_to_fit();
  tmpTris.shrink_to_fit();

  QVector<size_t> vertDims(1, tmpVerts.size());
  QVector<size_t> triDims(1, tmpTris.size());
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(m_TriangleDataContainerName);
  DataContainer::Pointer dc = getDataContainerArray()->getDataContainer(m_InternalTrianglesName);

  for(auto&& attr_mat : m_AttrMatList)
  {
    AttributeMatrix::Pointer tmpAttrMat = dc->getPrereqAttributeMatrix(this, attr_mat, -301);
    if(getErrorCondition() >= 0)
    {
      AttributeMatrix::Type tempAttrMatType = tmpAttrMat->getType();
      if(tempAttrMatType == AttributeMatrix::Type::Vertex)
      {
        tmpAttrMat->resizeAttributeArrays(vertDims);
        QList<QString> srcDataArrays = tmpAttrMat->getAttributeArrayNames();
        AttributeMatrix::Pointer srcAttrMat = m->getAttributeMatrix(tmpAttrMat->getName());
        assert(srcAttrMat);

        for(auto&& data_array : srcDataArrays)
        {
          IDataArray::Pointer src = srcAttrMat->getAttributeArray(data_array);
          IDataArray::Pointer dest = tmpAttrMat->getAttributeArray(data_array);

          assert(src);
          assert(dest);
          assert(src->getNumberOfComponents() == dest->getNumberOfComponents());

          EXECUTE_FUNCTION_TEMPLATE(this, copyData, src, src, dest, internalVertexMap)
        }
      }
      else if(tempAttrMatType == AttributeMatrix::Type::Face)
      {
        tmpAttrMat->resizeAttributeArrays(triDims);
        QList<QString> srcDataArrays = tmpAttrMat->getAttributeArrayNames();
        AttributeMatrix::Pointer srcAttrMat = m->getAttributeMatrix(tmpAttrMat->getName());
        assert(srcAttrMat);

        for(auto&& data_array : srcDataArrays)
        {
          IDataArray::Pointer src = srcAttrMat->getAttributeArray(data_array);
          IDataArray::Pointer dest = tmpAttrMat->getAttributeArray(data_array);

          assert(src);
          assert(dest);
          assert(src->getNumberOfComponents() == dest->getNumberOfComponents());

          EXECUTE_FUNCTION_TEMPLATE(this, copyData, src, src, dest, internalTriMap)
        }
      }
    }
  }

  TriangleGeom::Pointer internalTris = getDataContainerArray()->getDataContainer(m_InternalTrianglesName)->getGeometryAs<TriangleGeom>();
  internalTris->resizeVertexList(tmpVerts.size());
  internalTris->resizeTriList(tmpTris.size());
  float* internalVertices = internalTris->getVertexPointer(0);
  int64_t* internalTriangles = internalTris->getTriPointer(0);
  int64_t numInternalTris = internalTris->getNumberOfTris();
  int64_t numInternalVerts = internalTris->getNumberOfVertices();

  for(int64_t i = 0; i < numInternalVerts; i++)
  {
    internalVertices[3 * i + 0] = tmpVerts[i][0];
    internalVertices[3 * i + 1] = tmpVerts[i][1];
    internalVertices[3 * i + 2] = tmpVerts[i][2];
  }

  for(int64_t i = 0; i < numInternalTris; i++)
  {
    internalTriangles[3 * i + 0] = tmpTris[i][0];
    internalTriangles[3 * i + 1] = tmpTris[i][1];
    internalTriangles[3 * i + 2] = tmpTris[i][2];
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ExtractInternalSurfacesFromTriangleGeometry::newFilterInstance(bool copyFilterParameters) const
{
  ExtractInternalSurfacesFromTriangleGeometry::Pointer filter = ExtractInternalSurfacesFromTriangleGeometry::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractInternalSurfacesFromTriangleGeometry::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractInternalSurfacesFromTriangleGeometry::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractInternalSurfacesFromTriangleGeometry::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractInternalSurfacesFromTriangleGeometry::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid ExtractInternalSurfacesFromTriangleGeometry::getUuid()
{
  return QUuid("{52a069b4-6a46-5810-b0ec-e0693c636034}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractInternalSurfacesFromTriangleGeometry::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::GeometryFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ExtractInternalSurfacesFromTriangleGeometry::getHumanLabel() const
{
  return "Extract Internal Surfaces from Triangle Geometry";
}
