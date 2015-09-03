/* ============================================================================
 * Copyright (c) 2011 Michael A. Jackson (BlueQuartz Software)
 * Copyright (c) 2011 Dr. Michael A. Groeber (US Air Force Research Laboratories)
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
 * Neither the name of Michael A. Groeber, Michael A. Jackson, the US Air Force,
 * BlueQuartz Software nor the names of its contributors may be used to endorse
 * or promote products derived from this software without specific prior written
 * permission.
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
 *  This code was written under United States Air Force Contract number
 *                           FA8650-07-D-5800
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#include "ParaDisReader.h"

#include <QtCore/QtDebug>
#include <fstream>

#include <QtCore/QFileInfo>

#include "SIMPLib/Math/MatrixMath.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersWriter.h"

#include "SIMPLib/FilterParameters/InputFileFilterParameter.h"
#include "SIMPLib/FilterParameters/DoubleFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"

#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ParaDisReader::ParaDisReader() :
  FileReader(),
  m_EdgeDataContainerName(DREAM3D::Defaults::DataContainerName),
  m_VertexAttributeMatrixName(DREAM3D::Defaults::VertexAttributeMatrixName),
  m_EdgeAttributeMatrixName(DREAM3D::Defaults::EdgeAttributeMatrixName),
  m_InputFile(""),
  m_BurgersVector(2.5),
  m_NumberOfArmsArrayName(DREAM3D::VertexData::NumberOfArms),
  m_NodeConstraintsArrayName(DREAM3D::VertexData::NodeConstraints),
  m_BurgersVectorsArrayName(DREAM3D::EdgeData::BurgersVectors),
  m_SlipPlaneNormalsArrayName(DREAM3D::EdgeData::SlipPlaneNormals),
  m_DomainBoundsArrayName("DomainBounds"),
  m_NumberOfArms(NULL),
  m_NodeConstraints(NULL),
  m_BurgersVectors(NULL),
  m_SlipPlaneNormals(NULL)
{
  setupFilterParameters();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
ParaDisReader::~ParaDisReader()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ParaDisReader::setupFilterParameters()
{
  FilterParameterVector parameters;
  parameters.push_back(InputFileFilterParameter::New("Input File", "InputFile", getInputFile(), FilterParameter::Parameter, "*"));
  parameters.push_back(DoubleFilterParameter::New("Burgers Vector Length (Angstroms)", "BurgersVector", getBurgersVector(), FilterParameter::Parameter));
// parameters.push_back(SeparatorFilterParameter::New("Created Information", FilterParameter::Uncategorized));
  parameters.push_back(StringFilterParameter::New("Edge DataContainer Name", "EdgeDataContainerName", getEdgeDataContainerName(), FilterParameter::CreatedArray));
  parameters.push_back(StringFilterParameter::New("Vertex AttributeMatrix Name", "VertexAttributeMatrixName", getVertexAttributeMatrixName(), FilterParameter::CreatedArray));
  parameters.push_back(StringFilterParameter::New("Edge AttributeMatrix Name", "EdgeAttributeMatrixName", getEdgeAttributeMatrixName(), FilterParameter::CreatedArray));
  parameters.push_back(StringFilterParameter::New("Number Of Arms Array Name", "NumberOfArmsArrayName", getNumberOfArmsArrayName(), FilterParameter::CreatedArray));
  parameters.push_back(StringFilterParameter::New("Node Constraints Array Name", "NodeConstraintsArrayName", getNodeConstraintsArrayName(), FilterParameter::CreatedArray));
  parameters.push_back(StringFilterParameter::New("Burgers Vectors Array Name", "BurgersVectorsArrayName", getBurgersVectorsArrayName(), FilterParameter::CreatedArray));
  parameters.push_back(StringFilterParameter::New("Slip Plane Normals Array Name", "SlipPlaneNormalsArrayName", getSlipPlaneNormalsArrayName(), FilterParameter::CreatedArray));
  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
void ParaDisReader::readFilterParameters(AbstractFilterParametersReader* reader, int index)
{
  reader->openFilterGroup(this, index);
  setEdgeDataContainerName(reader->readString("EdgeDataContainerName", getEdgeDataContainerName() ) );
  setVertexAttributeMatrixName(reader->readString("VertexAttributeMatrixName", getVertexAttributeMatrixName() ) );
  setEdgeAttributeMatrixName(reader->readString("EdgeAttributeMatrixName", getEdgeAttributeMatrixName() ) );
  setSlipPlaneNormalsArrayName(reader->readString("SlipPlaneNormalsArrayName", getSlipPlaneNormalsArrayName() ) );
  setBurgersVectorsArrayName(reader->readString("BurgersVectorsArrayName", getBurgersVectorsArrayName() ) );
  setNodeConstraintsArrayName(reader->readString("NodeConstraintsArrayName", getNodeConstraintsArrayName() ) );
  setNumberOfArmsArrayName(reader->readString("NumberOfArmsArrayName", getNumberOfArmsArrayName() ) );
  setInputFile( reader->readString( "InputFile", getInputFile() ) );
  setBurgersVector( reader->readValue( "BurgersVector", getBurgersVector() ) );
  reader->closeFilterGroup();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int ParaDisReader::writeFilterParameters(AbstractFilterParametersWriter* writer, int index)
{
  writer->openFilterGroup(this, index);
  SIMPL_FILTER_WRITE_PARAMETER(EdgeDataContainerName)
  SIMPL_FILTER_WRITE_PARAMETER(VertexAttributeMatrixName)
  SIMPL_FILTER_WRITE_PARAMETER(EdgeAttributeMatrixName)
  SIMPL_FILTER_WRITE_PARAMETER(SlipPlaneNormalsArrayName)
  SIMPL_FILTER_WRITE_PARAMETER(BurgersVectorsArrayName)
  SIMPL_FILTER_WRITE_PARAMETER(NodeConstraintsArrayName)
  SIMPL_FILTER_WRITE_PARAMETER(NumberOfArmsArrayName)
  SIMPL_FILTER_WRITE_PARAMETER(InputFile)
  SIMPL_FILTER_WRITE_PARAMETER(BurgersVector)
  writer->closeFilterGroup();
  return ++index; // we want to return the next index that was just written to
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ParaDisReader::updateVertexInstancePointers()
{
  setErrorCondition(0);

  if( NULL != m_NumberOfArmsPtr.lock().get() ) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_NumberOfArms = m_NumberOfArmsPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */
  if( NULL != m_NodeConstraintsPtr.lock().get() ) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_NodeConstraints = m_NodeConstraintsPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ParaDisReader::updateEdgeInstancePointers()
{
  setErrorCondition(0);

  if( NULL != m_BurgersVectorsPtr.lock().get() ) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_BurgersVectors = m_BurgersVectorsPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */
  if( NULL != m_SlipPlaneNormalsPtr.lock().get() ) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_SlipPlaneNormals = m_SlipPlaneNormalsPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ParaDisReader::dataCheck()
{
  DataArrayPath tempPath;
  setErrorCondition(0);
  DataContainer::Pointer m = getDataContainerArray()->createNonPrereqDataContainer<AbstractFilter>(this, getEdgeDataContainerName());
  if(getErrorCondition() < 0) { return; }
  QVector<size_t> tDims(1, 0);
  AttributeMatrix::Pointer amV = m->createNonPrereqAttributeMatrix<AbstractFilter>(this, getVertexAttributeMatrixName(), tDims, DREAM3D::AttributeMatrixType::Vertex);
  if(getErrorCondition() < 0) { return; }
  AttributeMatrix::Pointer amE = m->createNonPrereqAttributeMatrix<AbstractFilter>(this, getEdgeAttributeMatrixName(), tDims, DREAM3D::AttributeMatrixType::Edge);
  if (getErrorCondition() < 0) { return; }
  tDims[0] = 1;
  AttributeMatrix::Pointer amMeta = m->createNonPrereqAttributeMatrix<AbstractFilter>(this, "_MetaData", tDims, DREAM3D::AttributeMatrixType::MetaData);
  if (getErrorCondition() < 0) { return; }

  QFileInfo fi(getInputFile());

  if (getInputFile().isEmpty() == true)
  {
    QString ss = QObject::tr("%1 needs the Input File Set and it was not.").arg(ClassName());
    setErrorCondition(-387);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }
  else if (fi.exists() == false)
  {
    QString ss = QObject::tr("The input file does not exist.");
    setErrorCondition(-388);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }
  QVector<size_t> dims(1, 1);
  tempPath.update(getEdgeDataContainerName(), getVertexAttributeMatrixName(), getNumberOfArmsArrayName() );
  m_NumberOfArmsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter, int32_t>(this,  tempPath, 0, dims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if( NULL != m_NumberOfArmsPtr.lock().get() ) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_NumberOfArms = m_NumberOfArmsPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */
  tempPath.update(getEdgeDataContainerName(), getVertexAttributeMatrixName(), getNodeConstraintsArrayName() );
  m_NodeConstraintsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter, int32_t>(this,  tempPath, 0, dims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if( NULL != m_NodeConstraintsPtr.lock().get() ) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_NodeConstraints = m_NodeConstraintsPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */
  dims[0] = 3;
  tempPath.update(getEdgeDataContainerName(), getEdgeAttributeMatrixName(), getBurgersVectorsArrayName() );
  m_BurgersVectorsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this,  tempPath, 0.0, dims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if( NULL != m_BurgersVectorsPtr.lock().get() ) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_BurgersVectors = m_BurgersVectorsPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */
  tempPath.update(getEdgeDataContainerName(), getEdgeAttributeMatrixName(), getSlipPlaneNormalsArrayName() );
  m_SlipPlaneNormalsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this,  tempPath, 0.0, dims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if( NULL != m_SlipPlaneNormalsPtr.lock().get() ) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_SlipPlaneNormals = m_SlipPlaneNormalsPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */
  dims[0] = 6;
  tempPath.update(getEdgeDataContainerName(), "_MetaData", getDomainBoundsArrayName());
  m_DomainBoundsPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, tempPath, 0.0, dims); /* Assigns the shared_ptr<> to an instance variable that is a weak_ptr<> */
  if (NULL != m_DomainBoundsPtr.lock().get()) /* Validate the Weak Pointer wraps a non-NULL pointer to a DataArray<T> object */
  { m_DomainBounds = m_DomainBoundsPtr.lock()->getPointer(0); } /* Now assign the raw pointer to data from the DataArray<T> object */


  if (m_InStream.isOpen() == true)
  {
    m_InStream.close();
  }

  if (getInputFile().isEmpty() == false && fi.exists() == true)
  {
    // We need to read the header of the input file to get the dimensions
    m_InStream.setFileName(getInputFile());
    if (!m_InStream.open(QIODevice::ReadOnly | QIODevice::Text))
    {
      QString ss = QObject::tr("ParaDisReader Input file could not be opened: %1").arg(getInputFile());
      setErrorCondition(-100);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      return;
    }

    //int error = readHeader();
    //add edges for preflight sake...they will get overwritten when actually reading the file
    SharedVertexList::Pointer vertices = EdgeGeom::CreateSharedVertexList(0);
    EdgeGeom::Pointer edgeGeom = EdgeGeom::CreateGeometry(0, vertices, DREAM3D::Geometry::EdgeGeometry);
    m->setGeometry(edgeGeom);

    //m_InStream.close();
    //if (error < 0)
    //{
    //  setErrorCondition(error);
    //  QString ss = QObject::tr("Error occurred trying to parse the dimensions from the input file. Is the input file a Dx file?");
    //  notifyErrorMessage(getHumanLabel(), ss, -11000);
    //}
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ParaDisReader::preflight()
{
  setInPreflight(true);
  emit preflightAboutToExecute();
  emit updateFilterParameters(this);
  dataCheck();
  emit preflightExecuted();
  setInPreflight(false);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void ParaDisReader::execute()
{
  int err = 0;

  dataCheck();
  if(getErrorCondition() < 0) { return; }

  m_InStream.setFileName(getInputFile());
  if (!m_InStream.open(QIODevice::ReadOnly | QIODevice::Text))
  {
    QString ss = QObject::tr("ParaDisReader Input file could not be opened: %1").arg(getInputFile());
    setErrorCondition(-100);
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }

  err = readHeader();
  if(err < 0)
  {
    m_InStream.close();
    return;
  }
  err = readFile();
  m_InStream.close();
  if(err < 0)
  {
    return;
  }

}


//-----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int ParaDisReader::readHeader()
{
  QString ss;
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getEdgeDataContainerName());
  AttributeMatrix::Pointer vertexAttrMat = m->getAttributeMatrix(getVertexAttributeMatrixName());

  //convert user input Burgers Vector to microns from angstroms
  float burgersVec = m_BurgersVector / 10000.0f;

  int error = 0;

  QByteArray buf;
  QList<QByteArray> tokens; /* vector to store the split data */

  bool ok = false;
  // Process the header information
  //bool done = false;
  //read Version line
  buf = m_InStream.readLine();
  buf = buf.trimmed();
  buf = buf.simplified();
  tokens = buf.split(' ');
  fileVersion = tokens[2].toInt(&ok, 10);

  int keepgoing = 1;
  //read until get to minCoordinates line
  while (keepgoing == 1)
  {
    buf = m_InStream.readLine();
    buf = buf.trimmed();
    buf = buf.simplified();
    tokens = buf.split(' ');
    QString word(tokens.at(0));
    if (word.compare("minCoordinates") == 0)
    {
      for (int32_t i = 0; i < 3; i++)
      {
        buf = m_InStream.readLine();
        buf = buf.trimmed();
        buf = buf.simplified();
        tokens = buf.split(' ');
        m_DomainBounds[i] = tokens[0].toFloat(&ok) * burgersVec;
      }
      keepgoing = 0;
    }
  }
  keepgoing = 1;
  //read until get to maxCoordinates line
  while (keepgoing == 1)
  {
    buf = m_InStream.readLine();
    buf = buf.trimmed();
    buf = buf.simplified();
    tokens = buf.split(' ');
    QString word(tokens.at(0));
    if (word.compare("maxCoordinates") == 0)
    {
      for (int32_t i = 3; i < 6; i++)
      {
        buf = m_InStream.readLine();
        buf = buf.trimmed();
        buf = buf.simplified();
        tokens = buf.split(' ');
        m_DomainBounds[i] = tokens[0].toFloat(&ok) * burgersVec;
      }
      keepgoing = 0;
    }
  }
  keepgoing = 1;
  //read until get to nodeCount line
  while (keepgoing == 1)
  {
    buf = m_InStream.readLine();
    buf = buf.trimmed();
    buf = buf.simplified();
    tokens = buf.split(' ');
    QString word(tokens.at(0));
    if (word.compare("nodeCount") == 0)
    {
      numVerts = tokens[2].toInt(&ok, 10);
      keepgoing = 0;
    }
  }
  keepgoing = 1;
  //read until get to nodalData line
  while(keepgoing == 1)
  {
    buf = m_InStream.readLine();
    buf = buf.trimmed();
    buf = buf.simplified();
    tokens = buf.split(' ');
    QString word(tokens.at(0));
    if(word.compare("nodalData") == 0)
    {
      //read nodal Data lines
      buf = m_InStream.readLine();
      buf = m_InStream.readLine();
      keepgoing = 0;
    }
  }

  EdgeGeom::Pointer edgeGeom = m->getGeometryAs<EdgeGeom>();
  edgeGeom->resizeVertexList(numVerts);

  QVector<size_t> tDims(1, numVerts);
  vertexAttrMat->resizeAttributeArrays(tDims);
  updateVertexInstancePointers();

  return error;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int ParaDisReader::readFile()
{
  DataContainer::Pointer m = getDataContainerArray()->getDataContainer(getEdgeDataContainerName());
  AttributeMatrix::Pointer edgeAttrMat = m->getAttributeMatrix(getEdgeAttributeMatrixName());

  QByteArray buf;
  QList<QByteArray> tokens; /* vector to store the split data */
  QList<QByteArray> subTokens; /* vector to store the split data */

  EdgeGeom::Pointer edgeGeom = m->getGeometryAs<EdgeGeom>();
  edgeGeom->resizeVertexList(numVerts);

  float* vertex = edgeGeom->getVertexPointer(0);
  numVerts = edgeGeom->getNumberOfVertices();

  bool ok = false;

  int nodeNum = 0;
  int neighborNode = 0;
  numEdges = 0;
  int nodeCounter = 0;
  uint64_t* ptr64;

  QMap<int64_t, int32_t> vertNumbers;
  struct
  {
    int32_t n1;
    int32_t n2;
  } uniqueNodeId;

  QVector<int> firstNodes;
  QVector<int> secondNodes;
  QVector<float> burgerXs;
  QVector<float> burgerYs;
  QVector<float> burgerZs;
  float burgVec[3];
  QVector<float> spnXs;
  QVector<float> spnYs;
  QVector<float> spnZs;
  float spNorm[3];

  //convert user input Burgers Vector to microns from angstroms
  float burgersVec = m_BurgersVector / 10000.0f;

  for(int j = 0; j < numVerts; j++)
  {
    buf = m_InStream.readLine();
    buf = buf.trimmed();
    buf = buf.simplified();
    tokens = buf.split(' ');
    subTokens = tokens[0].split(',');
    uniqueNodeId.n1 = subTokens[0].toInt(&ok, 10);
    uniqueNodeId.n2 = subTokens[1].toInt(&ok, 10);
    ptr64 = reinterpret_cast<uint64_t*>(&uniqueNodeId);
    nodeNum = vertNumbers.value(*ptr64, -1);
    if(nodeNum == -1)
    {
      nodeNum = nodeCounter;
      vertex[3 * nodeNum + 0] = tokens[1].toFloat(&ok) * burgersVec;
      vertex[3 * nodeNum + 1] = tokens[2].toFloat(&ok) * burgersVec;
      vertex[3 * nodeNum + 2] = tokens[3].toFloat(&ok) * burgersVec;
      m_NumberOfArms[nodeNum] = tokens[4].toInt(&ok, 10);
      m_NodeConstraints[nodeNum] = tokens[5].toInt(&ok, 10);
      vertNumbers.insert(*ptr64, nodeNum);
      nodeCounter++;
    }
    else
    {
      vertex[3 * nodeNum + 0] = tokens[1].toFloat(&ok) * burgersVec;
      vertex[3 * nodeNum + 1] = tokens[2].toFloat(&ok) * burgersVec;
      vertex[3 * nodeNum + 2] = tokens[3].toFloat(&ok) * burgersVec;
      m_NumberOfArms[nodeNum] = tokens[4].toInt(&ok, 10);
      m_NodeConstraints[nodeNum] = tokens[5].toInt(&ok, 10);
    }
    if(fileVersion >= 5)
    {
      buf = m_InStream.readLine();
    }
    for(int k = 0; k < m_NumberOfArms[nodeNum]; k++)
    {
      buf = m_InStream.readLine();
      buf = buf.trimmed();
      buf = buf.simplified();
      tokens = buf.split(' ');
      subTokens = tokens[0].split(',');
      uniqueNodeId.n1 = subTokens[0].toInt(&ok, 10);
      uniqueNodeId.n2 = subTokens[1].toInt(&ok, 10);
      ptr64 = reinterpret_cast<uint64_t*>(&uniqueNodeId);
      neighborNode = vertNumbers.value(*ptr64, -1);
      if(neighborNode == -1)
      {
        neighborNode = nodeCounter;
        vertNumbers.insert(*ptr64, neighborNode);
        nodeCounter++;
      }
      if(neighborNode > nodeNum)
      {
        numEdges++;
        firstNodes.push_back(nodeNum);
        secondNodes.push_back(neighborNode);
        burgVec[0] = tokens[1].toFloat(&ok);
        burgVec[1] = tokens[2].toFloat(&ok);
        burgVec[2] = tokens[3].toFloat(&ok);
        //burgVec[0] = tokens[1].toFloat(&ok) * burgersVec;
        //burgVec[1] = tokens[2].toFloat(&ok) * burgersVec;
        //burgVec[2] = tokens[3].toFloat(&ok) * burgersVec;
        burgerXs.push_back(burgVec[0]);
        burgerYs.push_back(burgVec[1]);
        burgerZs.push_back(burgVec[2]);
      }
      buf = m_InStream.readLine();
      buf = buf.trimmed();
      buf = buf.simplified();
      tokens = buf.split(' ');
      if(neighborNode > nodeNum)
      {
        spNorm[0] = tokens[0].toFloat(&ok);
        spNorm[1] = tokens[1].toFloat(&ok);
        spNorm[2] = tokens[2].toFloat(&ok);
        MatrixMath::Normalize3x1(spNorm);
        spnXs.push_back(spNorm[0]);
        spnYs.push_back(spNorm[1]);
        spnZs.push_back(spNorm[2]);
      }
    }
  }

  edgeGeom->resizeEdgeList(numEdges);
  int64_t* edge = edgeGeom->getEdgePointer(0);

  // Resize the edge attribute matrix to the number of vertices
  QVector<size_t> tDims (1, numEdges);
  edgeAttrMat->resizeAttributeArrays(tDims);
  updateEdgeInstancePointers();

  for(int i = 0; i < numEdges; i++)
  {
    edge[2 * i + 0] = firstNodes[i];
    edge[2 * i + 1] = secondNodes[i];
    m_BurgersVectors[3 * i + 0] = burgerXs[i];
    m_BurgersVectors[3 * i + 1] = burgerYs[i];
    m_BurgersVectors[3 * i + 2] = burgerZs[i];
    m_SlipPlaneNormals[3 * i + 0] = spnXs[i];
    m_SlipPlaneNormals[3 * i + 1] = spnYs[i];
    m_SlipPlaneNormals[3 * i + 2] = spnZs[i];
  }

  tokens.clear();
  m_InStream.close();

  notifyStatusMessage(getHumanLabel(), "Complete");
  return 0;
}


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer ParaDisReader::newFilterInstance(bool copyFilterParameters)
{
  ParaDisReader::Pointer filter = ParaDisReader::New();
  if(true == copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ParaDisReader::getCompiledLibraryName()
{ return DDDAnalysisToolboxConstants::DDDAnalysisToolboxBaseName; }


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ParaDisReader::getGroupName()
{ return DREAM3D::FilterGroups::Unsupported; }


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ParaDisReader::getSubGroupName()
{ return DREAM3D::FilterSubGroups::InputFilters; }


// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString ParaDisReader::getHumanLabel()
{ return "Read ParaDis File"; }

