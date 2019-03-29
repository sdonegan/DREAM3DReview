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

#include "FindArrayStatistics.h"

#include <cmath>
#include <cstring>
#include <functional>
#include <numeric>
#include <unordered_map>

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#include <tbb/task_scheduler_init.h>
#endif

#include "SIMPLib/Common/Constants.h"
#include "SIMPLib/Common/TemplateHelpers.h"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/FilterParameters/AbstractFilterParametersReader.h"
#include "SIMPLib/FilterParameters/AttributeMatrixSelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/BooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArrayCreationFilterParameter.h"
#include "SIMPLib/FilterParameters/DataArraySelectionFilterParameter.h"
#include "SIMPLib/FilterParameters/LinkedBooleanFilterParameter.h"
#include "SIMPLib/FilterParameters/SeparatorFilterParameter.h"
#include "SIMPLib/FilterParameters/StringFilterParameter.h"
#include "SIMPLib/Math/SIMPLibMath.h"

#include "DREAM3DReview/DREAM3DReviewConstants.h"
#include "DREAM3DReview/DREAM3DReviewVersion.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindArrayStatistics::FindArrayStatistics()
: m_FindLength(false)
, m_FindMin(false)
, m_FindMax(false)
, m_FindMean(false)
, m_FindMedian(false)
, m_FindStdDeviation(false)
, m_FindSummation(false)
, m_UseMask(false)
, m_StandardizeData(false)
, m_ComputeByIndex(false)
, m_DestinationAttributeMatrix("", "", "")
, m_MaskArrayPath("", "", "Mask")
, m_LengthArrayName("Length")
, m_MinimumArrayName("Minimum")
, m_MaximumArrayName("Maximum")
, m_MeanArrayName("Mean")
, m_MedianArrayName("Median")
, m_StdDeviationArrayName("StandardDeviation")
, m_SummationArrayName("Summation")
, m_StandardizedArrayName("Standardized")
, m_SelectedArrayPath("", "", "")
, m_FeatureIdsArrayPath("", "", "")
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
FindArrayStatistics::~FindArrayStatistics()
= default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayStatistics::setupFilterParameters()
{
  FilterParameterVector parameters;
  parameters.push_back(SeparatorFilterParameter::New("Statistics Options", FilterParameter::Parameter));
  QStringList linkedProps("LengthArrayName");
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Length", FindLength, FilterParameter::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "MinimumArrayName";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Minimum", FindMin, FilterParameter::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "MaximumArrayName";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Maximum", FindMax, FilterParameter::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "MeanArrayName";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Mean", FindMean, FilterParameter::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "MedianArrayName";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Median", FindMedian, FilterParameter::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "StdDeviationArrayName";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Standard Deviation", FindStdDeviation, FilterParameter::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "SummationArrayName";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Find Summation", FindSummation, FilterParameter::Parameter, FindArrayStatistics, linkedProps));
  AttributeMatrixSelectionFilterParameter::RequirementType amReq = AttributeMatrixSelectionFilterParameter::CreateRequirement(AttributeMatrix::Type::Any, IGeometry::Type::Any);
  linkedProps.clear();

  parameters.push_back(SeparatorFilterParameter::New("Algorithm Options", FilterParameter::Parameter));
  linkedProps << "MaskArrayPath";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Use Mask", UseMask, FilterParameter::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "FeatureIdsArrayPath";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Compute Statistics Per Feature/Ensemble", ComputeByIndex, FilterParameter::Parameter, FindArrayStatistics, linkedProps));
  linkedProps.clear();
  linkedProps << "StandardizedArrayName";
  parameters.push_back(SIMPL_NEW_LINKED_BOOL_FP("Standardize Data", StandardizeData, FilterParameter::Parameter, FindArrayStatistics, linkedProps));

  DataArraySelectionFilterParameter::RequirementType dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::Defaults::AnyPrimitive, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Attribute Array to Compute Statistics", SelectedArrayPath, FilterParameter::RequiredArray, FindArrayStatistics, dasReq));
  dasReq = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Int32, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Feature Ids", FeatureIdsArrayPath, FilterParameter::RequiredArray, FindArrayStatistics, dasReq));
  DataArraySelectionFilterParameter::RequirementType req = DataArraySelectionFilterParameter::CreateRequirement(SIMPL::TypeNames::Bool, 1, AttributeMatrix::Type::Any, IGeometry::Type::Any);
  parameters.push_back(SIMPL_NEW_DA_SELECTION_FP("Mask", MaskArrayPath, FilterParameter::RequiredArray, FindArrayStatistics, req));
  parameters.push_back(SIMPL_NEW_AM_SELECTION_FP("Destination Attribute Matrix", DestinationAttributeMatrix, FilterParameter::RequiredArray, FindArrayStatistics, amReq));

  parameters.push_back(SIMPL_NEW_STRING_FP("Length", LengthArrayName, FilterParameter::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_STRING_FP("Minimum", MinimumArrayName, FilterParameter::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_STRING_FP("Maximum", MaximumArrayName, FilterParameter::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_STRING_FP("Mean", MeanArrayName, FilterParameter::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_STRING_FP("Median", MedianArrayName, FilterParameter::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_STRING_FP("Standard Deviation", StdDeviationArrayName, FilterParameter::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_STRING_FP("Summation", SummationArrayName, FilterParameter::CreatedArray, FindArrayStatistics));
  parameters.push_back(SIMPL_NEW_STRING_FP("Standardized Data", StandardizedArrayName, FilterParameter::CreatedArray, FindArrayStatistics));

  setFilterParameters(parameters);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayStatistics::initialize()
{
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayStatistics::dataCheck()
{
  setErrorCondition(0);
  setWarningCondition(0);

  if(!getFindMin() && !getFindMax() && !getFindMean() && !getFindMedian() && !getFindStdDeviation() && !getFindSummation() && !getFindLength())
  {
    QString ss = QObject::tr("No statistics have been selected, so this filter will perform no operations");
    setWarningCondition(-701);
    notifyWarningMessage(getHumanLabel(), ss, getErrorCondition());
    return;
  }

  QVector<DataArrayPath> dataArrayPaths;

  m_InputArrayPtr = getDataContainerArray()->getPrereqIDataArrayFromPath<IDataArray, AbstractFilter>(this, getSelectedArrayPath());

  if(getErrorCondition() < 0)
  {
    return;
  }

  dataArrayPaths.push_back(getSelectedArrayPath());

  if(m_InputArrayPtr.lock()->getNumberOfComponents() != 1)
  {
    setErrorCondition(-11002);
    QString ss = QObject::tr("Input Attribute Array must be a scalar array");
    notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
  }

  if(!getComputeByIndex())
  {
    AttributeMatrix::Pointer destAttrMat = getDataContainerArray()->getPrereqAttributeMatrixFromPath<AbstractFilter>(this, getDestinationAttributeMatrix(), -301);
    if(getErrorCondition() < 0)
    {
      return;
    }
    QVector<size_t> tDims = destAttrMat->getTupleDimensions();
    if(tDims.size() != 1)
    {
      setErrorCondition(-11002);
      QString ss =
          QObject::tr("Since option \"Compute Statistics Per Feature/Ensemble Id\" is not selected, a single value, representative of the whole array, will be computed for each scalar statistic. "
                      "The selected destination Attribute Matrix (%1) must then have exactly 1 dimension, but the current selection has dimensions %2. "
                      "Consider creating a new Generic Attribute Matrix with scalar tuple dimensions to store the statistics.")
              .arg(getDestinationAttributeMatrix().getAttributeMatrixName())
              .arg(tDims.size());
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      return;
    }
    if(tDims[0] != 1)
    {
      setErrorCondition(-11002);
      QString ss =
          QObject::tr("Since option \"Compute Statistics Per Feature/Ensemble Id\" is not selected, a single value, representative of the whole array, will be computed for each scalar statistic. "
                      "The selected destination Attribute Matrix (%1) must then have an extent of 1 in its single dimension , but the current extent is %2. "
                      "Consider creating a new Generic Attribute Matrix with scalar tuple dimensions to store the statistics.")
              .arg(getDestinationAttributeMatrix().getAttributeMatrixName())
              .arg(tDims[0]);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      return;
    }
  }

  QVector<size_t> cDims(1, 1);

  if(getComputeByIndex())
  {
    m_FeatureIdsPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<int32_t>, AbstractFilter>(this, getFeatureIdsArrayPath(), cDims);
    if(m_FeatureIdsPtr.lock())
    {
      m_FeatureIds = m_FeatureIdsPtr.lock()->getPointer(0);
    }
    if(getErrorCondition() >= 0)
    {
      dataArrayPaths.push_back(getFeatureIdsArrayPath());
    }
  }

  EXECUTE_FUNCTION_TEMPLATE(this, createCompatibleArrays, m_InputArrayPtr.lock())

  if(m_FindLength)
  {
    DataArrayPath path(getDestinationAttributeMatrix().getDataContainerName(), getDestinationAttributeMatrix().getAttributeMatrixName(), getLengthArrayName());
    m_LengthPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<int64_t>, AbstractFilter, int64_t>(this, path, 0, cDims);
    if(m_LengthPtr.lock())
    {
      m_Length = m_LengthPtr.lock()->getPointer(0);
    }
  }

  if(getUseMask())
  {
    m_MaskPtr = getDataContainerArray()->getPrereqArrayFromPath<DataArray<bool>, AbstractFilter>(this, getMaskArrayPath(), cDims);
    if(m_MaskPtr.lock())
    {
      m_Mask = m_MaskPtr.lock()->getPointer(0);
    }
    if(getErrorCondition() >= 0)
    {
      dataArrayPaths.push_back(getMaskArrayPath());
    }
  }

  if(getStandardizeData())
  {
    if(!m_FindMean || !m_FindStdDeviation)
    {
      setErrorCondition(-11003);
      QString ss = QObject::tr("To standardize data, the \"Find Mean\" and \"Find Standard Deviation\" options must also be checked");
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
    }
    DataArrayPath path(getSelectedArrayPath().getDataContainerName(), getSelectedArrayPath().getAttributeMatrixName(), getStandardizedArrayName());
    m_StandardizedPtr = getDataContainerArray()->createNonPrereqArrayFromPath<DataArray<float>, AbstractFilter, float>(this, path, 0, cDims);
    if(m_StandardizedPtr.lock())
    {
      m_Standardized = m_StandardizedPtr.lock()->getPointer(0);
    }
    if(getErrorCondition() >= 0)
    {
      dataArrayPaths.push_back(path);
    }
  }

  getDataContainerArray()->validateNumberOfTuples<AbstractFilter>(this, dataArrayPaths);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayStatistics::preflight()
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
template <template <typename, typename...> class C, typename T, typename... Ts> int64_t findLength(C<T, Ts...>& source)
{
  return static_cast<int64_t>(source.size());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts> T findMin(C<T, Ts...>& source)
{
  if(source.empty())
  {
    return T(0);
  }
  return (*std::min_element(std::begin(source), std::end(source)));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts> T findMax(C<T, Ts...>& source)
{
  if(source.empty())
  {
    return T(0);
  }
  return (*std::max_element(std::begin(source), std::end(source)));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts> float findMean(C<T, Ts...>& source)
{
  if(source.empty())
  {
    return 0.0f;
  }
  float sum = std::accumulate(std::begin(source), std::end(source), 0.0f);
  return static_cast<float>(sum / source.size());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename... Ts> bool findMean(C<bool, Ts...>& source)
{
  if(source.empty())
  {
    return false;
  }
  size_t count = std::count(std::begin(source), std::end(source), true);
  return true ? count >= (source.size() - count) : false;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts> float findMedian(C<T, Ts...>& source)
{
  // Need a copy, not a reference, since we will be messing with the vector order
  std::vector<T> tmpList{std::begin(source), std::end(source)};
  if(tmpList.empty())
  {
    return 0.0f;
  }
  auto halfElements = tmpList.size() / 2;
  std::nth_element(std::begin(tmpList), tmpList.begin() + halfElements, std::end(tmpList));
  T medVal = tmpList[halfElements];
  if(tmpList.size() % 2 == 1)
  {
    return static_cast<float>(medVal);
  }
  
  
    return static_cast<float>((0.5f * (medVal + tmpList[halfElements - 1])));
  
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts> float findStdDeviation(C<T, Ts...>& source)
{
  if(source.empty())
  {
    return 0.0f;
  }
  std::vector<double> difference(source.size());
  float sum = std::accumulate(std::begin(source), std::end(source), 0.0f);
  float mean = static_cast<double>(sum / source.size());
  std::transform(std::begin(source), std::end(source), std::begin(difference), std::bind2nd(std::minus<float>(), mean));
  float squaredSum = std::inner_product(std::begin(difference), std::end(difference), std::begin(difference), 0.0f);
  return std::sqrt(squaredSum / source.size());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename... Ts> bool findStdDeviation(C<bool, Ts...>& source)
{
  if(source.empty())
  {
    return false;
  }
  size_t count = std::count(std::begin(source), std::end(source), true);
  return true ? count >= (source.size() - count) : false;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <template <typename, typename...> class C, typename T, typename... Ts> double findSummation(C<T, Ts...>& source)
{
  if(source.empty())
  {
    return 0.0f;
  }
  float sum = std::accumulate(std::begin(source), std::end(source), 0.0f);
  return sum;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void standardizeDataByIndex(IDataArray::Pointer dataPtr, FloatArrayType::Pointer standardized, bool useMask, bool* mask, int32_t* featureIds, IDataArray::Pointer meanPtr, IDataArray::Pointer stdPtr)
{
  typename DataArray<T>::Pointer inDataPtr = std::dynamic_pointer_cast<DataArray<T>>(dataPtr);
  FloatArrayType::Pointer muPtr = std::dynamic_pointer_cast<FloatArrayType>(meanPtr);
  FloatArrayType::Pointer sigPtr = std::dynamic_pointer_cast<FloatArrayType>(stdPtr);
  FloatArrayType::Pointer standardizedPtr = std::dynamic_pointer_cast<FloatArrayType>(standardized);

  T* dPtr = inDataPtr->getPointer(0);
  float* mPtr = muPtr->getPointer(0);
  float* sPtr = sigPtr->getPointer(0);
  float* stPtr = standardizedPtr->getPointer(0);
  size_t numTuples = inDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    if(useMask)
    {
      if(mask[i])
      {
        stPtr[i] = (static_cast<float>(dPtr[i]) - mPtr[featureIds[i]]) / sPtr[featureIds[i]];
      }
    }
    else
    {
      stPtr[i] = (static_cast<float>(dPtr[i]) - mPtr[featureIds[i]]) / sPtr[featureIds[i]];
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <>
void standardizeDataByIndex<bool>(IDataArray::Pointer dataPtr, FloatArrayType::Pointer standardized, bool useMask, bool* mask, int32_t* featureIds, IDataArray::Pointer meanPtr,
                                  IDataArray::Pointer stdPtr)
{
  // Standardization of a boolean array is a no-op
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> void standardizeData(IDataArray::Pointer dataPtr, FloatArrayType::Pointer standardized, bool useMask, bool* mask, IDataArray::Pointer meanPtr, IDataArray::Pointer stdPtr)
{
  typename DataArray<T>::Pointer inDataPtr = std::dynamic_pointer_cast<DataArray<T>>(dataPtr);
  FloatArrayType::Pointer muPtr = std::dynamic_pointer_cast<FloatArrayType>(meanPtr);
  FloatArrayType::Pointer sigPtr = std::dynamic_pointer_cast<FloatArrayType>(stdPtr);
  FloatArrayType::Pointer standardizedPtr = std::dynamic_pointer_cast<FloatArrayType>(standardized);

  T* dPtr = inDataPtr->getPointer(0);
  float* mPtr = muPtr->getPointer(0);
  float* sPtr = sigPtr->getPointer(0);
  float* stPtr = standardizedPtr->getPointer(0);
  size_t numTuples = inDataPtr->getNumberOfTuples();

  for(size_t i = 0; i < numTuples; i++)
  {
    if(useMask)
    {
      if(mask[i])
      {
        stPtr[i] = (static_cast<float>(dPtr[i]) - mPtr[0]) / sPtr[0];
      }
    }
    else
    {
      stPtr[i] = (static_cast<float>(dPtr[i]) - mPtr[0]) / sPtr[0];
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <> void standardizeData<bool>(IDataArray::Pointer dataPtr, FloatArrayType::Pointer standardized, bool useMask, bool* mask, IDataArray::Pointer meanPtr, IDataArray::Pointer stdPtr)
{
  // Standardization of a boolean array is a no-op
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T> class FindStatisticsByIndexImpl
{
public:
  FindStatisticsByIndexImpl(std::unordered_map<int32_t, std::list<T>>& featureDataMap, bool length, bool min, bool max, bool mean, bool median, bool stdDeviation, bool summation,
                            std::vector<IDataArray::Pointer>& arrays)
  : m_FeatureDataMap(featureDataMap)
  , m_Length(length)
  , m_Min(min)
  , m_Max(max)
  , m_Mean(mean)
  , m_Median(median)
  , m_StdDeviation(stdDeviation)
  , m_Summation(summation)
  , m_Arrays(arrays)
  {
  }

  virtual ~FindStatisticsByIndexImpl() = default;

  void compute(size_t start, size_t end) const
  {
    for(size_t i = start; i < end; i++)
    {
      if(m_Length)
      {
        if(m_Arrays[0])
        {
          int64_t val = findLength(m_FeatureDataMap[i]);
          m_Arrays[0]->initializeTuple(i, &val);
        }
      }
      if(m_Min)
      {
        if(m_Arrays[1])
        {
          T val = findMin(m_FeatureDataMap[i]);
          m_Arrays[1]->initializeTuple(i, &val);
        }
      }
      if(m_Max)
      {
        if(m_Arrays[2])
        {
          T val = findMax(m_FeatureDataMap[i]);
          m_Arrays[2]->initializeTuple(i, &val);
        }
      }
      if(m_Mean)
      {
        if(m_Arrays[3])
        {
          float val = findMean(m_FeatureDataMap[i]);
          m_Arrays[3]->initializeTuple(i, &val);
        }
      }
      if(m_Median)
      {
        if(m_Arrays[4])
        {
          float val = findMedian(m_FeatureDataMap[i]);
          m_Arrays[4]->initializeTuple(i, &val);
        }
      }
      if(m_StdDeviation)
      {
        if(m_Arrays[5])
        {
          float val = findStdDeviation(m_FeatureDataMap[i]);
          m_Arrays[5]->initializeTuple(i, &val);
        }
      }
      if(m_Summation)
      {
        if(m_Arrays[6])
        {
          float val = findSummation(m_FeatureDataMap[i]);
          m_Arrays[6]->initializeTuple(i, &val);
        }
      }
    }
  }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
  void operator()(const tbb::blocked_range<size_t>& r) const
  {
    compute(r.begin(), r.end());
  }
#endif

private:
  std::unordered_map<int32_t, std::list<T>>& m_FeatureDataMap;
  bool m_Length;
  bool m_Min;
  bool m_Max;
  bool m_Mean;
  bool m_Median;
  bool m_StdDeviation;
  bool m_Summation;
  std::vector<IDataArray::Pointer>& m_Arrays;
};

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void findStatisticsImpl(bool length, bool min, bool max, bool mean, bool median, bool stdDeviation, bool summation, std::vector<IDataArray::Pointer>& arrays, std::vector<T>& data)
{
  if(length)
  {
    if(arrays[0])
    {
      int64_t val = findLength(data);
      arrays[0]->initializeTuple(0, &val);
    }
  }
  if(min)
  {
    if(arrays[1])
    {
      T val = findMin(data);
      arrays[1]->initializeTuple(0, &val);
    }
  }
  if(max)
  {
    if(arrays[2])
    {
      T val = findMax(data);
      arrays[2]->initializeTuple(0, &val);
    }
  }
  if(mean)
  {
    if(arrays[3])
    {
      float val = findMean(data);
      arrays[3]->initializeTuple(0, &val);
    }
  }
  if(median)
  {
    if(arrays[4])
    {
      float val = findMedian(data);
      arrays[4]->initializeTuple(0, &val);
    }
  }
  if(stdDeviation)
  {
    if(arrays[5])
    {
      float val = findStdDeviation(data);
      arrays[5]->initializeTuple(0, &val);
    }
  }
  if(summation)
  {
    if(arrays[6])
    {
      float val = findSummation(data);
      arrays[6]->initializeTuple(0, &val);
    }
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
template <typename T>
void findStatistics(IDataArray::Pointer source, Int32ArrayType::Pointer featureIds, bool useMask, bool* mask, bool length, bool min, bool max, bool mean, bool median, bool stdDeviation,
                    bool summation, std::vector<IDataArray::Pointer>& arrays, int32_t numFeatures, bool computeByIndex)
{
  size_t numTuples = source->getNumberOfTuples();
  typename DataArray<T>::Pointer sourcePtr = std::dynamic_pointer_cast<DataArray<T>>(source);
  T* dataPtr = sourcePtr->getPointer(0);

  if(computeByIndex)
  {
    int32_t* featureIdsPtr = featureIds->getPointer(0);
    std::unordered_map<int32_t, std::list<T>> featureValueMap;

    for(size_t i = 0; i < numTuples; i++)
    {
      if(useMask)
      {
        if(mask[i])
        {
          featureValueMap[featureIdsPtr[i]].push_back(dataPtr[i]);
        }
      }
      else
      {
        featureValueMap[featureIdsPtr[i]].push_back(dataPtr[i]);
      }
    }

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    tbb::task_scheduler_init init;
    bool doParallel = true;
#endif

#ifdef SIMPL_USE_PARALLEL_ALGORITHMS
    if(doParallel)
    {
      tbb::parallel_for(tbb::blocked_range<size_t>(0, numFeatures), FindStatisticsByIndexImpl<T>(featureValueMap, length, min, max, mean, median, stdDeviation, summation, arrays),
                        tbb::auto_partitioner());
    }
    else
#endif
    {
      FindStatisticsByIndexImpl<T> serial(featureValueMap, length, min, max, mean, median, stdDeviation, summation, arrays);
      serial.compute(0, numFeatures);
    }
  }
  else
  {
    std::vector<T> data;
    data.reserve(numTuples);
    for(size_t i = 0; i < numTuples; i++)
    {
      if(useMask)
      {
        if(mask[i])
        {
          data.push_back(dataPtr[i]);
        }
      }
      else
      {
        data.push_back(dataPtr[i]);
      }
    }

    data.shrink_to_fit();

    findStatisticsImpl(length, min, max, mean, median, stdDeviation, summation, arrays, data);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void FindArrayStatistics::execute()
{
  setErrorCondition(0);
  setWarningCondition(0);
  dataCheck();
  if(getErrorCondition() < 0)
  {
    return;
  }

  if(!m_FindMin && !m_FindMax && !m_FindMean && !m_FindMedian && !m_FindStdDeviation && !m_FindSummation && !m_FindLength)
  {
    return;
  }

  int32_t numFeatures = 0;

  if(m_ComputeByIndex)
  {
    AttributeMatrix::Pointer attrMat =
        getDataContainerArray()->getDataContainer(m_DestinationAttributeMatrix.getDataContainerName())->getAttributeMatrix(m_DestinationAttributeMatrix.getAttributeMatrixName());
    numFeatures = static_cast<int32_t>(attrMat->getNumberOfTuples());
    bool mismatchedFeatures = false;
    int32_t largestFeature = 0;
    size_t totalPoints = m_FeatureIdsPtr.lock()->getNumberOfTuples();

    for(size_t i = 0; i < totalPoints; i++)
    {
      if(m_FeatureIds[i] > largestFeature)
      {
        largestFeature = m_FeatureIds[i];
        if(largestFeature >= numFeatures)
        {
          mismatchedFeatures = true;
          break;
        }
      }
    }

    if(mismatchedFeatures)
    {
      QString ss = QObject::tr("The number of objects in the selected Attribute Matrix destination (%1) is larger than the largest Id in the Feature/Ensemble Ids array").arg(numFeatures);
      setErrorCondition(-5555);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      return;
    }

    if(largestFeature != (numFeatures - 1))
    {
      QString ss = QObject::tr("The number of objects in the selected Attribute Matrix destination (%1) does not match the largest Id in the  Feature/Ensemble Ids array").arg(numFeatures);
      setErrorCondition(-5555);
      notifyErrorMessage(getHumanLabel(), ss, getErrorCondition());
      return;
    }
  }

  std::vector<IDataArray::Pointer> arrays(7);
  std::fill(std::begin(arrays), std::end(arrays), nullptr);

  for(size_t i = 0; i < 7; i++)
  {
    if(m_FindLength)
    {
      arrays[0] = m_LengthPtr.lock();
    }
    if(m_FindMin)
    {
      arrays[1] = m_MinimumPtr.lock();
    }
    if(m_FindMax)
    {
      arrays[2] = m_MaximumPtr.lock();
    }
    if(m_FindMean)
    {
      arrays[3] = m_MeanPtr.lock();
    }
    if(m_FindMedian)
    {
      arrays[4] = m_MedianPtr.lock();
    }
    if(m_FindStdDeviation)
    {
      arrays[5] = m_StandardDeviationPtr.lock();
    }
    if(m_FindSummation)
    {
      arrays[6] = m_SummationPtr.lock();
    }
  }

  EXECUTE_FUNCTION_TEMPLATE(this, findStatistics, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_FeatureIdsPtr.lock(), m_UseMask, m_Mask, m_FindLength, m_FindMin, m_FindMax, m_FindMean,
                            m_FindMedian, m_FindStdDeviation, m_FindSummation, arrays, numFeatures, m_ComputeByIndex);

  if(m_StandardizeData)
  {
    IDataArray::Pointer meanPtr = getDataContainerArray()->getAttributeMatrix(m_DestinationAttributeMatrix)->getAttributeArray(getMeanArrayName());
    IDataArray::Pointer stdPtr = getDataContainerArray()->getAttributeMatrix(m_DestinationAttributeMatrix)->getAttributeArray(getStdDeviationArrayName());

    if(m_ComputeByIndex)
    {
      EXECUTE_FUNCTION_TEMPLATE(this, standardizeDataByIndex, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_StandardizedPtr.lock(), m_UseMask, m_Mask, m_FeatureIds, meanPtr, stdPtr)
    }
    else
    {
      EXECUTE_FUNCTION_TEMPLATE(this, standardizeData, m_InputArrayPtr.lock(), m_InputArrayPtr.lock(), m_StandardizedPtr.lock(), m_UseMask, m_Mask, meanPtr, stdPtr)
    }
  }

}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
AbstractFilter::Pointer FindArrayStatistics::newFilterInstance(bool copyFilterParameters) const
{
  FindArrayStatistics::Pointer filter = FindArrayStatistics::New();
  if(copyFilterParameters)
  {
    copyFilterParameterInstanceVariables(filter.get());
  }
  return filter;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindArrayStatistics::getCompiledLibraryName() const
{
  return DREAM3DReviewConstants::DREAM3DReviewBaseName;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindArrayStatistics::getBrandingString() const
{
  return "DREAM3DReview";
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindArrayStatistics::getFilterVersion() const
{
  QString version;
  QTextStream vStream(&version);
  vStream << DREAM3DReview::Version::Major() << "." << DREAM3DReview::Version::Minor() << "." << DREAM3DReview::Version::Patch();
  return version;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindArrayStatistics::getGroupName() const
{
  return DREAM3DReviewConstants::FilterGroups::DREAM3DReviewFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QUuid FindArrayStatistics::getUuid()
{
  return QUuid("{bf35f515-294b-55ed-8c69-211b7e69cb56}");
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindArrayStatistics::getSubGroupName() const
{
  return DREAM3DReviewConstants::FilterSubGroups::StatisticsFilters;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
const QString FindArrayStatistics::getHumanLabel() const
{
  return "Find Attribute Array Statistics";
}
