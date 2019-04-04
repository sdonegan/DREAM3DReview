/* ============================================================================
* Software developed by US federal government employees (including military personnel) 
* as part of their official duties is not subject to copyright protection and is 
* considered “public domain” (see 17 USC Section 105). Public domain software can be used 
* by anyone for any purpose, and cannot be released under a copyright license 
* (including typical open source software licenses).
* 
* This source code file was originally written by United States DoD employees. The
* original source code files are released into the Public Domain.
* 
* Subsequent changes to the codes by others may elect to add a copyright and license
* for those changes.
* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
#ifndef _interpolatepointcloudtoregulargrid_h_
#define _interpolatepointcloudtoregulargrid_h_

#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/DataArrays/NeighborList.hpp"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/SIMPLib.h"

/**
 * @brief The InterpolatePointCloudToRegularGrid class. See [Filter documentation](@ref interpolatepointcloudtoregulargrid) for details.
 */
class InterpolatePointCloudToRegularGrid : public AbstractFilter
{
  Q_OBJECT

public:
  SIMPL_SHARED_POINTERS(InterpolatePointCloudToRegularGrid)
  SIMPL_FILTER_NEW_MACRO(InterpolatePointCloudToRegularGrid)
  SIMPL_TYPE_MACRO_SUPER_OVERRIDE(InterpolatePointCloudToRegularGrid, AbstractFilter)

  ~InterpolatePointCloudToRegularGrid() override;

  SIMPL_FILTER_PARAMETER(DataArrayPath, DataContainerName)
  Q_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)

  SIMPL_FILTER_PARAMETER(QVector<DataArrayPath>, ArraysToInterpolate)
  Q_PROPERTY(QVector<DataArrayPath> ArraysToInterpolate READ getArraysToInterpolate WRITE setArraysToInterpolate)

  SIMPL_FILTER_PARAMETER(QVector<DataArrayPath>, ArraysToCopy)
  Q_PROPERTY(QVector<DataArrayPath> ArraysToCopy READ getArraysToCopy WRITE setArraysToCopy)

  SIMPL_FILTER_PARAMETER(DataArrayPath, VoxelIndicesArrayPath)
  Q_PROPERTY(DataArrayPath VoxelIndicesArrayPath READ getVoxelIndicesArrayPath WRITE setVoxelIndicesArrayPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, InterpolatedDataContainerName)
  Q_PROPERTY(DataArrayPath InterpolatedDataContainerName READ getInterpolatedDataContainerName WRITE setInterpolatedDataContainerName)

  SIMPL_FILTER_PARAMETER(QString, InterpolatedAttributeMatrixName)
  Q_PROPERTY(QString InterpolatedAttributeMatrixName READ getInterpolatedAttributeMatrixName WRITE setInterpolatedAttributeMatrixName)

  SIMPL_FILTER_PARAMETER(int, InterpolationTechnique)
  Q_PROPERTY(int InterpolationTechnique READ getInterpolationTechnique WRITE setInterpolationTechnique)

  SIMPL_FILTER_PARAMETER(FloatVec3Type, KernelSize)
  Q_PROPERTY(FloatVec3Type KernelSize READ getKernelSize WRITE setKernelSize)

  SIMPL_FILTER_PARAMETER(FloatVec3Type, Sigmas)
  Q_PROPERTY(FloatVec3Type Sigmas READ getSigmas WRITE setSigmas)

  SIMPL_FILTER_PARAMETER(bool, UseMask)
  Q_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)

  SIMPL_FILTER_PARAMETER(DataArrayPath, MaskArrayPath)
  Q_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)

  SIMPL_FILTER_PARAMETER(bool, StoreKernelDistances)
  Q_PROPERTY(bool StoreKernelDistances READ getStoreKernelDistances WRITE setStoreKernelDistances)

  SIMPL_FILTER_PARAMETER(QString, KernelDistancesArrayName)
  Q_PROPERTY(QString KernelDistancesArrayName READ getKernelDistancesArrayName WRITE setKernelDistancesArrayName)

  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  const QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
  */
  const QString getBrandingString() const override;

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  const QString getFilterVersion() const override;

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  const QString getGroupName() const override;

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  const QString getSubGroupName() const override;

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  const QString getHumanLabel() const override;

  /**
   * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
   */
  void setupFilterParameters() override;

  /**
   * @brief readFilterParameters Reimplemented from @see AbstractFilter class
   */
  void readFilterParameters(AbstractFilterParametersReader* reader, int index) override;

  /**
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
  * @brief preflight Reimplemented from @see AbstractFilter class
  */
  void preflight() override;

  /**
  * @brief getUuid Return the unique identifier for this filter.
  * @return A QUuid object.
  */
  const QUuid getUuid() override;

signals:
  /**
   * @brief updateFilterParameters Emitted when the Filter requests all the latest Filter parameters
   * be pushed from a user-facing control (such as a widget)
   * @param filter Filter instance pointer
   */
  void updateFilterParameters(AbstractFilter* filter);

  /**
   * @brief parametersChanged Emitted when any Filter parameter is changed internally
   */
  void parametersChanged();

  /**
   * @brief preflightAboutToExecute Emitted just before calling dataCheck()
   */
  void preflightAboutToExecute();

  /**
   * @brief preflightExecuted Emitted just after calling dataCheck()
   */
  void preflightExecuted();

protected:
  InterpolatePointCloudToRegularGrid();

  /**
   * @brief determineKernel Determines the kernel to be used for interpolation
   * @param kernelNumVoxels Voxel extents of the kernel
   * @return Float vector for the kernel values
   */
  void determineKernel(int64_t kernelNumVoxels[3]);

  /**
   * @brief determineKernelDistances Determines the Euclidean distance between each
   * voxel in the kernel and the kernel center
   * @param kernelNumVoxels Voxel extents of the kernel
   * @return Float vector for the kernel distances
   */
  void determineKernelDistances(int64_t kernelNumVoxels[3], FloatVec3Type res);

  /**
   * @brief mapKernelDistances Maps kernel point distances during interpolation
   * @param kernel Voxel dimensions of the kernel
   * @param dims Total dimensions of the interpolation grid
   * @param curX Current x position
   * @param curY Current y position
   * @param curZ Current z position
   */
  void mapKernelDistances(int64_t kernel[3], size_t dims[3], size_t curX, size_t curY, size_t curZ);
  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  DEFINE_DATAARRAY_VARIABLE(int64_t, VoxelIndices)
  DEFINE_DATAARRAY_VARIABLE(bool, Mask)

  NeighborList<float>::WeakPointer m_KernelDistances;

  QList<QString> m_AttrMatList;
  std::vector<IDataArray::WeakPointer> m_SourceArraysToInterpolate;
  std::vector<IDataArray::WeakPointer> m_SourceArraysToCopy;
  std::vector<IDataArray::WeakPointer> m_DynamicArraysToInterpolate;
  std::vector<IDataArray::WeakPointer> m_DynamicArraysToCopy;
  std::vector<float> m_Kernel;
  std::vector<float> m_KernelValDistances;

  InterpolatePointCloudToRegularGrid(const InterpolatePointCloudToRegularGrid&) = delete; // Copy Constructor Not Implemented
  InterpolatePointCloudToRegularGrid(InterpolatePointCloudToRegularGrid&&) = delete;      // Move Constructor Not Implemented
  void operator=(const InterpolatePointCloudToRegularGrid&) = delete;                     // Operator '=' Not Implemented
};

#endif /* _InterpolatePointCloudToRegularGrid_H_ */
