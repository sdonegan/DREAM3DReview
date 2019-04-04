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
#ifndef _mappointcloudtoregulargrid_h_
#define _mappointcloudtoregulargrid_h_

#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/FilterParameters/IntVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

/**
 * @brief The MapPointCloudToRegularGrid class. See [Filter documentation](@ref mappointcloudtoregulargrid) for details.
 */
class MapPointCloudToRegularGrid : public AbstractFilter
{
  Q_OBJECT

public:
  SIMPL_SHARED_POINTERS(MapPointCloudToRegularGrid)
  SIMPL_FILTER_NEW_MACRO(MapPointCloudToRegularGrid)
  SIMPL_TYPE_MACRO_SUPER_OVERRIDE(MapPointCloudToRegularGrid, AbstractFilter)

  ~MapPointCloudToRegularGrid() override;

  SIMPL_FILTER_PARAMETER(DataArrayPath, DataContainerName)
  Q_PROPERTY(DataArrayPath DataContainerName READ getDataContainerName WRITE setDataContainerName)

  SIMPL_FILTER_PARAMETER(QString, ImageDataContainerName)
  Q_PROPERTY(QString ImageDataContainerName READ getImageDataContainerName WRITE setImageDataContainerName)

  SIMPL_FILTER_PARAMETER(DataArrayPath, ImageDataContainerPath)
  Q_PROPERTY(DataArrayPath ImageDataContainerPath READ getImageDataContainerPath WRITE setImageDataContainerPath)

  SIMPL_FILTER_PARAMETER(DataArrayPath, VoxelIndicesArrayPath)
  Q_PROPERTY(DataArrayPath VoxelIndicesArrayPath READ getVoxelIndicesArrayPath WRITE setVoxelIndicesArrayPath)

  SIMPL_FILTER_PARAMETER(IntVec3Type, GridDimensions)
  Q_PROPERTY(IntVec3Type GridDimensions READ getGridDimensions WRITE setGridDimensions)

  SIMPL_FILTER_PARAMETER(bool, UseMask)
  Q_PROPERTY(bool UseMask READ getUseMask WRITE setUseMask)

  SIMPL_FILTER_PARAMETER(int, CreateDataContainer)
  Q_PROPERTY(int CreateDataContainer READ getCreateDataContainer WRITE setCreateDataContainer)

  SIMPL_FILTER_PARAMETER(DataArrayPath, MaskArrayPath)
  Q_PROPERTY(DataArrayPath MaskArrayPath READ getMaskArrayPath WRITE setMaskArrayPath)

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
  MapPointCloudToRegularGrid();

  /**
   * @brief createRegularGrid Creates the structured rectilinear grid to interpolate
   * the point cloud onto
   */
  void createRegularGrid();

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

  std::vector<float> m_MeshMinExtents;
  std::vector<float> m_MeshMaxExtents;

  MapPointCloudToRegularGrid(const MapPointCloudToRegularGrid&) = delete; // Copy Constructor Not Implemented
  MapPointCloudToRegularGrid(MapPointCloudToRegularGrid&&) = delete;      // Move Constructor Not Implemented
  void operator=(const MapPointCloudToRegularGrid&) = delete;             // Operator '=' Not Implemented
};

#endif /* _MapPointCloudToRegularGrid_H_ */
