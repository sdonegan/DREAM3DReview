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
#ifndef _SliceTriangleGeometry_h_
#define _SliceTriangleGeometry_h_

#include "SIMPLib/Common/SIMPLibSetGetMacros.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/SIMPLib.h"

/**
 * @brief The SliceTriangleGeometry class. See [Filter documentation](@ref SliceTriangleGeometry) for details.
 */
class SliceTriangleGeometry : public AbstractFilter
{
  Q_OBJECT

public:
  SIMPL_SHARED_POINTERS(SliceTriangleGeometry)
  SIMPL_FILTER_NEW_MACRO(SliceTriangleGeometry)
  SIMPL_TYPE_MACRO_SUPER_OVERRIDE(SliceTriangleGeometry, AbstractFilter)

  ~SliceTriangleGeometry() override;

  SIMPL_FILTER_PARAMETER(DataArrayPath, CADDataContainerName)
  Q_PROPERTY(DataArrayPath CADDataContainerName READ getCADDataContainerName WRITE setCADDataContainerName)

  SIMPL_FILTER_PARAMETER(QString, SliceDataContainerName)
  Q_PROPERTY(QString SliceDataContainerName READ getSliceDataContainerName WRITE setSliceDataContainerName)

  SIMPL_FILTER_PARAMETER(QString, EdgeAttributeMatrixName)
  Q_PROPERTY(QString EdgeAttributeMatrixName READ getEdgeAttributeMatrixName WRITE setEdgeAttributeMatrixName)

  SIMPL_FILTER_PARAMETER(QString, SliceAttributeMatrixName)
  Q_PROPERTY(QString SliceAttributeMatrixName READ getSliceAttributeMatrixName WRITE setSliceAttributeMatrixName)

  SIMPL_FILTER_PARAMETER(QString, SliceIdArrayName)
  Q_PROPERTY(QString SliceIdArrayName READ getSliceIdArrayName WRITE setSliceIdArrayName)

  SIMPL_FILTER_PARAMETER(QString, AreasArrayName)
  Q_PROPERTY(QString AreasArrayName READ getAreasArrayName WRITE setAreasArrayName)

  SIMPL_FILTER_PARAMETER(QString, PerimetersArrayName)
  Q_PROPERTY(QString PerimetersArrayName READ getPerimetersArrayName WRITE setPerimetersArrayName)

  SIMPL_FILTER_PARAMETER(bool, HaveRegionIds)
  Q_PROPERTY(bool HaveRegionIds READ getHaveRegionIds WRITE setHaveRegionIds)

  SIMPL_FILTER_PARAMETER(DataArrayPath, RegionIdArrayPath)
  Q_PROPERTY(DataArrayPath RegionIdArrayPath READ getRegionIdArrayPath WRITE setRegionIdArrayPath)

  SIMPL_FILTER_PARAMETER(FloatVec3Type, SliceDirection)
  Q_PROPERTY(FloatVec3Type SliceDirection READ getSliceDirection WRITE setSliceDirection)

  SIMPL_FILTER_PARAMETER(float, SliceResolution)
  Q_PROPERTY(float SliceResolution READ getSliceResolution WRITE setSliceResolution)

  SIMPL_FILTER_PARAMETER(float, Zstart)
  Q_PROPERTY(float Zstart READ getZstart WRITE setZstart)

  SIMPL_FILTER_PARAMETER(float, Zend)
  Q_PROPERTY(float Zend READ getZend WRITE setZend)

  SIMPL_FILTER_PARAMETER(int, SliceRange)
  Q_PROPERTY(int SliceRange READ getSliceRange WRITE setSliceRange)

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
  SliceTriangleGeometry();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck();

  /**
   * @brief rotateVertices either rotates CAD triangles and/or edges into or out of the sectioning orientation so the sectioning is always along 001
   */
  void rotateVertices(unsigned int direction, float* n, int64_t numVerts, float* verts);

  /**
   * @brief updateEdgeInstancePointers
   */
  void determineBoundsAndNumSlices(float& minDim, float& maxDim, int64_t numTris, int64_t* tris, float* triVerts);

  /**
  * @brief updateEdgeInstancePointers
  */
  void updateEdgeInstancePointers();

  /**
  * @brief updateSliceInstancePointers
  */
  void updateSliceInstancePointers();

private:
  enum RotationDirection
  {
    rotForward,
    rotBackward
  };

  int32_t m_NumberOfSlices;

  DEFINE_DATAARRAY_VARIABLE(int32_t, SliceId)
  DEFINE_DATAARRAY_VARIABLE(float, Area)
  DEFINE_DATAARRAY_VARIABLE(float, Perimeter)
  DEFINE_DATAARRAY_VARIABLE(int32_t, RegionId)
  DEFINE_DATAARRAY_VARIABLE(int32_t, TriRegionId)

  SliceTriangleGeometry(const SliceTriangleGeometry&) = delete;   // Copy Constructor Not Implemented
  SliceTriangleGeometry(SliceTriangleGeometry&&) = delete;        // Move Constructor Not Implemented
  void operator=(const SliceTriangleGeometry&) = delete; // Operator '=' Not Implemented
};

#endif /* _SliceTriangleGeometry_H_ */
