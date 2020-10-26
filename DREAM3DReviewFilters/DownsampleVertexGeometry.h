/* ============================================================================
 * Software developed by US federal government employees (including military personnel)
 * as part of their official duties is not subject to copyright protection and is
 * considered "public domain" (see 17 USC Section 105). Public domain software can be used
 * by anyone for any purpose, and cannot be released under a copyright license
 * (including typical open source software licenses).
 *
 * This source code file was originally written by United States DoD employees. The
 * original source code files are released into the Public Domain.
 *
 * Subsequent changes to the codes by others may elect to add a copyright and license
 * for those changes.
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

#pragma once

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/FilterParameters/FloatVec3FilterParameter.h"
#include "SIMPLib/Filtering/AbstractFilter.h"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The DownsampleVertexGeometry class. See [Filter documentation](@ref DownsampleVertexGeometry) for details.
 */
class DREAM3DReview_EXPORT DownsampleVertexGeometry : public AbstractFilter
{
  Q_OBJECT

  PYB11_BEGIN_BINDINGS(DownsampleVertexGeometry SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(DownsampleVertexGeometry)
  PYB11_FILTER_NEW_MACRO(DownsampleVertexGeometry)
  PYB11_PROPERTY(DataArrayPath VertexAttrMatPath READ getVertexAttrMatPath WRITE setVertexAttrMatPath)
  PYB11_PROPERTY(int DownsampleType READ getDownsampleType WRITE setDownsampleType)
  PYB11_PROPERTY(int DecimationFreq READ getDecimationFreq WRITE setDecimationFreq)
  PYB11_PROPERTY(float DecimationFraction READ getDecimationFraction WRITE setDecimationFraction)
  PYB11_PROPERTY(FloatVec3Type GridResolution READ getGridResolution WRITE setGridResolution)
  PYB11_END_BINDINGS()

public:
  using Self = DownsampleVertexGeometry;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<DownsampleVertexGeometry> New();

  /**
   * @brief Returns the name of the class for DownsampleVertexGeometry
   */
  QString getNameOfClass() const override;

  /**
   * @brief Returns the name of the class for DownsampleVertexGeometry
   */
  static QString ClassName();

  ~DownsampleVertexGeometry() override;

  /**
   * @brief Setter property for VertexAttrMatPath
   */
  void setVertexAttrMatPath(const DataArrayPath& value);

  /**
   * @brief Getter property for VertexAttrMatPath
   * @return Value of VertexAttrMatPath
   */
  DataArrayPath getVertexAttrMatPath() const;
  Q_PROPERTY(DataArrayPath VertexAttrMatPath READ getVertexAttrMatPath WRITE setVertexAttrMatPath)

  /**
   * @brief Setter property for DownsampleType
   */
  void setDownsampleType(const int& value);

  /**
   * @brief Getter property for DownsampleType
   * @return Value of DownsampleType
   */
  int getDownsampleType() const;
  Q_PROPERTY(int DownsampleType READ getDownsampleType WRITE setDownsampleType)

  /**
   * @brief Setter property for DecimationFreq
   */
  void setDecimationFreq(const int& value);

  /**
   * @brief Getter property for DecimationFreq
   * @return Value of DecimationFreq
   */
  int getDecimationFreq() const;
  Q_PROPERTY(int DecimationFreq READ getDecimationFreq WRITE setDecimationFreq)

  /**
   * @brief Setter property for DecimationFraction
   */
  void setDecimationFraction(const float& value);

  /**
   * @brief Getter property for DecimationFraction
   * @return Value of DecimationFraction
   */
  float getDecimationFraction() const;
  Q_PROPERTY(float DecimationFraction READ getDecimationFraction WRITE setDecimationFraction)

  /**
   * @brief Setter property for GridResolution
   */
  void setGridResolution(const FloatVec3Type& value);

  /**
   * @brief Getter property for GridResolution
   * @return Value of GridResolution
   */
  FloatVec3Type getGridResolution() const;
  Q_PROPERTY(FloatVec3Type GridResolution READ getGridResolution WRITE setGridResolution)

  /**
   * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
   */
  QString getCompiledLibraryName() const override;

  /**
   * @brief getBrandingString Returns the branding string for the filter, which is a tag
   * used to denote the filter's association with specific plugins
   * @return Branding string
   */
  QString getBrandingString() const override;

  /**
   * @brief getFilterVersion Returns a version string for this filter. Default
   * value is an empty string.
   * @return
   */
  QString getFilterVersion() const override;

  /**
   * @brief newFilterInstance Reimplemented from @see AbstractFilter class
   */
  AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters) const override;

  /**
   * @brief getGroupName Reimplemented from @see AbstractFilter class
   */
  QString getGroupName() const override;

  /**
   * @brief getSubGroupName Reimplemented from @see AbstractFilter class
   */
  QString getSubGroupName() const override;

  /**
   * @brief getHumanLabel Reimplemented from @see AbstractFilter class
   */
  QString getHumanLabel() const override;

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
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

protected:
  DownsampleVertexGeometry();

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  DataArrayPath m_VertexAttrMatPath = {"", "", ""};
  int m_DownsampleType = {0};
  int m_DecimationFreq = {2};
  float m_DecimationFraction = {0.5};
  FloatVec3Type m_GridResolution = {1.0f, 1.0f, 1.0f};

  void removeNthPoint();

  void removeFractionPoints();

  void gridDownsample();

public:
  DownsampleVertexGeometry(const DownsampleVertexGeometry&) = delete;            // Copy Constructor Not Implemented
  DownsampleVertexGeometry(DownsampleVertexGeometry&&) = delete;                 // Move Constructor Not Implemented
  DownsampleVertexGeometry& operator=(const DownsampleVertexGeometry&) = delete; // Copy Assignment Not Implemented
  DownsampleVertexGeometry& operator=(DownsampleVertexGeometry&&) = delete;      // Move Assignment Not Implemented
};
