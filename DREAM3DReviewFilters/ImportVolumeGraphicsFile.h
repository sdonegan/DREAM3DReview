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

#include <memory>

#include <QtCore/QFile>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The ImportVolumeGraphicsFile class. See [Filter documentation](@ref ImportVolumeGraphicsFile) for details.
 */
class DREAM3DReview_EXPORT ImportVolumeGraphicsFile : public AbstractFilter
{
  Q_OBJECT
  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(ImportVolumeGraphicsFile SUPERCLASS AbstractFilter)
  PYB11_SHARED_POINTERS(ImportVolumeGraphicsFile)
  PYB11_FILTER_NEW_MACRO(ImportVolumeGraphicsFile)
  PYB11_PROPERTY(QString VGHeaderFile READ getVGHeaderFile WRITE setVGHeaderFile)
  PYB11_PROPERTY(QString DataContainerName READ getDataContainerName WRITE setDataContainerName)
  PYB11_PROPERTY(QString CellAttributeMatrixName READ getCellAttributeMatrixName WRITE setCellAttributeMatrixName)
  PYB11_PROPERTY(QString DensityArrayName READ getDensityArrayName WRITE setDensityArrayName)
  PYB11_METHOD(QString getVGDataFile)
  PYB11_END_BINDINGS()
  // End Python bindings declarations
public:
  using Self = ImportVolumeGraphicsFile;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ImportVolumeGraphicsFile> New();

  /**
   * @brief Returns the name of the class for ImportVolumeGraphicsFile
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ImportVolumeGraphicsFile
   */
  static QString ClassName();

  ~ImportVolumeGraphicsFile() override;

  struct RepresentationBlock
  {
    SizeVec3Type size = {0, 0, 0};
    std::string resamplemode;
    FloatVec4Type mirror = {0.0F, 0.0F, 0.0F, 0.0F};
    std::string dataType;
    FloatVec2Type dataRange = {0.0F, 0.0F};
    int8_t bitsPerElement = {0};
  };
  struct FileBlock
  {
    SizeVec3Type RegionOfInterestStart = {0, 0, 0};
    SizeVec3Type RegionOfInterestEnd = {0, 0, 0};
    std::string FileFormat = "raw";
    SizeVec3Type Size = {0, 0, 0};
    std::string Name = "AS_N610_02.vol";
    std::string Datatype = "float";
    FloatVec2Type datarange = {0.0F, 0.0F};
    int8_t BitsPerElement = 32;
  };
  struct GeometryBlock
  {
    std::string status = "visible";
    FloatVec3Type relativeposition = {0.0F, 0.0F, 0.0F};
    FloatVec3Type position = {0.0F, 0.0F, 0.0F};
    FloatVec3Type resolution = {0.0F, 0.0F, 0.0F};
    FloatVec3Type scale = {1.0F, 1.0F, 1.0F};
    FloatVec3Type center = {0.0F, 0.0F, 0.0F};
    std::array<int32_t, 9> rotate = {1, 0, 0, 0, 1, 0, 0, 0, 1};
    std::string unit = "mm";
  };

  /**
   * @brief Setter property for VGDataFile
   */
  void setVGDataFile(const QString& value);
  /**
   * @brief Getter property for VGDataFile
   * @return Value of VGDataFile
   */
  QString getVGDataFile() const;
  Q_PROPERTY(QString VGDataFile READ getVGDataFile WRITE setVGDataFile)

  /**
   * @brief Setter property for VGHeaderFile
   */
  void setVGHeaderFile(const QString& value);
  /**
   * @brief Getter property for VGHeaderFile
   * @return Value of VGHeaderFile
   */
  QString getVGHeaderFile() const;
  Q_PROPERTY(QString VGHeaderFile READ getVGHeaderFile WRITE setVGHeaderFile)

  /**
   * @brief Setter property for DataContainerName
   */
  void setDataContainerName(const QString& value);
  /**
   * @brief Getter property for DataContainerName
   * @return Value of DataContainerName
   */
  QString getDataContainerName() const;
  Q_PROPERTY(QString DataContainerName READ getDataContainerName WRITE setDataContainerName)

  /**
   * @brief Setter property for CellAttributeMatrixName
   */
  void setCellAttributeMatrixName(const QString& value);
  /**
   * @brief Getter property for CellAttributeMatrixName
   * @return Value of CellAttributeMatrixName
   */
  QString getCellAttributeMatrixName() const;
  Q_PROPERTY(QString CellAttributeMatrixName READ getCellAttributeMatrixName WRITE setCellAttributeMatrixName)

  /**
   * @brief Setter property for DensityArrayName
   */
  void setDensityArrayName(const QString& value);
  /**
   * @brief Getter property for DensityArrayName
   * @return Value of DensityArrayName
   */
  QString getDensityArrayName() const;
  Q_PROPERTY(QString DensityArrayName READ getDensityArrayName WRITE setDensityArrayName)

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
   * @brief execute Reimplemented from @see AbstractFilter class
   */
  void execute() override;

  /**
   * @brief getUuid Return the unique identifier for this filter.
   * @return A QUuid object.
   */
  QUuid getUuid() const override;

protected:
  ImportVolumeGraphicsFile();

  /**
   * @brief readBinaryCTFile Reads the raw binary CT file
   * @return Integer error code
   */
  int32_t readVolFile();

  /**
   * @brief readHeaderMetaData Reads the number of voxels and voxel extents
   * from the NSI header file
   * @return Integer error code
   */
  void readHeaderMetaData(ImageGeom::Pointer image);

  /**
   * @brief parseDimensions
   * @param in
   * @return
   */
  SizeVec3Type parseDimensions(QFile& in);

  /**
   * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
   */
  void dataCheck() override;

  /**
   * @brief Initializes all the private instance variables.
   */
  void initialize();

private:
  std::weak_ptr<DataArray<float>> m_DensityPtr;
  float* m_Density = nullptr;

  QString m_VGDataFile = {};
  QString m_VGHeaderFile = {};
  QString m_DataContainerName = {"VolumeGraphics"};
  QString m_CellAttributeMatrixName = {"CT Data"};
  QString m_DensityArrayName = {"Density"};

  QFile m_InHeaderStream;
  QFile m_InStream;

public:
  ImportVolumeGraphicsFile(const ImportVolumeGraphicsFile&) = delete;            // Copy Constructor Not Implemented
  ImportVolumeGraphicsFile(ImportVolumeGraphicsFile&&) = delete;                 // Move Constructor Not Implemented
  ImportVolumeGraphicsFile& operator=(const ImportVolumeGraphicsFile&) = delete; // Copy Assignment Not Implemented
  ImportVolumeGraphicsFile& operator=(ImportVolumeGraphicsFile&&) = delete;      // Move Assignment Not Implemented
};
