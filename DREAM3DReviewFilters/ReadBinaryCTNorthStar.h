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
/* Revisions past APRIL 5 2019 use the below license */
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
#pragma once

#include <memory>

#include <QtCore/QFile>

#include "SIMPLib/Filtering/AbstractFilter.h"
#include "SIMPLib/Geometry/ImageGeom.h"
#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/DataArrays/DataArray.hpp"

#include "DREAM3DReview/DREAM3DReviewDLLExport.h"

/**
 * @brief The ReadBinaryCTNorthStar class. See [Filter documentation](@ref readbinaryctnorthstar) for details.
 */
class DREAM3DReview_EXPORT ReadBinaryCTNorthStar : public AbstractFilter
{
  Q_OBJECT

  // clang-format off
  // Start Python bindings declarations
  PYB11_BEGIN_BINDINGS(ReadBinaryCTNorthStar SUPERCLASS AbstractFilter)
  PYB11_FILTER()
  PYB11_SHARED_POINTERS(ReadBinaryCTNorthStar)
  PYB11_FILTER_NEW_MACRO(ReadBinaryCTNorthStar)
  PYB11_PROPERTY(QString InputHeaderFile READ getInputHeaderFile WRITE setInputHeaderFile)
  PYB11_PROPERTY(QString DataContainerName READ getDataContainerName WRITE setDataContainerName)
  PYB11_PROPERTY(QString CellAttributeMatrixName READ getCellAttributeMatrixName WRITE setCellAttributeMatrixName)
  PYB11_PROPERTY(QString DensityArrayName READ getDensityArrayName WRITE setDensityArrayName)
  PYB11_PROPERTY(int32_t LengthUnit READ getLengthUnit WRITE setLengthUnit)
  PYB11_PROPERTY(bool ImportSubvolume READ getImportSubvolume WRITE setImportSubvolume)
  PYB11_PROPERTY(IntVec3Type StartVoxelCoord READ getStartVoxelCoord WRITE setStartVoxelCoord)
  PYB11_PROPERTY(IntVec3Type EndVoxelCoord READ getEndVoxelCoord WRITE setEndVoxelCoord)
  PYB11_PROPERTY(QString VolumeDescription READ getVolumeDescription)
  PYB11_PROPERTY(QString DataFileInfo READ getDataFileInfo)
  PYB11_PROPERTY(QString ImportedVolumeDescription READ getImportedVolumeDescription)
  PYB11_END_BINDINGS()
  // End Python bindings declarations
  // clang-format on

public:
  using Self = ReadBinaryCTNorthStar;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  static std::shared_ptr<ReadBinaryCTNorthStar> New();

  /**
   * @brief Returns the name of the class for ReadBinaryCTNorthStar
   */
  QString getNameOfClass() const override;
  /**
   * @brief Returns the name of the class for ReadBinaryCTNorthStar
   */
  static QString ClassName();

  ~ReadBinaryCTNorthStar() override;

  /**
   * @brief Setter property for InputHeaderFile
   */
  void setInputHeaderFile(const QString& value);
  /**
   * @brief Getter property for InputHeaderFile
   * @return Value of InputHeaderFile
   */
  QString getInputHeaderFile() const;
  Q_PROPERTY(QString InputHeaderFile READ getInputHeaderFile WRITE setInputHeaderFile)

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
   * @brief Setter property for LengthUnit
   */
  void setLengthUnit(int32_t value);
  /**
   * @brief Getter property for LengthUnit
   * @return Value of LengthUnit
   */
  int32_t getLengthUnit() const;
  Q_PROPERTY(int32_t LengthUnit READ getLengthUnit WRITE setLengthUnit)

  /**
   * @brief Setter property for ImportSubvolume
   */
  void setImportSubvolume(bool value);
  /**
   * @brief Getter property for ImportSubvolume
   * @return Value of ImportSubvolume
   */
  bool getImportSubvolume() const;
  Q_PROPERTY(bool ImportSubvolume READ getImportSubvolume WRITE setImportSubvolume)

  /**
   * @brief Setter property for StartVoxelCoord
   */
  void setStartVoxelCoord(const IntVec3Type& value);
  /**
   * @brief Getter property for StartVoxelCoord
   * @return Value of StartVoxelCoord
   */
  IntVec3Type getStartVoxelCoord() const;
  Q_PROPERTY(IntVec3Type StartVoxelCoord READ getStartVoxelCoord WRITE setStartVoxelCoord)

  /**
   * @brief Setter property for EndVoxelCoord
   */
  void setEndVoxelCoord(const IntVec3Type& value);
  /**
   * @brief Getter property for EndVoxelCoord
   * @return Value of EndVoxelCoord
   */
  IntVec3Type getEndVoxelCoord() const;
  Q_PROPERTY(IntVec3Type EndVoxelCoord READ getEndVoxelCoord WRITE setEndVoxelCoord)

  /**
   * @brief getNewBoxDimensions
   * @return
   */
  QString getVolumeDescription();
  Q_PROPERTY(QString VolumeDescription READ getVolumeDescription)

  /**
   * @brief getNewBoxDimensions
   * @return
   */
  QString getDataFileInfo();
  Q_PROPERTY(QString DataFileInfo READ getDataFileInfo)

  /**
   * @brief getNewBoxDimensions
   * @return
   */
  QString getImportedVolumeDescription();
  Q_PROPERTY(QString ImportedVolumeDescription READ getImportedVolumeDescription)

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
  ReadBinaryCTNorthStar();

  /**
   * @brief sanityCheckFileSizeVersusAllocatedSize Ensures the allocated array and the raw
   * binary file have the same number of bytes
   * @param allocatedBytes Number of bytes allocated
   * @param fileSize Size of the raw binary file
   * @return Integer error code
   */
  int32_t sanityCheckFileSizeVersusAllocatedSize(size_t allocatedBytes, size_t fileSize);

  /**
   * @brief readBinaryCTFile Reads the raw binary CT file
   * @return Integer error code
   */
  int32_t readBinaryCTFiles();

  /**
   * @brief readHeaderMetaData Reads the number of voxels and voxel extents
   * from the NSI header file
   * @return Integer error code
   */
  int32_t readHeaderMetaData();
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

  bool m_ImportSubvolume = {false};
  IntVec3Type m_StartVoxelCoord = {0, 0, 0};
  IntVec3Type m_EndVoxelCoord = {1, 1, 1};

  std::vector<std::pair<QString, int64_t>> m_DataFiles;
  QString m_InputHeaderFile = {};
  QString m_DataContainerName = {"CT DataContainer"};
  QString m_CellAttributeMatrixName = {"CT Scan Data"};
  QString m_DensityArrayName = {"Density"};

  QFile m_InHeaderStream;
  QFile m_InStream;

  ImageGeom::Pointer m_OriginalVolume;
  ImageGeom::Pointer m_ImportedVolume;

  int32_t m_LengthUnit = {7}; // millimeter length units

public:
  ReadBinaryCTNorthStar(const ReadBinaryCTNorthStar&) = delete;            // Copy Constructor Not Implemented
  ReadBinaryCTNorthStar(ReadBinaryCTNorthStar&&) = delete;                 // Move Constructor Not Implemented
  ReadBinaryCTNorthStar& operator=(const ReadBinaryCTNorthStar&) = delete; // Copy Assignment Not Implemented
  ReadBinaryCTNorthStar& operator=(ReadBinaryCTNorthStar&&) = delete;      // Move Assignment Not Implemented
};
