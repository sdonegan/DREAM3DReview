/*
 * Your License or Copyright can go here
 */

#ifndef _importmassifdata_h_
#define _importmassifdata_h_

#include "SIMPLib/SIMPLib.h"
#include "SIMPLib/Common/AbstractFilter.h"
#include "SIMPLib/Common/SIMPLibSetGetMacros.h"

/**
 * @brief The ImportMASSIFData class. See [Filter documentation](@ref importmassifdata) for details.
 */
class ImportMASSIFData : public AbstractFilter
{
  Q_OBJECT

  public:
    SIMPL_SHARED_POINTERS(ImportMASSIFData)
    SIMPL_STATIC_NEW_MACRO(ImportMASSIFData)
    SIMPL_TYPE_MACRO_SUPER(ImportMASSIFData, AbstractFilter)

    virtual ~ImportMASSIFData();

    SIMPL_FILTER_PARAMETER(QString, MassifInputFilePath)
    Q_PROPERTY(QString MassifInputFilePath READ getMassifInputFilePath WRITE setMassifInputFilePath)

    SIMPL_FILTER_PARAMETER(QString, FilePrefix)
    Q_PROPERTY(QString FilePrefix READ getFilePrefix WRITE setFilePrefix)

    SIMPL_FILTER_PARAMETER(int, StepNumber)
    Q_PROPERTY(int StepNumber READ getStepNumber WRITE setStepNumber)

    /**
     * @brief getCompiledLibraryName Reimplemented from @see AbstractFilter class
     */
    virtual const QString getCompiledLibraryName();

    /**
     * @brief getBrandingString Returns the branding string for the filter, which is a tag
     * used to denote the filter's association with specific plugins
     * @return Branding string
    */
    virtual const QString getBrandingString();

    /**
     * @brief getFilterVersion Returns a version string for this filter. Default
     * value is an empty string.
     * @return
     */
    virtual const QString getFilterVersion();

    /**
     * @brief newFilterInstance Reimplemented from @see AbstractFilter class
     */
    virtual AbstractFilter::Pointer newFilterInstance(bool copyFilterParameters);

    /**
     * @brief getGroupName Reimplemented from @see AbstractFilter class
     */
    virtual const QString getGroupName();

    /**
     * @brief getSubGroupName Reimplemented from @see AbstractFilter class
     */
    virtual const QString getSubGroupName();

    /**
     * @brief getHumanLabel Reimplemented from @see AbstractFilter class
     */
    virtual const QString getHumanLabel();

    /**
     * @brief setupFilterParameters Reimplemented from @see AbstractFilter class
     */
    virtual void setupFilterParameters();

    /**
     * @brief execute Reimplemented from @see AbstractFilter class
     */
    virtual void execute();

    /**
    * @brief preflight Reimplemented from @see AbstractFilter class
    */
    virtual void preflight();

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
    ImportMASSIFData();

    /**
    * @brief dataCheck Checks for the appropriate parameter values and availability of arrays
    */
    void dataCheck();

    /**
    * @brief Initializes all the private instance variables.
    */
    void initialize();

  private:
    QString           m_PaddedStep = "";

    DEFINE_DATAARRAY_VARIABLE(float, DField)
    DEFINE_DATAARRAY_VARIABLE(float, EField)
    DEFINE_DATAARRAY_VARIABLE(float, SField)

    /**
     * @brief createDataArrayPaths
     * @return
     */
    QVector<DataArrayPath> createDataArrayPaths();

    /**
     * @brief createHDF5DatasetPaths
     * @return
     */
    QVector<QString> createHDF5DatasetPaths();

    /**
     * @brief generateIndexString
     * @param index
     * @param maxIndex
     * @return
     */
    QString generateIndexString(int index, int maxIndex);

    /**
     * @brief getDataContainerGeometry
     * @param tDims
     * @param origin
     * @param res
     */
    void getDataContainerGeometry(QVector<size_t> &tDims, QVector<float> &origin, QVector<float> &res);

    /**
     * @brief readIDataArray
     * @param gid
     * @param name
     * @param metaDataOnly
     * @return
     */
    IDataArray::Pointer readIDataArray(hid_t gid, const QString& name, QVector<size_t> geoDims, bool metaDataOnly);

    ImportMASSIFData(const ImportMASSIFData&); // Copy Constructor Not Implemented
    void operator=(const ImportMASSIFData&); // Operator '=' Not Implemented
};

#endif /* _ImportMASSIFData_H_ */
