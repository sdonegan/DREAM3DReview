# Based on (03) SmallIN100 Mesh Statistics (from EBSD Surface Meshing)
# Tests the Normalize Attribute Arrays filter
#
import os
import simpl
import simplpy
import simpl_helpers as sh
import simpl_test_dirs as sd
import orientationanalysispy
import surfacemeshingpy
import dream3dreviewpy

def start_test():
    # Create Data Container Array
    dca = simpl.DataContainerArray()

    # Read DREAM3D File
    dcap = simpl.DataContainerArrayProxy()
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'Confidence Index').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'CriticalField').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'EulerAngles').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'FeatureIds').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'FeatureReferenceMisorientations').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'Fit').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'GBManhattanDistances').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'IPFColor').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'Image Quality').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'KernelAverageMisorientations').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'Mask').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'ParentIds').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'Phases').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'QPManhattanDistances').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'Quats').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'SEM Signal').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('EBSD Scan Data').getDataArrayProxy(
        'TJManhattanDistances').Flag = 2

    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Phase Data').getDataArrayProxy(
        'CrystalStructures').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Phase Data').getDataArrayProxy(
        'LatticeConstants').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Phase Data').getDataArrayProxy(
        'MaterialName').Flag = 2

    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy('Active').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'AspectRatios').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'AvgEuler').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'AvgEulerAngles').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'AvgQuats').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'AxisEulerAngles').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'AxisLengths').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'Centroids').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'CriticalFields').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'EquivalentDiameters').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'FeatureAvgMisorientations').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'MisorientationList').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'NeighborList').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'NeighborhoodList').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'Neighborhoods').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'NumElements').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'NumNeighbors').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'NumNeighbors2').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'Omega3s').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'ParentIds').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy('Poles').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy('Phases').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'Schmids').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'Shape Volumes').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'SharedSurfaceAreaList').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'Size Volumes').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'SlipSystems').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'Sphericity').Flag = 2
    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('Grain Data').getDataArrayProxy(
        'SurfaceAreaVolumeRatio').Flag = 2

    dcap.getDataContainerProxy('Small IN100').getAttributeMatrixProxy('NewGrain Data').getDataArrayProxy(
        'Active').Flag = 2

    dcap.getDataContainerProxy('TriangleDataContainer').getAttributeMatrixProxy('FaceData').getDataArrayProxy(
        'FaceLabels').Flag = 2
    dcap.getDataContainerProxy('TriangleDataContainer').getAttributeMatrixProxy('FaceFeatureData').Flag = 2
    dcap.getDataContainerProxy('TriangleDataContainer').getAttributeMatrixProxy('VertexData').getDataArrayProxy(
        'NodeType').Flag = 2

    err = simplpy.data_container_reader(dca,
                                        sd.GetBuildDirectory() +
                                        '/Data/Output/SurfaceMesh/SmallIN100_Smoothed.dream3d',
                                        False, dcap)
    assert err == 0, f'DataContainerReader ErrorCondition {err}'

    # Generate Triangle Areas
    err = surfacemeshingpy.triangle_area_filter(dca,
                                                simpl.DataArrayPath('TriangleDataContainer', 'FaceData', 'FaceAreas'))
    assert err == 0, f'TriangleAreaFilter ErrorCondition {err}'

    # Generate Triangle Normals
    err = surfacemeshingpy.triangle_normal_filter(dca,
                                                  simpl.DataArrayPath('TriangleDataContainer', 'FaceData',
                                                                      'FaceNormals'))
    assert err == 0, f'TriangleNormalFilter ErrorCondition {err}'

    # Find Minimum Triangle Dihedral Angle
    err = surfacemeshingpy.triangle_dihedral_angle_filter(dca, simpl.DataArrayPath('TriangleDataContainer', 'FaceData',
                                                                                   'FaceDihedralAngles'))
    assert err == 0, f'TriangleDihedralAngleFilter ErrorCondition {err}'

    # Generate IPF Colors (Face)
    err = orientationanalysispy.generate_face_ipf_coloring(dca,
                                                           simpl.DataArrayPath('TriangleDataContainer', 'FaceData',
                                                                               'FaceLabels'),
                                                           simpl.DataArrayPath('TriangleDataContainer', 'FaceData',
                                                                               'FaceNormals'),
                                                           simpl.DataArrayPath('Small IN100', 'Grain Data',
                                                                               'AvgEulerAngles'),
                                                           simpl.DataArrayPath('Small IN100', 'Grain Data', 'Phases'),
                                                           simpl.DataArrayPath('Small IN100', 'Phase Data',
                                                                               'CrystalStructures'),
                                                           'SurfaceMeshFaceIPFColors')
    assert err == 0, f'GenerateFaceIPFColoring ErrorCondition {err}'

    # Generate Misorientation Colors (Face)
    err = orientationanalysispy.generate_face_misorientation_coloring(dca,
                                                                      simpl.DataArrayPath('TriangleDataContainer',
                                                                                          'FaceData', 'FaceLabels'),
                                                                      simpl.DataArrayPath('Small IN100', 'Grain Data',
                                                                                          'AvgQuats'),
                                                                      simpl.DataArrayPath('Small IN100', 'Grain Data',
                                                                                          'Phases'),
                                                                      simpl.DataArrayPath('Small IN100', 'Phase Data',
                                                                                          'CrystalStructures'),
                                                                      'SurfaceMeshFaceMisorientationColors')
    assert err == 0, f'GenerateFaceMisorientationColoring ErrorCondition {err}'

    # Normalize Attribute Arrays
    selected_data_arrays = [simpl.DataArrayPath('TriangleDataContainer', 'FaceData', 'FaceDihedralAngles'),
                            simpl.DataArrayPath('TriangleDataContainer', 'FaceData', 'FaceAreas')]
    err = dream3dreviewpy.normalize_arrays(dca, selected_data_arrays, 0, 0, 1, '_Normalized', False,
                                           simpl.DataArrayPath('', '', ''), 0)
    assert err == 0, f'NormalizeArrays ErrorCondition {err}'

    # Write to DREAM3D file
    err = sh.WriteDREAM3DFile(sd.GetBuildDirectory() + '/Data/Output/DREAM3DReview/NormalizeArrays.dream3d',
                              dca)
    assert err == 0, f'WriteDREAM3DFile ErrorCondition: {err}'

if __name__ == '__main__':
    print('Starting Test %s ' % os.path.basename(__file__))
    start_test()
    print('Ending Test %s ' % os.path.basename(__file__))