# Based on CreateQuadGeometry example
# Tests the K Distance Graph filter

# These are the simpl_py python modules
import os
import simpl
import simplpy
import simpl_helpers as sh
import simpl_test_dirs as sd
import dream3dreviewpy

def start_test():
    # Create Data Container Array
    dca = simpl.DataContainerArray()

    # Create the Data Container
    err = simplpy.create_data_container(dca, 'DataContainer')
    assert err == 0, f'DataContainer ErrorCondition: {err}'

    # Read vertices
    err = simplpy.create_attribute_matrix(dca, simpl.DataArrayPath('DataContainer', 'Bounds', ''), simpl.AttributeMatrix.Type.Cell, sh.CreateDynamicTableData([[144]]))
    assert err == 0, f'CreateAttributeMatrix - Error: {err}'

    err = simplpy.import_asci_data_array(dca, simpl.DataArrayPath('DataContainer', 'Bounds', 'Vertices'), simpl.NumericTypes.Float, 3, 1, sd.GetBuildDirectory() + '/Data/SIMPL/VertexCoordinates.csv', 0)
    assert err == 0, f'ImportAsciDataArray - Error {err}'

    # Read quads
    err = simplpy.create_attribute_matrix(dca, simpl.DataArrayPath('DataContainer', 'QuadList', ''), simpl.AttributeMatrix.Type.Cell, sh.CreateDynamicTableData([[121]]))
    assert err == 0, f'CreateAttributeMatrix - Error: {err}'

    err = simplpy.import_asci_data_array(dca, simpl.DataArrayPath('DataContainer', 'QuadList', 'Quads'), simpl.NumericTypes.SizeT, 4, 1, sd.GetBuildDirectory() + '/Data/SIMPL/QuadConnectivity.csv', 0)
    assert err == 0, f'ImportAsciDataArray - Error {err}'

    # Create Geometry
    err = sh.CreateGeometry(dca, 0, simpl.IGeometry.Type.Quad, 'DataContainer', False,
                            shared_vertex_list_array_path=simpl.DataArrayPath('DataContainer', 'Bounds', 'Vertices'),
                            shared_quad_list_array_path=simpl.DataArrayPath('DataContainer', 'QuadList', 'Quads'),
                            vertex_attribute_matrix_name='VertexData',
                            face_attribute_matrix_name='FaceData')
    assert err == 0, f'Create Geometry -  ErrorCondition: {err}'

    # K Distance Graph
    err = dream3dreviewpy.k_distance_graph(dca, simpl.DataArrayPath('DataContainer', 'QuadList',
                                                                    'Quads'),
                                           False, simpl.DataArrayPath('', '', ''),
                                           simpl.DataArrayPath('DataContainer', 'QuadList',
                                                               'KDistance'), 7, 3)
    assert err == 0, f'KDistanceGraph ErrorCondition: {err}'

    # Write DREAM3D File
    err = simplpy.data_container_writer(dca, sd.GetBuildDirectory() +
                                        '/Data/Output/DREAM3DReview/' +
                                        'KDistanceGraph.dream3d', True, False)
    assert err == 0, f'DataContainerWriter ErrorCondition: {err}'

if __name__ == '__main__':
    print('Starting Test %s ' % os.path.basename(__file__))
    start_test()
    print('Ending Test %s ' % os.path.basename(__file__))