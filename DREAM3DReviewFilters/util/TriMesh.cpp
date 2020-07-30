#include "TriMesh.h"

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TriMesh::TriMesh(VertexCoordList vertices)
{
  initialize(vertices);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TriMesh::~TriMesh() = default;

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TriMesh::initialize(VertexCoordList& vertices)
{
  m_Triangles = TriList();
  m_Vertices = VertexList();
  for(auto&& vert : vertices)
  {
    m_Vertices.emplace_back(TriMeshPrimitives::Vertex(vert));
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TriMesh::buildTriangleLinks()
{
  for(auto&& vert : m_Vertices)
  {
    vert.triangleLinks.clear();
  }

  for(auto t = 0; t < m_Triangles.size(); t++)
  {
    addLinkToTriangle(m_Triangles[t].vert0, t);
    addLinkToTriangle(m_Triangles[t].vert1, t);
    addLinkToTriangle(m_Triangles[t].vert2, t);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TriMesh::addLinkToTriangle(int64_t vertex, int64_t triangle)
{
  m_Vertices[vertex].triangleLinks.push_back(triangle);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TriMesh::removeLinkFromTriangle(int64_t vertex, int64_t triangle)
{
  m_Vertices[vertex].triangleLinks.remove(triangle);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
std::vector<int64_t> TriMesh::getTrianglesToVertex(int64_t vertex)
{
  std::vector<int64_t> triangles;
  triangles.reserve(m_Vertices[vertex].triangleLinks.size());
  std::copy(std::begin(m_Vertices[vertex].triangleLinks), std::end(m_Vertices[vertex].triangleLinks), std::back_inserter(triangles));

  return triangles;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int64_t TriMesh::getTriangleEdgeNeighbor(int64_t vertex0, int64_t vertex1, int64_t triangle)
{
  std::list<int64_t> vert0tris = m_Vertices[vertex0].triangleLinks;
  std::list<int64_t> vert1tris = m_Vertices[vertex1].triangleLinks;
  vert0tris.sort();
  vert1tris.sort();

  std::list<int64_t> sharedTriangles;
  std::set_intersection(vert0tris.begin(), vert0tris.end(), vert1tris.begin(), vert1tris.end(), std::back_inserter(sharedTriangles));

  if(sharedTriangles.size() > 2)
  {
    return -1;
  }                                 // non-manifold
  sharedTriangles.remove(triangle); // remove the source triangle
  if(sharedTriangles.empty())
  {
    return -1;
  } // boundary

  return (*sharedTriangles.begin());
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TriMesh::replaceTriangleVertices(int64_t vertex0, int64_t vertex1, int64_t vertex2, int64_t triangle)
{
  m_Triangles[triangle].vert0 = vertex0;
  m_Triangles[triangle].vert1 = vertex1;
  m_Triangles[triangle].vert2 = vertex2;
  m_Triangles[triangle].updateEdges();
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int64_t TriMesh::addTriangle(int64_t vertex0, int64_t vertex1, int64_t vertex2)
{
  m_Triangles.emplace_back(TriMeshPrimitives::Triangle(vertex0, vertex1, vertex2));
  addLinkToTriangle(vertex0, m_Triangles.size() - 1);
  addLinkToTriangle(vertex1, m_Triangles.size() - 1);
  addLinkToTriangle(vertex2, m_Triangles.size() - 1);

  return (m_Triangles.size() - 1);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TriMesh::removeVertices(std::vector<int64_t> vertices)
{
  for(auto&& index : vertices)
  {
    VertexList::iterator iter = m_Vertices.begin() + index;
    m_Vertices.erase(iter);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TriMesh::removeTriangles(std::vector<int64_t> triangles)
{
  for(auto&& index : triangles)
  {
    TriList::iterator iter = m_Triangles.begin() + index;
    m_Triangles.erase(iter);
  }
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TriMesh::getTriangleVertices(int64_t triangle, int64_t vertices[3])
{
  vertices[0] = m_Triangles[triangle].vert0;
  vertices[1] = m_Triangles[triangle].vert1;
  vertices[2] = m_Triangles[triangle].vert2;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TriMesh::getVertexCoordinates(int64_t vertex, float coordinates[3])
{
  coordinates[0] = m_Vertices[vertex].vert[0];
  coordinates[1] = m_Vertices[vertex].vert[1];
  coordinates[2] = m_Vertices[vertex].vert[2];
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
void TriMesh::getVertexCoordinates(int64_t vertex, double coordinates[3])
{
  coordinates[0] = static_cast<double>(m_Vertices[vertex].vert[0]);
  coordinates[1] = static_cast<double>(m_Vertices[vertex].vert[1]);
  coordinates[2] = static_cast<double>(m_Vertices[vertex].vert[2]);
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int64_t TriMesh::getOppositeVertex(int64_t vertex0, int64_t vertex1, int64_t triangle)
{
  TriMeshPrimitives::Edge edge(vertex0, vertex1);
  return (m_Triangles[triangle].getOppositeVertex(edge));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
int64_t TriMesh::getNumberOfTriangles()
{
  return (static_cast<int64_t>(m_Triangles.size()));
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TriangleGeom::Pointer TriMesh::convertToTriangleGeometry()
{
  SharedVertexList::Pointer vertices = TriangleGeom::CreateSharedVertexList(m_Vertices.size());
  TriangleGeom::Pointer triangles = TriangleGeom::CreateGeometry(m_Triangles.size(), vertices, SIMPL::Geometry::TriangleGeometry);
  float* verts = triangles->getVertexPointer(0);
  size_t* tris = triangles->getTriPointer(0);

  for(auto i = 0; i < m_Vertices.size(); i++)
  {
    verts[3 * i + 0] = m_Vertices[i].vert[0];
    verts[3 * i + 1] = m_Vertices[i].vert[1];
    verts[3 * i + 2] = m_Vertices[i].vert[2];
  }

  for(auto i = 0; i < m_Triangles.size(); i++)
  {
    tris[3 * i + 0] = m_Triangles[i].vert0;
    tris[3 * i + 1] = m_Triangles[i].vert1;
    tris[3 * i + 2] = m_Triangles[i].vert2;
  }

  return triangles;
}

// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------
TriMesh::Pointer TriMesh::NullPointer()
{
  return Pointer(static_cast<Self*>(nullptr));
}
