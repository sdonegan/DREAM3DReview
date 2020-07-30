#pragma once

#include "SIMPLib/Geometry/TriangleGeom.h"

#include "DREAM3DReviewFilters/util/TriMeshPrimitives.hpp"

/**
 * @brief The TriMesh class
 */
class TriMesh
{
public:
  using Self = TriMesh;
  using Pointer = std::shared_ptr<Self>;
  using ConstPointer = std::shared_ptr<const Self>;
  using WeakPointer = std::weak_ptr<Self>;
  using ConstWeakPointer = std::weak_ptr<const Self>;
  static Pointer NullPointer();

  virtual ~TriMesh();

  typedef std::vector<std::vector<float>> VertexCoordList;
  typedef std::vector<TriMeshPrimitives::Vertex> VertexList;
  typedef std::vector<TriMeshPrimitives::Triangle> TriList;

  static Pointer New(VertexCoordList vertices)
  {
    Pointer sharedPtr(new TriMesh(vertices));
    return sharedPtr;
  }

  VertexList getVertices()
  {
    return m_Vertices;
  }

  TriList getTriangles()
  {
    return m_Triangles;
  }

  void buildTriangleLinks();

  void addLinkToTriangle(int64_t vertex, int64_t triangle);

  void removeLinkFromTriangle(int64_t vertex, int64_t triangle);

  std::vector<int64_t> getTrianglesToVertex(int64_t vertex);

  int64_t getTriangleEdgeNeighbor(int64_t vertex0, int64_t vertex1, int64_t triangle);

  void replaceTriangleVertices(int64_t vertex0, int64_t vertex1, int64_t vertex2, int64_t triangle);

  int64_t addTriangle(int64_t vertex0, int64_t vertex1, int64_t vertex2);

  void removeVertices(std::vector<int64_t> vertices);

  void removeTriangles(std::vector<int64_t> triangles);

  void getTriangleVertices(int64_t triangle, int64_t vertices[3]);

  void getVertexCoordinates(int64_t vertex, float coordinates[3]);

  void getVertexCoordinates(int64_t vertex, double coordinates[3]);

  int64_t getOppositeVertex(int64_t vertex0, int64_t vertex1, int64_t triangle);

  int64_t getNumberOfTriangles();

  TriangleGeom::Pointer convertToTriangleGeometry();

protected:
  explicit TriMesh(VertexCoordList vertices);

private:
  VertexList m_Vertices;
  TriList m_Triangles;

  void initialize(VertexCoordList& vertices);

  TriMesh(const TriMesh&) = delete;        // Copy Constructor Not Implemented
  void operator=(const TriMesh&) = delete; // Operator '=' Not Implemented
};
