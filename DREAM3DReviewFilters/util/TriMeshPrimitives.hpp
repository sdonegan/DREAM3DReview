#pragma once

#include <list>
#include <vector>

#include "SIMPLib/SIMPLib.h"

namespace TriMeshPrimitives
{
class Vertex
{
public:
  Vertex(const std::vector<float>& vert)
  : vert(vert)
  {
  }

  std::vector<float> vert;
  std::list<int64_t> triangleLinks;
};

inline bool operator==(const Vertex& vert0, const Vertex& vert1)
{
  return (vert0.vert == vert1.vert);
}

class Edge
{
public:
  Edge(const int64_t& vert0, const int64_t& vert1)
  : vert0(vert0)
  , vert1(vert1)
  {
  }

  int64_t vert0;
  int64_t vert1;
};

inline bool operator==(const Edge& edge0, const Edge& edge1)
{
  return (edge0.vert0 == edge1.vert0 && edge0.vert1 == edge1.vert1) || (edge0.vert0 == edge1.vert1 && edge0.vert1 == edge1.vert0);
}

class Triangle
{
public:
  Triangle(const int64_t& vert0, const int64_t& vert1, const int64_t& vert2)
  : vert0(vert0)
  , vert1(vert1)
  , vert2(vert2)
  , edge0(vert0, vert1)
  , edge1(vert1, vert2)
  , edge2(vert2, vert0)
  {
  }

  int64_t vert0;
  int64_t vert1;
  int64_t vert2;
  Edge edge0;
  Edge edge1;
  Edge edge2;

  int64_t getOppositeVertex(Edge& edge)
  {
    if(edge == edge0)
    {
      return vert2;
    }
    else if(edge == edge1)
    {
      return vert0;
    }
    else if(edge == edge2)
    {
      return vert1;
    }
    else
    {
      return -1;
    }
  }

  void updateEdges()
  {
    edge0 = Edge(vert0, vert1);
    edge1 = Edge(vert1, vert2);
    edge2 = Edge(vert2, vert0);
  }
};

inline bool operator==(const Triangle& tri0, const Triangle& tri1)
{
  return (tri0.vert0 == tri1.vert0 || tri0.vert0 == tri1.vert1 || tri0.vert0 == tri1.vert2) && (tri0.vert1 == tri1.vert0 || tri0.vert1 == tri1.vert1 || tri0.vert1 == tri1.vert2) &&
         (tri0.vert2 == tri1.vert0 || tri0.vert2 == tri1.vert1 || tri0.vert2 == tri1.vert2);
}
} // namespace TriMeshPrimitives
