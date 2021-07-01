#pragma once

#include <stdint.h>
#include "LinAlg.h"

/// Face-Vertex Triangle Mesh (one array for vertices, one for indexing)
/// Intended only for static geometry: Cannot easily modify the TriMesh once instantiated
struct TriMeshI{
  /// The number of vertices in this mesh
  int num_vertices;
  /// The number of triangles in this mesh
  int num_faces;
  /// The vertices of this mesh
  float (* vs)[3];
  /// The indices for each triangle (counter clock wise)
  int (* idx)[3];
};

/// AoS Face-Vertex Triangle Mesh: directly store vertices of the triangles
/// Intended only for static geometry: Cannot easily modify the TriMesh once instantiated
struct TriMesh {
  /// The number of triangles in this mesh
  int num_faces;
  /// The indices for each triangle (counter clock wise)
  struct Vec3f (* const faces)[3];
};
