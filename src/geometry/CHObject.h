#pragma once

/// The convex hull of a set of vertices.
/// Points that lie in its interior are ignored,
/// but should be avoided for optimal performance.
struct CHObject{
  /// The number of vertices comprising the convex hull of this object.
  int num_points;
  const float (* points)[3];
};
