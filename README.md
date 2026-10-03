# Geom2D - Objects and Algorithms for 2D Geometry in Ruby

This library implements objects for 2D geometry, like points, lines, line segments, arcs, curves and
so on, as well as algorithms for these objects, like line-line intersections and arc approximation
by Bézier curves.


## License

Copyright (C) 2018-2023 Thomas Leitner <t_leitner@gmx.at>, licensed under the MIT - see the
**LICENSE** file.


## Features

* Objects
  * Point
  * Segment
  * Polygon
  * PolygonSet
  * Rectangle
  * Polyline (TODO)
  * Rectangle (TODO)
  * QuadraticCurve (TODO)
  * QubicCurve (TODO)
  * Arc (TODO)
  * Circle (TODO)
  * Path (TODO)
* Algorithms
  * Segment-Segment Intersection
  * Boolean Operations on PolygonSets
  * Constrained Delaunay Triangulation of PolygonSets (multiple outer contours and holes)

## Usage

~~~ ruby
require 'geom2d'

# Point, can also be interpreted as vector
point1 = Geom2D::Point(2, 2)
point2 = Geom2D::Point([2, 2])   # arrays are fine but not as efficient
point3 = Geom2D::Point(point2)   # copy constructor

# Segment defined by two points or a point and a vector
line1 = Geom2D::Segment(point1, point2)
line2 = Geom2D::Segment(point1, vector: point2)
line3 = Geom2D::Segment([3, 4], [9, 6])   # arrays are also possible

# Segment intersection
line1.intersect(line3)  # => intersection_point

# Constrained Delaunay triangulation: outer contours counterclockwise, holes clockwise
outer = Geom2D::Polygon([0, 0], [10, 0], [10, 10], [0, 10])
hole  = Geom2D::Polygon([3, 3], [3, 7], [7, 7], [7, 3])
mesh = Geom2D::Algorithms::ConstrainedDelaunayTriangulation.run(
  Geom2D::PolygonSet(outer, hole)
)
mesh.triangles.each do |triangle|
  # triangle.vertices  -> indices into mesh.points
  # triangle.neighbors -> indices into mesh.triangles (nil on the boundary)
end
mesh.area                 # => total area of the triangulated region
mesh.constraint_edges     # => all preserved constraint edges

# The result of the boolean operations can directly be triangulated
subject = Geom2D::PolygonSet(Geom2D::Polygon([0, 0], [10, 0], [10, 10], [0, 10]))
clipping = Geom2D::PolygonSet(Geom2D::Polygon([5, 5], [15, 5], [15, 15], [5, 15]))
result = Geom2D::Algorithms::PolygonOperation.run(subject, clipping, :difference)
Geom2D::Algorithms::ConstrainedDelaunayTriangulation.run(result).triangles
~~~
