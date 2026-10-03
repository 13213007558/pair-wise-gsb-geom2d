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
  * Clipped Power Diagram Area Optimization
  * Boolean Operations on PolygonSets

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
~~~

Power diagram weights can be fitted so that the cells clipped to a bounded polygon match prescribed
areas. The returned weights have mean zero, which removes the constant-weight gauge freedom:

~~~ ruby
domain = Geom2D::Polygon([0, 0], [1, 0], [1, 1], [0, 1])
sites = [[0.25, 0.25], [0.75, 0.25], [0.75, 0.75], [0.25, 0.75]]
result = Geom2D::Algorithms::PowerDiagram.solve(domain, sites, [0.7, 0.1, 0.1, 0.1])

result.status             # => :optimal
result.weights            # => mean-zero power weights
result.areas              # => clipped cell areas
result.area_residuals     # => achieved area - target area, per site
result.coverage_residual  # => absolute total-area mismatch
~~~
