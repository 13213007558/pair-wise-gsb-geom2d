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

### Guillotine cutting

For small rectangular cutting-stock instances, the exact guillotine optimizer returns a binary
cutting tree rather than just piece coordinates. Every internal cut runs completely through its
current sub-board, and the integer saw kerf is included in the panel geometry. Per-piece rotation
is opt-in and unused panels may be left as waste:

~~~ ruby
orders = [
  { width: 10, height: 4, quantity: 3, value: 12, rotatable: true },
  { width: 4, height: 4, quantity: 2, value: 8 }
]
solution = Geom2D::Algorithms::GuillotinePacking.solve(
  24, 12, orders, kerf_width: 1, time_limit: 5, node_limit: 10_000
)

solution.total_value       # total value of ordered leaves
solution.optimal?          # false if the time or node budget stopped exhaustive proof
solution.used_quantities   # counts for every order
solution.validate!         # checks geometry, quantities, kerfs, and leaf value sums

solution.cuts.each do |cut|
  # cut.direction is :vertical or :horizontal
  # cut.position is the width or height of the first child
  # cut.first_child and cut.second_child form the executable binary tree
end
~~~
