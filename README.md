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
  * Static and Moving Convex Hulls
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

# Static convex hull; vertices are returned in counter-clockwise order.
points = [Geom2D.Point(0, 0), Geom2D.Point(2, 0), Geom2D.Point(1, 1)]
Geom2D::Algorithms.convex_hull(points)

# Event-driven analysis of points moving with constant velocity.
MovingPoint = Geom2D::Algorithms::MovingConvexHull::MovingPoint
moving_points = [
  MovingPoint.new([0, 0], [0, 0]),
  MovingPoint.new([2, 0], [0, 1]),
  MovingPoint.new([2, 2], [0, 0]),
  MovingPoint.new([0, 2], [0, 0])
]
analysis = Geom2D::Algorithms.moving_convex_hull(moving_points, 0, 2)
analysis[:events]    # Interior times where the ordered hull changes
analysis[:intervals] # Open intervals between consecutive event times
~~~

### Moving convex hull

`moving_convex_hull` enumerates the finite set of candidate changes for small inputs:
every pair-coincidence time and every triple-collinearity time strictly inside the analysis
interval. Candidates at the same time are processed as one simultaneous event. After grouping,
the static monotone-chain hull is evaluated once in each adjacent open interval, so sampling any
time in a returned interval gives the same ordered hull as computing a static hull at that time.
Short intervals are never removed based on their duration.

Results distinguish event times and adjacent open intervals:

* `events` contains interior times at which the cyclic hull order of neighboring open intervals
  changes, plus times where distinct points coincide on a hull vertex only at the event itself.
  Other candidate collinearity times are omitted.
* Every entry in `intervals` has `open: true`, endpoint times, ordered counter-clockwise
  representative indices/IDs, and the corresponding moving points.
* At an event itself, a point lying in the interior of a hull edge is not part of the strict
  convex hull. Query the static hull at that exact time if the event-time boundary hull is needed.
* The first point of an always-coincident trajectory represents that group in the ordered hull;
  `coincident_groups` and each interval's `hull_vertex_groups` expose all IDs sharing that
  position throughout the interval.
* Custom point IDs must be unique; event pair/triple indices identify representative points.

`epsilon` (default `Geom2D::Utils.precision`, `1e-10`) is an absolute geometric tolerance for
cross products and coincident coordinates/trajectories. It does not impose a minimum event
spacing. Root finding uses high-precision exact sign tests internally; `time_tolerance` defaults
to zero and should be set only when known floating-point input errors make one simultaneous event
appear as several numerically nearby events.
