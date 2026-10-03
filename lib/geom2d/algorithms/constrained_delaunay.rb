# -*- frozen_string_literal: true -*-
#
#--
# geom2d - 2D Geometric Objects and Algorithms
# Copyright (C) 2018-2023 Thomas Leitner <t_leitner@gmx.at>
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
#++

require 'set'
require 'geom2d/algorithms'
require 'geom2d/utils'
require 'geom2d/polygon_set'

module Geom2D
  module Algorithms

    # Computes a constrained Delaunay triangulation of a Geom2D::PolygonSet.
    #
    # The entry method is ConstrainedDelaunay.run which takes a PolygonSet and returns a
    # Triangulation object containing the triangles as well as their adjacency relations.
    #
    # The polygon set may contain multiple outer contours as well as holes (and islands within
    # holes). Which regions are filled is determined with the even-odd rule over all contours, so
    # the orientation of the contours does not matter. Only the filled regions are triangulated,
    # holes stay empty.
    #
    # Degenerate and special cases are handled with the following stable rules:
    #
    # * Points that are closer together than Utils.precision times the bounding box diagonal are
    #   treated as identical and merged.
    # * Points lying (within tolerance) on a constraint edge split that edge into sub-edges which
    #   are all marked as constrained, so constraint edges are never lost.
    # * For cocircular point quadruples the existing diagonal is kept (no edge flip), making the
    #   result deterministic.
    # * Constraint edges that intersect each other (except at shared endpoints) are rejected by
    #   raising a ConstraintConflictError that describes the concrete conflict.
    #
    # Example:
    #
    #   result = Geom2D::Algorithms::ConstrainedDelaunay.run(polygon_set)
    #   result.triangles  # => [[vertex_index, ...], ...]
    #   result.neighbors  # => [[triangle_index_or_nil, ...], ...]
    #   result.valid?(polygon_set)  # => true
    class ConstrainedDelaunay

      include Utils

      # Error raised when the input constraint edges intersect each other.
      class ConstraintConflictError < StandardError

        # The location of the first conflicting edge as [contour_index, edge_index] pair.
        attr_reader :location1

        # The first conflicting edge as Geom2D::Segment.
        attr_reader :segment1

        # The location of the second conflicting edge as [contour_index, edge_index] pair.
        attr_reader :location2

        # The second conflicting edge as Geom2D::Segment.
        attr_reader :segment2

        # The intersection of the two edges, either a Geom2D::Point or a Geom2D::Segment (in case
        # of overlapping collinear edges).
        attr_reader :intersection

        def initialize(location1, segment1, location2, segment2, intersection)
          @location1 = location1
          @segment1 = segment1
          @location2 = location2
          @segment2 = segment2
          @intersection = intersection
          super("Intersecting constraints: edge #{location1[1]} of contour #{location1[0]} " +
                "#{segment1} conflicts with edge #{location2[1]} of contour #{location2[0]} " +
                "#{segment2} (intersection: #{intersection})")
        end

      end

      # The result of a constrained Delaunay triangulation: the triangles of the filled regions
      # of a polygon set together with their adjacency relations.
      class Triangulation

        include Utils

        # The vertices of the triangulation as Geom2D::Point objects.
        attr_reader :vertices

        # The triangles of the triangulation. Each triangle is an array of three vertex indices
        # into #vertices, ordered counterclockwise.
        attr_reader :triangles

        # The adjacency relations of the triangles: +neighbors[t][k]+ is the index of the
        # triangle sharing the edge opposite vertex +k+ of triangle +t+, or +nil+ if there is no
        # such triangle (i.e. the edge lies on the boundary of the triangulation).
        attr_reader :neighbors

        # The constraint edges of the triangulation as [vertex_index, vertex_index] pairs
        # (sorted ascending).
        attr_reader :constraint_edges

        def initialize(vertices, triangles, neighbors, constraint_edges)
          @vertices = vertices
          @triangles = triangles
          @neighbors = neighbors
          @constraint_edges = constraint_edges
          @constraint_edge_set = Set.new(constraint_edges)
        end

        # Calls the given block once for each triangle, passing an array of the three corner
        # points in counterclockwise order.
        #
        # If no block is given, an Enumerator is returned.
        def each_triangle
          return to_enum(__method__) unless block_given?
          @triangles.each {|a, b, c| yield([@vertices[a], @vertices[b], @vertices[c]]) }
        end

        # Returns the total area of all triangles.
        def area
          each_triangle.sum do |(a, b, c)|
            0.5 * ((b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y))
          end
        end

        # Returns +true+ if the edge between the vertices with the given indices is a
        # constraint edge.
        def constrained_edge?(index1, index2)
          @constraint_edge_set.include?([index1, index2].minmax)
        end

        # Returns +true+ if this triangulation is a valid constrained Delaunay triangulation
        # of the given polygon set. See #validate for the performed checks.
        def valid?(polygon_set)
          validate(polygon_set).empty?
        end

        # Validates this triangulation against the given polygon set and returns an array of
        # problem descriptions (empty if the triangulation is valid).
        #
        # The following checks are performed:
        #
        # * The sum of the triangle areas equals the filled area of the polygon set.
        # * The interiors of the triangles do not overlap each other.
        # * No triangle lies inside a hole or outside the polygon set.
        # * All constraint edges of the polygon set are preserved (possibly split into
        #   collinear sub-edges).
        # * All flippable, non-constraint interior edges satisfy the Delaunay empty-circle
        #   condition.
        def validate(polygon_set)
          compute_tolerances
          problems = []
          problems.concat(area_problems(polygon_set))
          problems.concat(overlap_problems)
          problems.concat(hole_problems(polygon_set))
          problems.concat(constraint_problems(polygon_set))
          problems.concat(delaunay_problems)
          problems
        end

        # Returns the crossing number of the point (x, y) with respect to the given contours
        # (arrays of Geom2D::Point objects). Regions with an odd crossing number are filled
        # when using the even-odd rule.
        def self.crossing_number(contours, x, y)
          number = 0
          contours.each do |contour|
            contour.each_with_index do |p1, i|
              p2 = contour[(i + 1) % contour.length]
              next if (p1.y > y) == (p2.y > y)
              x_intersect = p1.x + (y - p1.y) * (p2.x - p1.x) / (p2.y - p1.y)
              number += 1 if x_intersect > x
            end
          end
          number
        end

        private

        def compute_tolerances
          if @vertices.empty?
            scale = 1.0
          else
            minx, maxx = @vertices.map(&:x).minmax
            miny, maxy = @vertices.map(&:y).minmax
            scale = Math.hypot(maxx - minx, maxy - miny)
            scale = 1.0 if scale < 1.0
          end
          @point_tol = Utils.precision * scale
          @orient_tol = Utils.precision * scale * scale
          @incircle_tol = Utils.precision * scale**4
        end

        def orient(a, b, c)
          det = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)
          if det > @orient_tol then 1 elsif det < -@orient_tol then -1 else 0 end
        end

        def incircle(a, b, c, d)
          ax = a.x - d.x
          ay = a.y - d.y
          bx = b.x - d.x
          by = b.y - d.y
          cx = c.x - d.x
          cy = c.y - d.y
          (ax * ax + ay * ay) * (bx * cy - cx * by) -
            (bx * bx + by * by) * (ax * cy - cx * ay) +
            (cx * cx + cy * cy) * (ax * by - bx * ay)
        end

        def area_problems(polygon_set)
          expected = filled_area(polygon_set)
          if (area - expected).abs > 1e-8 * [expected, 1.0].max
            ["Total triangle area #{area} differs from filled polygon set area #{expected}"]
          else
            []
          end
        end

        # Computes the filled area of the polygon set using the even-odd rule: the absolute
        # area of each contour is added (even nesting depth) or subtracted (odd depth).
        def filled_area(polygon_set)
          contours = polygon_set.polygons.map(&:to_a)
          total = 0.0
          contours.each_with_index do |contour, ci|
            next if contour.length < 3
            contour_area = 0.0
            contour.each_with_index do |p1, i|
              p2 = contour[(i + 1) % contour.length]
              contour_area += p1.x * p2.y - p2.x * p1.y
            end
            vertex = contour[0]
            depth = contours.each_with_index.count do |other, oi|
              oi != ci && Triangulation.crossing_number([other], vertex.x, vertex.y).odd?
            end
            total += (depth.even? ? 1 : -1) * contour_area.abs / 2.0
          end
          total
        end

        def overlap_problems
          problems = []
          points = @triangles.map {|tri| tri.map {|vi| @vertices[vi] } }
          points.each_with_index do |t1, i|
            ((i + 1)...points.length).each do |j|
              t2 = points[j]
              next unless bbox_overlap?(t1, t2)
              if triangles_overlap?(t1, t2)
                problems << "Triangles #{i} and #{j} overlap each other"
              end
            end
          end
          problems
        end

        def bbox_overlap?(t1, t2)
          minx1, maxx1 = t1.map(&:x).minmax
          minx2, maxx2 = t2.map(&:x).minmax
          miny1, maxy1 = t1.map(&:y).minmax
          miny2, maxy2 = t2.map(&:y).minmax
          minx1 < maxx2 + @point_tol && minx2 < maxx1 + @point_tol &&
            miny1 < maxy2 + @point_tol && miny2 < maxy1 + @point_tol
        end

        def triangles_overlap?(t1, t2)
          3.times do |k|
            p1 = t1[k]
            p2 = t1[(k + 1) % 3]
            3.times do |l|
              q1 = t2[l]
              q2 = t2[(l + 1) % 3]
              if orient(p1, p2, q1) * orient(p1, p2, q2) < 0 &&
                  orient(q1, q2, p1) * orient(q1, q2, p2) < 0
                return true
              end
            end
          end
          t2.any? {|q| strictly_inside?(q, t1) } || t1.any? {|p| strictly_inside?(p, t2) }
        end

        def strictly_inside?(p, triangle)
          orient(triangle[0], triangle[1], p) > 0 && orient(triangle[1], triangle[2], p) > 0 &&
            orient(triangle[2], triangle[0], p) > 0
        end

        def hole_problems(polygon_set)
          problems = []
          contours = polygon_set.polygons.map(&:to_a)
          @triangles.each_with_index do |(a, b, c), ti|
            cx = (@vertices[a].x + @vertices[b].x + @vertices[c].x) / 3.0
            cy = (@vertices[a].y + @vertices[b].y + @vertices[c].y) / 3.0
            if Triangulation.crossing_number(contours, cx, cy).even?
              problems << "Triangle #{ti} lies inside a hole or outside the polygon set"
            end
          end
          problems
        end

        def constraint_problems(polygon_set)
          problems = []
          polygon_set.polygons.each_with_index do |polygon, ci|
            points = polygon.to_a
            points.each_with_index do |p, ei|
              q = points[(ei + 1) % points.length]
              next if p == q
              unless constraint_covered?(p, q)
                problems << "Constraint edge #{ei} of contour #{ci} (#{p}-#{q}) is not preserved"
              end
            end
          end
          problems
        end

        # Checks whether the constraint edge from p to q is covered by (collinear) constraint
        # edges of the triangulation.
        def constraint_covered?(p, q)
          dx = q.x - p.x
          dy = q.y - p.y
          length_squared = dx * dx + dy * dy
          intervals = []
          @constraint_edges.each do |(i, j)|
            a = @vertices[i]
            b = @vertices[j]
            next unless collinear_on?(p, q, a) && collinear_on?(p, q, b)
            ta = ((a.x - p.x) * dx + (a.y - p.y) * dy) / length_squared
            tb = ((b.x - p.x) * dx + (b.y - p.y) * dy) / length_squared
            lo, hi = [ta, tb].minmax
            lo = [lo, 0.0].max
            hi = [hi, 1.0].min
            intervals << [lo, hi] if hi > lo
          end
          param_tol = @point_tol / Math.sqrt(length_squared)
          covered = 0.0
          intervals.sort.each do |lo, hi|
            return false if lo > covered + param_tol
            covered = hi if hi > covered
          end
          covered >= 1.0 - param_tol
        end

        def collinear_on?(p, q, a)
          orient(p, q, a) == 0 &&
            (a.x - p.x) * (a.x - q.x) + (a.y - p.y) * (a.y - q.y) <= @orient_tol
        end

        def delaunay_problems
          problems = []
          edge_map = Hash.new {|h, k| h[k] = [] }
          @triangles.each_with_index do |tri, ti|
            3.times do |k|
              i = tri[(k + 1) % 3]
              j = tri[(k + 2) % 3]
              edge_map[[i, j].minmax] << [ti, tri[k]]
            end
          end
          edge_map.each do |(i, j), list|
            next unless list.length == 2
            next if @constraint_edge_set.include?([i, j])
            (t1, c), (_, d) = list
            pc = @vertices[c]
            pd = @vertices[d]
            # the edge is only flippable if the two triangles form a convex quadrilateral
            next unless orient(pc, pd, @vertices[i]) * orient(pc, pd, @vertices[j]) < 0
            a, b, cc = @triangles[t1].map {|vi| @vertices[vi] }
            if incircle(a, b, cc, pd) > @incircle_tol
              problems << "Edge (#{i}, #{j}) violates the Delaunay empty-circle condition"
            end
          end
          problems
        end

      end

      # Computes the constrained Delaunay triangulation of the given Geom2D::PolygonSet and
      # returns a Triangulation object.
      #
      # Raises a ConstraintConflictError if constraint edges of the polygon set intersect each
      # other (except at shared endpoints).
      def self.run(polygon_set)
        new(polygon_set).run
      end

      def initialize(polygon_set) # :nodoc:
        @polygon_set = polygon_set
      end

      def run # :nodoc:
        prepare_input
        return Triangulation.new([], [], [], []) if @nr_input_points < 3

        check_constraint_conflicts
        create_supertriangle
        @nr_input_points.times {|p| insert_point(p) }
        @constraints.each {|u, v,| insert_constraint(u, v) }
        build_result
      end

      private

      # Extracts the contours of the polygon set, merges (nearly) identical points and builds
      # the list of constraint edges [u, v, contour_index, edge_index].
      def prepare_input
        raw = @polygon_set.polygons.map do |polygon|
          polygon.to_a.map {|pt| [pt.x.to_f, pt.y.to_f] }
        end
        all = raw.flatten(1)
        @points_x = []
        @points_y = []
        @contours = []
        @constraints = []
        @nr_input_points = 0
        return if all.empty?

        minx, maxx = all.map(&:first).minmax
        miny, maxy = all.map(&:last).minmax
        @scale = Math.hypot(maxx - minx, maxy - miny)
        @scale = 1.0 if @scale < 1.0
        @point_tol = Utils.precision * @scale
        @orient_tol = Utils.precision * @scale * @scale
        @incircle_tol = Utils.precision * @scale**4

        grid = Hash.new {|h, k| h[k] = [] }
        index_maps = raw.map do |contour|
          contour.map {|x, y| deduped_point_index(grid, x, y) }
        end

        index_maps.each_with_index do |map, ci|
          contour = []
          map.each {|idx| contour << idx if contour.empty? || contour.last != idx }
          contour.pop if contour.length > 1 && contour.first == contour.last
          if contour.length < 3
            raise ArgumentError, "Contour #{ci} has less than three distinct vertices"
          end
          @contours << contour
          contour.each_with_index do |u, ei|
            @constraints << [u, contour[(ei + 1) % contour.length], ci, ei]
          end
        end
        @nr_input_points = @points_x.length
      end

      # Returns the index of the given point, merging it with an existing (nearly) identical
      # point if there is one.
      def deduped_point_index(grid, x, y)
        gx = (x / @point_tol).floor
        gy = (y / @point_tol).floor
        (gx - 1).upto(gx + 1) do |ix|
          (gy - 1).upto(gy + 1) do |iy|
            grid[[ix, iy]].each do |idx|
              if (@points_x[idx] - x).abs < @point_tol && (@points_y[idx] - y).abs < @point_tol
                return idx
              end
            end
          end
        end
        index = @points_x.length
        @points_x << x
        @points_y << y
        grid[[gx, gy]] << index
        index
      end

      #--
      # Constraint conflict detection
      #++

      def check_constraint_conflicts
        @constraints.each_with_index do |c1, i1|
          ((i1 + 1)...@constraints.length).each do |i2|
            check_constraint_pair(c1, @constraints[i2])
          end
        end
      end

      def check_constraint_pair(c1, c2)
        u1, v1, = c1
        u2, v2, = c2
        o_u2 = orient_sign(u1, v1, u2)
        o_v2 = orient_sign(u1, v1, v2)

        if u1 == u2 || u1 == v2 || v1 == u2 || v1 == v2
          # the edges share an endpoint, only a collinear overlap is a conflict
          if o_u2 == 0 && o_v2 == 0
            overlap = collinear_overlap(u1, v1, u2, v2)
            constraint_conflict(c1, c2, overlap) if overlap
          end
          return
        end

        if o_u2 == 0 && o_v2 == 0
          overlap = collinear_overlap(u1, v1, u2, v2)
          constraint_conflict(c1, c2, overlap) if overlap
          return
        end

        o_u1 = orient_sign(u2, v2, u1)
        o_v1 = orient_sign(u2, v2, v1)
        if o_u2 * o_v2 < 0 && o_u1 * o_v1 < 0
          constraint_conflict(c1, c2, segment(u1, v1).intersect(segment(u2, v2)))
        elsif o_u2 == 0 && on_segment?(u1, v1, u2)
          constraint_conflict(c1, c2, point(u2))
        elsif o_v2 == 0 && on_segment?(u1, v1, v2)
          constraint_conflict(c1, c2, point(v2))
        elsif o_u1 == 0 && on_segment?(u2, v2, u1)
          constraint_conflict(c1, c2, point(u1))
        elsif o_v1 == 0 && on_segment?(u2, v2, v1)
          constraint_conflict(c1, c2, point(v1))
        end
      end

      # Returns the overlap of the two collinear segments as Geom2D::Segment, or +nil+ if they
      # do not overlap (touching at a single endpoint is not an overlap).
      def collinear_overlap(u1, v1, u2, v2)
        dx = @points_x[v1] - @points_x[u1]
        dy = @points_y[v1] - @points_y[u1]
        if dx.abs >= dy.abs
          param = lambda {|w| (@points_x[w] - @points_x[u1]) / dx }
          length = dx.abs
        else
          param = lambda {|w| (@points_y[w] - @points_y[u1]) / dy }
          length = dy.abs
        end
        lo, hi = [param.call(u2), param.call(v2)].minmax
        from = [0.0, lo].max
        to = [1.0, hi].min
        return nil if to - from <= @point_tol / length

        Geom2D::Segment.new([@points_x[u1] + from * dx, @points_y[u1] + from * dy],
                            [@points_x[u1] + to * dx, @points_y[u1] + to * dy])
      end

      def constraint_conflict(c1, c2, intersection)
        u1, v1, ci1, ei1 = c1
        u2, v2, ci2, ei2 = c2
        raise ConstraintConflictError.new([ci1, ei1], segment(u1, v1),
                                          [ci2, ei2], segment(u2, v2), intersection)
      end

      def point(i)
        Geom2D::Point(@points_x[i], @points_y[i])
      end

      def segment(i, j)
        Geom2D::Segment.new(point(i), point(j))
      end

      #--
      # Triangulation data structure and geometric predicates
      #++

      # Sign of the orientation of the points a, b, c (given as indices): +1 for
      # counterclockwise, -1 for clockwise and 0 for (nearly) collinear.
      def orient_sign(a, b, c)
        det = (@points_x[b] - @points_x[a]) * (@points_y[c] - @points_y[a]) -
              (@points_x[c] - @points_x[a]) * (@points_y[b] - @points_y[a])
        if det > @orient_tol then 1 elsif det < -@orient_tol then -1 else 0 end
      end

      # Sign of the incircle test of point d with respect to the counterclockwise triangle
      # (a, b, c): +1 if d lies inside the circumcircle, -1 if outside and 0 if (nearly)
      # cocircular.
      def incircle_sign(a, b, c, d)
        ax = @points_x[a] - @points_x[d]
        ay = @points_y[a] - @points_y[d]
        bx = @points_x[b] - @points_x[d]
        by = @points_y[b] - @points_y[d]
        cx = @points_x[c] - @points_x[d]
        cy = @points_y[c] - @points_y[d]
        det = (ax * ax + ay * ay) * (bx * cy - cx * by) -
              (bx * bx + by * by) * (ax * cy - cx * ay) +
              (cx * cx + cy * cy) * (ax * by - bx * ay)
        if det > @incircle_tol then 1 elsif det < -@incircle_tol then -1 else 0 end
      end

      # Returns +true+ if point p lies (within tolerance) on the segment from a to b. The
      # points have to be collinear already.
      def on_segment?(a, b, p)
        (@points_x[p] - @points_x[a]) * (@points_x[p] - @points_x[b]) +
          (@points_y[p] - @points_y[a]) * (@points_y[p] - @points_y[b]) <= @orient_tol
      end

      def create_supertriangle
        minx, maxx = @points_x[0...@nr_input_points].minmax
        miny, maxy = @points_y[0...@nr_input_points].minmax
        cx = (minx + maxx) / 2.0
        cy = (miny + maxy) / 2.0
        m = 10.0 * @scale
        @tris = []
        @edges = {}
        @constrained = {}
        @vertex_tri = {}
        s0 = add_vertex(cx - 2.0 * m, cy - m)
        s1 = add_vertex(cx + 2.0 * m, cy - m)
        s2 = add_vertex(cx, cy + 2.0 * m)
        @last_tri = add_triangle(s0, s1, s2)
      end

      def add_vertex(x, y)
        @points_x << x
        @points_y << y
        @points_x.length - 1
      end

      def add_triangle(a, b, c)
        b, c = c, b if orient_sign(a, b, c) < 0
        index = @tris.length
        @tris << [a, b, c]
        @vertex_tri[a] = index
        @vertex_tri[b] = index
        @vertex_tri[c] = index
        each_triangle_edge([a, b, c]) do |i, j|
          key = i < j ? [i, j] : [j, i]
          (@edges[key] ||= []) << index
        end
        @last_tri = index
        index
      end

      def remove_triangle(index)
        tri = @tris[index]
        @tris[index] = nil
        each_triangle_edge(tri) do |i, j|
          key = i < j ? [i, j] : [j, i]
          list = @edges[key]
          list.delete(index)
          @edges.delete(key) if list.empty?
        end
      end

      def each_triangle_edge(tri)
        yield(tri[0], tri[1])
        yield(tri[1], tri[2])
        yield(tri[2], tri[0])
      end

      # Returns the other triangle adjacent to the edge (i, j), or +nil+.
      def neighbor_triangle(index, i, j)
        key = i < j ? [i, j] : [j, i]
        @edges[key]&.find {|ti| ti != index }
      end

      #--
      # Point insertion (incremental, with Lawson edge flips)
      #++

      # Inserts the point with the given index into the triangulation. If the point coincides
      # with an existing vertex, nothing is inserted and the vertex index is returned.
      def insert_point(p)
        t = locate_triangle(p)
        a, b, c = @tris[t]
        [a, b, c].each do |v|
          return v if (@points_x[v] - @points_x[p]).abs < @point_tol &&
            (@points_y[v] - @points_y[p]).abs < @point_tol
        end
        edge = nil
        [[a, b], [b, c], [c, a]].each do |i, j|
          if orient_sign(i, j, p) == 0 && on_segment?(i, j, p)
            edge = [i, j]
            break
          end
        end
        edge ? split_edge(t, edge[0], edge[1], p) : split_triangle(t, p)
        p
      end

      # Returns the index of the triangle containing the point p, using a visibility walk
      # starting from the last created triangle.
      def locate_triangle(p)
        t = @last_tri
        t = @tris.index {|tri| tri } unless t && @tris[t]
        iterations = 0
        loop do
          crossed = nil
          each_triangle_edge(@tris[t]) do |i, j|
            if orient_sign(i, j, p) < 0
              crossed = [i, j]
              break
            end
          end
          break unless crossed
          nt = neighbor_triangle(t, *crossed)
          break unless nt
          t = nt
          iterations += 1
          break if iterations > @tris.length
        end
        unless @tris[t] && point_in_triangle?(@tris[t], p)
          t = @tris.index {|tri| tri && point_in_triangle?(tri, p) }
        end
        raise "Unable to locate a triangle containing point #{p}" unless t
        @last_tri = t
        t
      end

      def point_in_triangle?(tri, p)
        orient_sign(tri[0], tri[1], p) >= 0 && orient_sign(tri[1], tri[2], p) >= 0 &&
          orient_sign(tri[2], tri[0], p) >= 0
      end

      def split_triangle(t, p)
        a, b, c = @tris[t]
        remove_triangle(t)
        t1 = add_triangle(a, b, p)
        t2 = add_triangle(b, c, p)
        t3 = add_triangle(c, a, p)
        legalize([[a, b, t1], [b, c, t2], [c, a, t3]], p)
      end

      def split_edge(t, i, j, p)
        opp = (@tris[t] - [i, j]).first
        t2 = neighbor_triangle(t, i, j)
        d = t2 ? (@tris[t2] - [i, j]).first : nil
        remove_triangle(t)
        remove_triangle(t2) if t2
        stack = []
        ta = add_triangle(i, p, opp)
        tb = add_triangle(p, j, opp)
        stack.push([i, opp, ta], [j, opp, tb])
        if t2
          tc = add_triangle(j, p, d)
          td = add_triangle(p, i, d)
          stack.push([j, d, tc], [i, d, td])
        end
        legalize(stack, p)
      end

      # Legalizes the edges [i, j, triangle] on the stack by flipping them if the opposite
      # vertex lies strictly inside the circumcircle. Cocircular configurations keep the
      # existing diagonal (stable rule).
      def legalize(stack, p)
        iterations = 0
        until stack.empty?
          i, j, t = stack.pop
          tri = @tris[t]
          next unless tri
          t2 = neighbor_triangle(t, i, j)
          next unless t2
          d = (@tris[t2] - [i, j]).first
          next unless incircle_sign(tri[0], tri[1], tri[2], d) > 0

          x = y = nil
          3.times do |k|
            if tri[k] == i && tri[(k + 1) % 3] == j
              x, y = i, j
              break
            elsif tri[k] == j && tri[(k + 1) % 3] == i
              x, y = j, i
              break
            end
          end
          remove_triangle(t)
          remove_triangle(t2)
          nt1 = add_triangle(x, p, d)
          nt2 = add_triangle(p, y, d)
          stack.push([x, d, nt1], [y, d, nt2])
          iterations += 1
          raise "Edge legalization did not terminate" if iterations > 1_000_000
        end
      end

      #--
      # Constraint edge insertion
      #++

      # Inserts the constraint edge from u to v. Existing vertices lying on the edge split it
      # into sub-edges (stable collinear rule) which are inserted one by one.
      def insert_constraint(u, v)
        on = []
        @nr_input_points.times do |w|
          next if w == u || w == v
          on << w if orient_sign(u, v, w) == 0 && on_segment?(u, v, w)
        end
        dx = @points_x[v] - @points_x[u]
        dy = @points_y[v] - @points_y[u]
        on.sort_by! {|w| (@points_x[w] - @points_x[u]) * dx + (@points_y[w] - @points_y[u]) * dy }
        ([u] + on + [v]).each_cons(2) {|a, b| insert_constraint_edge(a, b) }
      end

      # Inserts the constraint edge from u to v (no other vertex lies on the segment) by
      # removing the crossed triangles and retriangulating the two resulting cavity sides.
      def insert_constraint_edge(u, v)
        key = u < v ? [u, v] : [v, u]
        if @edges.key?(key)
          @constrained[key] = true
          return
        end
        crossed = crossed_triangles(u, v)
        upper, lower = cavity_chains(crossed, u, v)
        crossed.each {|t| remove_triangle(t) }
        triangulate_cavity(upper)
        triangulate_cavity(lower)
        @constrained[key] = true
      end

      # Returns the ordered list of triangles crossed by the segment from u to v.
      def crossed_triangles(u, v)
        crossed = []
        t = find_first_crossed_triangle(u, v)
        entry = nil
        loop do
          crossed << t
          tri = @tris[t]
          break if tri.include?(v)
          exit_edge = nil
          each_triangle_edge(tri) do |i, j|
            key = i < j ? [i, j] : [j, i]
            next if key == entry
            oi = orient_sign(u, v, i)
            oj = orient_sign(u, v, j)
            if oi != 0 && oj != 0 && oi != oj
              exit_edge = key
              break
            end
          end
          raise "Unable to insert constraint edge" unless exit_edge
          t = neighbor_triangle(t, *exit_edge)
          raise "Unable to insert constraint edge" unless t
          entry = exit_edge
          raise "Unable to insert constraint edge" if crossed.length > @tris.length
        end
        crossed
      end

      # Returns the triangle incident to vertex u that is crossed first by the segment from u
      # to v, walking around the fan of triangles at u.
      def find_first_crossed_triangle(u, v)
        t = @vertex_tri[u]
        unless t && @tris[t] && @tris[t].include?(u)
          t = @tris.index {|tri| tri && tri.include?(u) }
        end
        raise "Unable to insert constraint edge" unless t
        iterations = 0
        loop do
          tri = @tris[t]
          k = tri.index(u)
          a = tri[(k + 1) % 3]
          b = tri[(k + 2) % 3]
          oa = orient_sign(u, a, v)
          ob = orient_sign(u, b, v)
          return t if oa > 0 && ob < 0
          t = neighbor_triangle(t, u, oa < 0 ? a : b)
          raise "Unable to insert constraint edge" unless t
          iterations += 1
          raise "Unable to insert constraint edge" if iterations > @tris.length
        end
      end

      # Computes the two boundary chains (upper/left and lower/right, each from u to v) of the
      # cavity left by the crossed triangles.
      def cavity_chains(crossed, u, v)
        count = Hash.new(0)
        crossed.each do |t|
          each_triangle_edge(@tris[t]) do |i, j|
            count[i < j ? [i, j] : [j, i]] += 1
          end
        end
        adjacency = Hash.new {|h, k| h[k] = [] }
        count.each do |(i, j), c|
          next unless c == 1
          adjacency[i] << j
          adjacency[j] << i
        end
        n1, n2 = adjacency[u]
        chain1 = walk_cavity_chain(adjacency, u, n1, v)
        chain2 = walk_cavity_chain(adjacency, u, n2, v)
        orient_sign(u, v, n1) > 0 ? [chain1, chain2] : [chain2, chain1]
      end

      def walk_cavity_chain(adjacency, u, first, v)
        chain = [u, first]
        prev = u
        cur = first
        until cur == v
          nxt = adjacency[cur].find {|w| w != prev }
          chain << nxt
          prev = cur
          cur = nxt
        end
        chain
      end

      # Triangulates one side of a constraint cavity (a polygonal chain from the constraint
      # start to its end, with all intermediate vertices on the same side) with Delaunay
      # triangles.
      def triangulate_cavity(chain)
        return if chain.length < 3
        if chain.length == 3
          add_triangle(chain[0], chain[1], chain[2])
          return
        end
        first = chain.first
        last = chain.last
        best = nil
        (1...(chain.length - 1)).each do |ci|
          cand = chain[ci]
          ccw = orient_sign(first, cand, last) > 0
          good = true
          (1...(chain.length - 1)).each do |wi|
            next if wi == ci
            w = chain[wi]
            inside = ccw ? incircle_sign(first, cand, last, w) > 0 :
              incircle_sign(first, last, cand, w) > 0
            if inside
              good = false
              break
            end
          end
          if good
            best = ci
            break
          end
        end
        best ||= 1
        add_triangle(first, chain[best], last)
        triangulate_cavity(chain[0..best])
        triangulate_cavity(chain[best..-1])
      end

      #--
      # Region classification and result construction
      #++

      # Returns the crossing number of the point (x, y) with respect to all contours.
      def crossing_number(x, y)
        number = 0
        @contours.each do |contour|
          contour.each_with_index do |u, ui|
            v = contour[(ui + 1) % contour.length]
            y1 = @points_y[u]
            y2 = @points_y[v]
            next if (y1 > y) == (y2 > y)
            x_intersect = @points_x[u] + (y - y1) * (@points_x[v] - @points_x[u]) / (y2 - y1)
            number += 1 if x_intersect > x
          end
        end
        number
      end

      # Builds the final Triangulation object: triangles touching the supertriangle are
      # removed and the remaining triangles are kept iff their centroid lies in a filled
      # region of the polygon set (odd crossing number, even-odd rule).
      def build_result
        kept = []
        @tris.each_with_index do |tri, ti|
          next unless tri
          next if tri.any? {|w| w >= @nr_input_points }
          cx = (@points_x[tri[0]] + @points_x[tri[1]] + @points_x[tri[2]]) / 3.0
          cy = (@points_y[tri[0]] + @points_y[tri[1]] + @points_y[tri[2]]) / 3.0
          kept << [ti, tri] if crossing_number(cx, cy).odd?
        end

        vertex_map = {}
        vertices = []
        kept.each do |(_, tri)|
          tri.each do |w|
            unless vertex_map.key?(w)
              vertex_map[w] = vertices.length
              vertices << point(w)
            end
          end
        end

        tri_map = {}
        kept.each_with_index {|(ti, _), ni| tri_map[ti] = ni }
        triangles = kept.map {|(_, tri)| tri.map {|w| vertex_map[w] } }
        neighbors = kept.each_with_index.map do |(ti, tri), _|
          tri_neighbors = [nil, nil, nil]
          3.times do |k|
            i = tri[(k + 1) % 3]
            j = tri[(k + 2) % 3]
            key = i < j ? [i, j] : [j, i]
            @edges[key]&.each do |t2|
              tri_neighbors[k] = tri_map[t2] if t2 != ti && tri_map.key?(t2)
            end
          end
          tri_neighbors
        end

        constraint_edges = []
        @constrained.each_key do |(i, j)|
          constraint_edges << [vertex_map[i], vertex_map[j]].minmax if vertex_map.key?(i) &&
            vertex_map.key?(j)
        end

        Triangulation.new(vertices, triangles, neighbors, constraint_edges)
      end

    end

  end
end
