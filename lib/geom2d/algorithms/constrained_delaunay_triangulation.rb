# -*- frozen_string_literal: true -*-
#
#--
# geom2d - 2D Geometric Objects and Algorithms
# Copyright (C) 2018-2023 Thomas Leitner <t_leitner@gmx.at>
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
#++

require 'geom2d/algorithms'
require 'geom2d/polygon_set'

module Geom2D
  module Algorithms

    # Performs a constrained Delaunay triangulation (CDT) of a planar straight-line graph that may
    # consist of multiple outer contours and holes.
    #
    # The input is a Geom2D::PolygonSet (or a single Geom2D::Polygon). Outer contours have to be
    # oriented counterclockwise and holes clockwise, following the convention used by the boolean
    # polygon operations (see Geom2D::Algorithms::PolygonOperation). Contours may not cross each
    # other; if they do, a ConstrainedDelaunayTriangulation::ConstraintError is raised that lists
    # all conflicting edges.
    #
    # The algorithm first builds an unconstrained Delaunay triangulation of all vertices using the
    # Bowyer-Watson algorithm, then restores the constraint edges by removing the triangles
    # crossed by each constraint edge and re-triangulating the resulting cavities with ear
    # clipping. Afterwards all non-constraint edges are made Delaunay using edge flips. Finally the
    # triangles outside the planar region (i.e. outside all outer contours and inside holes) are
    # classified and removed using a winding-number test.
    #
    # Stability rules (the input vertex order is never changed):
    #
    # * Collinear vertices on a constraint edge are handled by splitting the edge at such
    #   vertices; collinear points never cause an edge flip.
    # * Co-circular quadrilaterals are not flipped: the edge created first (i.e. the one induced
    #   by the input vertex order) is kept. This makes the result deterministic.
    # * All geometric predicates (orientation and in-circle) are evaluated using exact rational
    #   arithmetic, so there are no floating-point precision related failures.
    class ConstrainedDelaunayTriangulation

      include Utils

      # Raised if the input constraints are invalid, e.g. if two constraint edges properly cross
      # each other.
      #
      # The attribute #conflicts contains an array of hashes describing the detected conflicts,
      # each with the keys :edges (the two Geom2D::Segment objects that conflict) and :point (the
      # Geom2D::Point of the proper intersection, if any).
      class ConstraintError < ArgumentError

        # The array of detected conflicts.
        attr_reader :conflicts

        # Creates a new ConstraintError for the given conflicts.
        def initialize(conflicts)
          @conflicts = conflicts
          super("Invalid input constraints: #{conflicts.size} conflicting edge pair(s) found")
        end

      end

      # A triangle of the triangulation result.
      #
      # The vertex indices refer to Triangulation#points, and the neighbor indices to
      # Triangulation#triangles. The entry +neighbors[i]+ is the triangle adjacent across the edge
      # opposite +vertices[i]+, or +nil+ if the edge is on the boundary of the triangulated
      # region.
      Triangle = Struct.new(:vertices, :neighbors) do

        # The points array of the owning triangulation; set when the Triangulation result is
        # created.
        attr_accessor :points

        # Returns the signed area of the triangle (positive for counterclockwise orientation).
        def area
          pts = vertices.map {|index| points[index] }
          0.5 * ((pts[1].x - pts[0].x) * (pts[2].y - pts[0].y) -
                  (pts[2].x - pts[0].x) * (pts[1].y - pts[0].y))
        end

        def inspect # :nodoc:
          "#<#{self.class} vertices=#{vertices} neighbors=#{neighbors}>"
        end
        alias to_s inspect

      end

      # The result of a constrained Delaunay triangulation.
      class Triangulation

        # The array of used points (Geoem2D::Point objects), indexed by the triangle vertex
        # indices.
        attr_reader :points

        # The array of Triangle objects of the triangulated region.
        attr_reader :triangles

        # The set of canonical, unordered constraint edges as arrays of two vertex indices.
        attr_reader :constraint_edges

        def initialize(points, triangles, constraint_edges) # :nodoc:
          @points = points
          @triangles = triangles
          @constraint_edges = constraint_edges
          @triangles.each {|triangle| triangle.points = @points }
        end

        # Returns +true+ if the edge between the two given vertex indices is a constraint edge.
        def constraint_edge?(node1, node2)
          @constraint_edges.include?([node1, node2].minmax)
        end

        # Returns the sum of the (positive) areas of all triangles.
        def area
          triangles.sum(&:area)
        end

        # Calls the block once for each triangle, yielding it and the triangle's index.
        def each_triangle(&block)
          @triangles.each_with_index(&block)
        end

        def inspect # :nodoc:
          "#<#{self.class}: #{triangles.size} triangles, #{points.size} points>"
        end
        alias to_s inspect

      end

      # Internal mesh triangle. +v+ holds the three vertex node indices in counterclockwise order,
      # +n+ the three neighboring mesh triangles (nil if there is none), where +n[i]+ is opposite
      # +v[i]+. +alive+ marks triangles that have been removed by the algorithm.
      class Tri # :nodoc: all
        attr_accessor :v, :n, :alive

        def initialize(v0, v1, v2)
          @v = [v0, v1, v2]
          @n = [nil, nil, nil]
          @alive = true
        end

        def vertex_index(node)
          @v.index(node)
        end
      end

      # Performs the constrained Delaunay triangulation of the given Geom2D::PolygonSet or
      # Geom2D::Polygon and returns a Triangulation result.
      def self.run(polygon_set)
        new(polygon_set).run
      end

      # Creates a new CDT algorithm object for the given Geom2D::PolygonSet or Geom2D::Polygon.
      def initialize(polygon_set)
        polygon_set = Geom2D::PolygonSet(polygon_set) if polygon_set.kind_of?(Geom2D::Polygon)
        @polygons = polygon_set.polygons

        @nodes = []       # [x, y] float coordinate arrays, in first-seen order
        @rpoints = {}     # canonical [x, y] => node index
        @coord = []       # node index => [Rational x, Rational y]
        @tris = []        # all (possibly dead) mesh triangles
        @constraints = [] # canonical atomic constraint edges as sorted [node1, node2]
      end

      # Performs the triangulation and returns the Triangulation result.
      def run
        @contours = extract_contours
        validate_constraints
        build_initial_triangulation
        recover_constraints
        delaunay_flip
        classify_and_output(@contours)
      end

      private

      # Exact geometric predicates

      # Returns the sign of the 2D orientation determinant of the three nodes, using exact
      # rational arithmetic.
      def orient(a, b, c)
        pa = @coord[a]
        pb = @coord[b]
        pc = @coord[c]
        value = (pb[0] - pa[0]) * (pc[1] - pa[1]) -
                (pc[0] - pa[0]) * (pb[1] - pa[1])
        value.positive? ? 1 : (value.negative? ? -1 : 0)
      end

      # Returns the sign of the in-circle determinant for the counterclockwise triangle a-b-c and
      # the point d: positive if d is inside the circumcircle, zero if co-circular, negative if it
      # is outside.
      def incircle(a, b, c, d)
        pa = @coord[a]; pb = @coord[b]; pc = @coord[c]; pd = @coord[d]
        adx = pa[0] - pd[0]; ady = pa[1] - pd[1]
        bdx = pb[0] - pd[0]; bdy = pb[1] - pd[1]
        cdx = pc[0] - pd[0]; cdy = pc[1] - pd[1]
        det = (adx * adx + ady * ady) * (bdx * cdy - cdx * bdy) -
              (bdx * bdx + bdy * bdy) * (adx * cdy - cdx * ady) +
              (cdx * cdx + cdy * cdy) * (adx * bdy - bdx * ady)
        det.positive? ? 1 : (det.negative? ? -1 : 0)
      end

      # Returns +true+ if point p lies on the closed segment between a and b.
      def on_segment?(p, a, b)
        orient(a, b, p) == 0 &&
          ((@coord[p][0] - @coord[a][0]) * (@coord[p][0] - @coord[b][0]) <= 0) &&
          ((@coord[p][1] - @coord[a][1]) * (@coord[p][1] - @coord[b][1]) <= 0)
      end

      # Returns the vector difference b - a of two node coordinate pairs.
      def diff(b, a)
        [@coord[b][0] - @coord[a][0], @coord[b][1] - @coord[a][1]]
      end

      # Returns the signed doubled area of the node polygon.
      def polygon_area(nodes)
        area = Rational(0)
        nodes.each_with_index do |b, i|
          a = @coord[nodes[(i - 1) % nodes.size]]
          area += a[0] * @coord[b][1] - @coord[b][0] * a[1]
        end
        area
      end

      # Input processing

      # Registers a coordinate pair as a mesh node (deduplicating equal coordinates) and returns
      # its node index.
      def node_for(point)
        point = Geom2D::Point(point)
        key = [Float(point.x), Float(point.y)]
        @rpoints.fetch(key) do
          index = @nodes.size
          @nodes << key
          @rpoints[key] = index
          @coord << [Rational(key[0]), Rational(key[1])]
          index
        end
      end

      # Extracts the contours from the input polygons, registers all nodes and fills the list of
      # atomic constraint edges. Returns an array of contour descriptors: hashes with :nodes (the
      # node indices) and :hole (true/false).
      def extract_contours
        contours = @polygons.map do |polygon|
          nodes = polygon.each_vertex.map {|point| node_for(point) }

          # remove consecutive duplicate nodes, including a duplicated closing node
          cleaned = []
          nodes.each {|node| cleaned << node if cleaned.empty? || cleaned.last != node }
          cleaned.pop if cleaned.size > 1 && cleaned.first == cleaned.last

          if cleaned.size < 3
            raise ArgumentError, "Polygon with fewer than 3 distinct vertices cannot be " \
                                "triangulated: #{polygon.inspect}"
          end

          area = polygon_area(cleaned)
          if area == 0
            raise ArgumentError, "Degenerate polygon with zero area cannot be triangulated: " \
                                "#{polygon.inspect}"
          end

          { nodes: cleaned, hole: area.negative? }
        end

        raw_edges = contours.flat_map do |contour|
          contour[:nodes].each_with_index.map do |node, i|
            [node, contour[:nodes][(i + 1) % contour[:nodes].size]].minmax
          end
        end.uniq.sort

        @raw_edges = raw_edges
        @constraints = atomic_constraint_edges(raw_edges).sort.uniq
        contours
      end

      # Splits the given unique, canonical constraint edges at any registered node lying in the
      # interior of an edge, returning the resulting list of canonical atomic edges.
      def atomic_constraint_edges(raw_edges)
        raw_edges.flat_map do |a, b|
          direction = diff(b, a)
          inner = @nodes.each_index.select do |node|
            node != a && node != b && on_segment?(node, a, b)
          end
          ordered = ([a, b] + inner).sort_by do |node|
            p = diff(node, a)
            p[0] * direction[0] + p[1] * direction[1]
          end
          ordered.each_cons(2).reject {|c, d| c == d }.map {|c, d| [c, d].minmax }
        end
      end

      # Validates the atomic constraint edges: proper intersections (including partial collinear
      # overlaps, which are not expressible as atomic edges) are rejected with all conflicts
      # listed in a ConstraintError.
      def validate_constraints
        conflicts = []
        @raw_edges.each_with_index do |(a, b), i|
          (i + 1).upto(@raw_edges.size - 1) do |j|
            c, d = @raw_edges[j]
            next unless orient(a, b, c) == 0 && orient(a, b, d) == 0
            # Project all endpoints onto the common line direction: a non-empty overlap that
            # consists of more than a single point is a partial/full overlap and invalid.
            direction = diff(b, a)
            length_sq = direction[0] * direction[0] + direction[1] * direction[1]
            t = lambda do |node|
              p = diff(node, a)
              (p[0] * direction[0] + p[1] * direction[1]) / length_sq
            end
            t1 = [t.call(a), t.call(b)].minmax
            t2 = [t.call(c), t.call(d)].minmax
            lo = [t1[0], t2[0]].max
            hi = [t1[1], t2[1]].min
            if hi > lo
              # Coincident intervals are a shared boundary. If one edge fully contains the
              # other, the longer edge is atomized at the contained endpoints (T-junction).
              # Any remaining non-trivial overlap is a dangling partial overlap and invalid.
              contained1 = t1[0] >= t2[0] && t1[1] <= t2[1]
              contained2 = t2[0] >= t1[0] && t2[1] <= t1[1]
              conflicts << conflict(a, b, c, d, nil) unless contained1 || contained2
            end
          end
        end
        @constraints.each_with_index do |(a, b), i|
          (i + 1).upto(@constraints.size - 1) do |j|
            c, d = @constraints[j]
            next if [a, b, c, d].uniq.size < 4 # shared vertices are fine

            o1 = orient(a, b, c)
            o2 = orient(a, b, d)
            o3 = orient(c, d, a)
            o4 = orient(c, d, b)
            if o1 == 0 && on_segment?(c, a, b) && o2 == 0 && on_segment?(d, a, b)
              conflicts << conflict(a, b, c, d, nil)
            elsif (o1 * o2).negative? && (o3 * o4).negative?
              conflicts << conflict(a, b, c, d, proper_intersection(a, b, c, d))
            end
          end
        end
        raise ConstraintError.new(conflicts) unless conflicts.empty?
      end

      # Builds a conflict descriptor for the two atomic edges a-b and c-d.
      def conflict(a, b, c, d, point)
        {
          edges: [Geom2D::Segment.new(@nodes[a], @nodes[b]),
                  Geom2D::Segment.new(@nodes[c], @nodes[d])],
          point: point,
        }
      end

      # Computes the exact intersection point of the two properly intersecting segments a-b and
      # c-d.
      def proper_intersection(a, b, c, d)
        pa = @coord[a]; pc = @coord[c]
        rvec = diff(b, a)
        svec = diff(d, c)
        denominator = rvec[0] * svec[1] - rvec[1] * svec[0]
        t = ((pc[0] - pa[0]) * svec[1] - (pc[1] - pa[1]) * svec[0]) / denominator
        Geom2D::Point.new((pa[0] + t * rvec[0]).to_f,
                          (pa[1] + t * rvec[1]).to_f)
      end

      # Mesh maintenance

      # Rebuilds the adjacency information of all alive triangles: the edge map (canonical edge =>
      # triangles sharing it) and the node-to-triangle map.
      def rebuild_structure
        @edge_map = Hash.new {|hash, key| hash[key] = [] }
        @node_tris = Hash.new {|hash, key| hash[key] = [] }
        @tris.each do |tri|
          next unless tri.alive
          3.times do |i|
            edge = [tri.v[(i + 1) % 3], tri.v[(i + 2) % 3]].minmax
            @edge_map[edge] << [tri, i]
            @node_tris[tri.v[i]] << tri
          end
        end
      end

      # Rebuilds the per-triangle neighbor arrays from the current edge map.
      def relink_triangles
        @tris.each do |tri|
          next unless tri.alive
          3.times do |i|
            edge = [tri.v[(i + 1) % 3], tri.v[(i + 2) % 3]].minmax
            other = @edge_map[edge].map(&:first).find {|candidate| candidate != tri }
            tri.n[i] = other
          end
        end
      end

      # Initial Delaunay triangulation

      # Builds the Delaunay triangulation of all input vertices using the Bowyer-Watson algorithm,
      # starting from a sufficiently large super triangle.
      def build_initial_triangulation
        xs = @nodes.map(&:first)
        ys = @nodes.map(&:last)
        min_x, max_x = xs.minmax
        min_y, max_y = ys.minmax
        span = [[max_x - min_x, max_y - min_y].max, 1].max
        mid_x = (min_x + max_x) / 2.0
        mid_y = (min_y + max_y) / 2.0
        s0 = node_for([mid_x - 20 * span, mid_y - span])
        s1 = node_for([mid_x + 20 * span, mid_y - span])
        s2 = node_for([mid_x, mid_y + 20 * span])
        @super_nodes = [s0, s1, s2]
        @tris << Tri.new(s0, s1, s2)
        rebuild_structure

        # only insert the actual input vertices; the super triangle nodes have the highest indices
        @nodes[0...-3].each_with_index {|_, index| insert_node(index) }
      end

      # Inserts the given node into the current Delaunay triangulation using Bowyer-Watson.
      def insert_node(node)
        cavity = find_cavity(node)
        boundary = []
        cavity.each_key do |tri|
          3.times do |i|
            other = tri.n[i]
            next if other && cavity.include?(other)
            e0 = tri.v[(i + 1) % 3]
            e1 = tri.v[(i + 2) % 3]
            boundary << [e0, e1, other, tri]
          end
        end

        cavity.each_key {|tri| tri.alive = false }
        boundary.sort_by {|e0, e1,| [e0, e1] }.each do |e0, e1, outside, removed|
          # The boundary edge e0->e1 is counterclockwise around the removed triangle, so the
          # removed triangle lies on its left; the new node therefore lies on its right and the
          # new triangle is CCW when ordered as e0-e1-node.
          tri = Tri.new(e0, e1, node)
          tri.n[tri.vertex_index(node)] = outside
          @tris << tri
        end

        rebuild_structure
        relink_triangles
      end

      # Finds the Bowyer-Watson cavity for the node, i.e. all alive triangles whose open
      # circumcircle contains the node. Co-circular triangles are not part of the cavity.
      def find_cavity(node)
        tri = locate_triangle(node)
        cavity = {}
        queue = [tri]
        until queue.empty?
          current = queue.shift
          next if cavity.key?(current)
          if incircle(current.v[0], current.v[1], current.v[2], node) > 0
            cavity[current] = true
            current.n.each {|neighbor| queue << neighbor if neighbor }
          end
        end
        cavity
      end

      # Locates an alive triangle containing the node, walking towards it through the mesh.
      def locate_triangle(node, start = nil)
        tri = start || @tris.find(&:alive)
        visited = {}
        loop do
          raise "Mesh walk did not terminate while locating node #{node}" if visited[tri]
          visited[tri] = true
          next_tri = nil
          3.times do |i|
            if orient(tri.v[i], tri.v[(i + 1) % 3], node) < 0
              next_tri = tri.n[(i + 2) % 3]
              raise "Node #{node} outside of triangulation domain during walk" unless next_tri
              break
            end
          end
          next_tri ? tri = next_tri : (return tri)
        end
      end

      # Constraint recovery

      # Restores all atomic constraint edges by re-triangulating the mesh corridors crossed by
      # each missing edge.
      def recover_constraints
        @constraints.each do |edge|
          next if @edge_map.include?(edge)
          recover_edge(*edge)
        end
      end

      # Restores the single constraint edge a-b.
      #
      # Because the constraint edges were atomized, no other mesh vertex can lie in the interior
      # of a-b (every input vertex is already in the mesh). The ray a->b therefore only exits
      # triangles through the interior of edges, which makes the corridor walk unambiguous.
      def recover_edge(a, b)
        tri = start_triangle_for_ray(a, b)
        corridor = [tri]
        entry_t = Rational(0)
        until triangle_contains_node?(tri, b)
          edge_data = crossed_edge(tri, a, b, entry_t)
          raise "Corridor walk for constraint edge #{a}-#{b} left the mesh" unless edge_data
          t = edge_data[0]
          other = edge_data[1]
          entry_t = t
          raise "Corridor walk for constraint edge #{a}-#{b} reached the outer boundary" unless other
          break if corridor.include?(other)
          tri = other
          corridor << tri
        end
        retriangulate_corridor(corridor, a, b)
      end

      # Returns a triangle incident to node a whose interior (or an incident edge) is crossed by
      # the ray a->b.
      def start_triangle_for_ray(a, b)
        straddling = @node_tris[a].select do |tri|
          k = tri.vertex_index(a)
          p1 = tri.v[(k + 1) % 3]
          p2 = tri.v[(k + 2) % 3]
          o1 = orient(a, b, p1)
          o2 = orient(a, b, p2)
          (o1 * o2).negative?
        end
        forward = straddling.select do |tri|
          k = tri.vertex_index(a)
          p1 = tri.v[(k + 1) % 3]
          p2 = tri.v[(k + 2) % 3]
          t = ray_edge_parameter(a, b, p1, p2)
          t.positive?
        end
        forward.min_by do |tri|
          k = tri.vertex_index(a)
          ray_edge_parameter(a, b, tri.v[(k + 1) % 3], tri.v[(k + 2) % 3])
        end || straddling.first || @node_tris[a].first
      end

      # Returns +true+ if the closed counterclockwise triangle contains node p.
      def triangle_contains_node?(tri, p)
        3.times {|i| return false if orient(tri.v[i], tri.v[(i + 1) % 3], p) < 0 }
        true
      end

      # Finds the edge of the triangle crossed first by the ray a->b (with t in (0, 1]), ignoring
      # edges containing node a, and returns [neighbor_triangle, e0, e1]. Returns +nil+ if the ray
      # does not exit through any edge (i.e. node b is inside the triangle).
      def crossed_edge(tri, a, b, entry_t)
        best = nil
        3.times do |i|
          e0 = tri.v[i]
          e1 = tri.v[(i + 1) % 3]
          if e0 == b || e1 == b
            # b is a vertex of the edge: the ray hits it at t = 1
            t = Rational(1)
            next unless t > entry_t
            neighbor = tri.n[(i + 2) % 3]
            best = [t, neighbor, e0, e1] if best.nil? || t < best[0]
          elsif e0 == a || e1 == a
            next # edges through the start node are handled in the first triangle
          else
            o0 = orient(a, b, e0)
            o1 = orient(a, b, e1)
            next unless (o0 * o1).negative?
            t = ray_edge_parameter(a, b, e0, e1)
            next unless t > entry_t && t < 1
            if best.nil? || t < best[0]
              best = [t, tri.n[(i + 2) % 3], e0, e1]
            end
          end
        end
        best
      end

      # Returns the exact parameter t with which the ray a->b properly intersects edge e0-e1.
      def ray_edge_parameter(a, b, e0, e1)
        pa = @coord[a]; pe = @coord[e0]
        rvec = diff(b, a)
        svec = diff(e1, e0)
        denominator = rvec[0] * svec[1] - rvec[1] * svec[0]
        ((pe[0] - pa[0]) * svec[1] - (pe[1] - pa[1]) * svec[0]) / denominator
      end

      # Removes the corridor triangles and re-triangulates the resulting polygon split by the
      # constraint edge a-b using ear clipping.
      def retriangulate_corridor(corridor, a, b)
        directed = Hash.new(0)
        corridor.each do |tri|
          3.times do |i|
            e0 = tri.v[i]
            e1 = tri.v[(i + 1) % 3]
            other = tri.n[(i + 2) % 3]
            if other && corridor.include?(other)
              # interior edges appear twice in opposite directions
            else
              directed[[e0, e1]] += 1
            end
          end
        end

        boundary = boundary_loop(directed, a)
        raise "Corridor boundary is not a single closed loop" unless boundary

        # build the two boundary chains from b back to a
        chain = boundary[0...-1]
        bpos = chain.index(b)
        raise "Corridor boundary does not contain both constraint edge end nodes" unless bpos
        polygon1 = [a] + chain[1..bpos]
        polygon2 = [a] + chain[(bpos + 1)..].reverse + [b]

        new_triangles = []
        ear_triangulate(polygon1) {|tri| new_triangles << tri }
        ear_triangulate(polygon2) {|tri| new_triangles << tri }

        corridor.each {|tri| tri.alive = false }
        new_triangles.each {|tri| @tris << tri }
        rebuild_structure
        relink_triangles
      end

      # Builds the directed boundary loop starting at node a from the map of directed boundary
      # edges (each edge occurring exactly once).
      def boundary_loop(directed, a)
        outgoing = directed.keys.group_by(&:first)
        loop_nodes = [a]
        seen = {}
        loop do
          edges = outgoing[loop_nodes.last]
          if edges.nil? || edges.empty?
            return loop_nodes if loop_nodes.first == loop_nodes.last && loop_nodes.size > 1
            return nil
          end
          edge = edges.first
          return nil if seen[edge]
          seen[edge] = true
          next_node = edge.last
          loop_nodes << next_node
          return loop_nodes if next_node == a
          return nil if loop_nodes.size > @nodes.size + 1
        end
      end

      # Ear-clips the simple node polygon (reorienting it counterclockwise if needed), creating
      # Tri objects and yielding them.
      def ear_triangulate(polygon)
        polygon = polygon.reverse if polygon_area(polygon) < 0
        return if polygon.size < 3

        remaining = polygon.dup
        while remaining.size > 3
          ear_index = find_ear(remaining)
          raise "No ear found during corridor re-triangulation" unless ear_index
          nodes = three_consecutive(remaining, ear_index)
          yield make_corridor_triangle(nodes)
          remaining.delete_at(ear_index)
        end
        yield make_corridor_triangle(remaining)
      end

      # Returns the index of a valid ear of the counterclockwise polygon, or +nil+ if none exists.
      # Collinear vertices are never clipped (they become middle vertices of ears instead), which
      # guarantees that every boundary node is used.
      def find_ear(nodes)
        nodes.each_with_index do |_, i|
          p0 = nodes[(i - 1) % nodes.size]
          p1 = nodes[i]
          p2 = nodes[(i + 1) % nodes.size]
          next unless orient(p0, p1, p2) > 0
          empty = nodes.each_with_index.all? do |node, j|
            next true if j == (i - 1) % nodes.size || j == i || j == (i + 1) % nodes.size
            !strictly_inside_triangle?(node, p0, p1, p2)
          end
          return i if empty
        end
        nil
      end

      # Returns the three consecutive nodes starting at index i in counterclockwise order.
      def three_consecutive(nodes, i)
        [nodes[i], nodes[(i + 1) % nodes.size], nodes[(i - 1) % nodes.size]]
      end

      # Returns +true+ if node p is strictly inside the counterclockwise triangle a-b-c.
      def strictly_inside_triangle?(p, a, b, c)
        orient(a, b, p) > 0 && orient(b, c, p) > 0 && orient(c, a, p) > 0
      end

      # Creates a CCW triangle from the three given nodes.
      def make_corridor_triangle(nodes)
        Tri.new(nodes[0], nodes[1], nodes[2])
      end

      # Delaunay edge flipping

      # Flips flippable non-constraint edges until every such edge is locally Delaunay. The
      # in-circle test uses a strict comparison, so co-circular edges are never flipped and the
      # triangulation stays deterministic.
      def delaunay_flip
        loop do
          rebuild_structure
          relink_triangles
          flipped = flip_one_edge
          break unless flipped
        end
      end

      # Finds the first flippable, non-constraint edge that violates the empty circle condition
      # and flips it. Returns +true+ if an edge was flipped.
      def flip_one_edge
        @edge_map.each do |edge, entries|
          next if entries.size != 2 || @constraints.include?(edge)
          tri1, slot1 = entries[0]
          tri2, slot2 = entries[1]
          next unless tri1.alive && tri2.alive

          op1 = tri1.v[slot1]
          op2 = tri2.v[slot2]
          e0, e1 = edge
          o1 = orient(e0, e1, op1)
          o2 = orient(e0, e1, op2)
          next unless (o1.positive? && o2.negative?) || (o1.negative? && o2.positive?)
          # both diagonals have to lie inside the quadrilateral
          next unless (orient(op1, op2, e0) * orient(op1, op2, e1)).negative?

          if o1.positive? # op1 is left of e0->e1, op2 is right
            ccw_triangle = [e0, e1, op1]
            new_triangles = [[op1, e0, op2], [op1, op2, e1]]
          else
            ccw_triangle = [e1, e0, op1]
            new_triangles = [[op1, e1, op2], [op1, op2, e0]]
          end
          next unless incircle(ccw_triangle[0], ccw_triangle[1], ccw_triangle[2], op2).positive?

          do_flip(tri1, tri2, new_triangles)
          return true
        end
        false
      end

      # Replaces the two triangles sharing the edge e0-e1 by two triangles sharing the new
      # diagonal op1-op2.
      def do_flip(tri1, tri2, new_triangles)
        tri1.alive = false
        tri2.alive = false
        new_triangles.each do |vertices|
          @tris << Tri.new(*vertices)
        end
      end

      # Classification and output

      # Removes triangles incident to super-triangle nodes and triangles outside the planar
      # region (including the interiors of holes), and builds the public Triangulation result.
      def classify_and_output(_contours)
        rebuild_structure
        relink_triangles

        alive = @tris.select {|tri| tri.alive && (tri.v & @super_nodes).empty? }
        inside = alive.select {|tri| region_winding(centroid(tri)).nonzero? }

        used_nodes = inside.flat_map(&:v).uniq.sort
        node_map = {}
        used_nodes.each_with_index {|node, index| node_map[node] = index }
        tri_map = {}
        inside.each_with_index {|tri, index| tri_map[tri] = index }

        points = used_nodes.map {|node| Geom2D::Point.new(@nodes[node][0], @nodes[node][1]) }
        triangles = inside.map do |tri|
          vertices = tri.v.map {|node| node_map[node] }
          neighbors = tri.n.map do |neighbor|
            index = tri_map[neighbor]
            index unless index.nil? || neighbor.nil?
          end
          Triangle.new(vertices, neighbors)
        end
        constraints = @constraints.select {|a, b| node_map[a] && node_map[b] }.map do |a, b|
          [node_map[a], node_map[b]]
        end

        Triangulation.new(points, triangles, constraints)
      end

      # Returns the exact arithmetic centroid of the triangle as a pair of Rational values.
      def centroid(tri)
        p0 = @coord[tri.v[0]]
        p1 = @coord[tri.v[1]]
        p2 = @coord[tri.v[2]]
        [(p0[0] + p1[0] + p2[0]) / 3, (p0[1] + p1[1] + p2[1]) / 3]
      end

      # Returns the signed winding number of the point with respect to the whole planar region
      # using exact arithmetic. Counterclockwise outer contours contribute +1, clockwise holes
      # -1. Non-zero winding means the point is inside the region; this also handles nested
      # islands-in-holes correctly.
      def region_winding(point)
        @contours.sum do |contour|
          sign = contour[:hole] ? -1 : 1
          point_in_contour?(point, contour[:nodes]) ? sign : 0
        end
      end

      # Returns +true+ if the point lies inside (or on the boundary of) the closed contour given
      # by its node indices, using the exact ray-crossing test.
      def point_in_contour?(point, nodes)
        inside = false
        nodes.each_with_index do |a, i|
          b = nodes[(i + 1) % nodes.size]
          pa = @coord[a]
          pb = @coord[b]
          next unless (pa[1] > point[1]) != (pb[1] > point[1])
          cross_x = pa[0] + (point[1] - pa[1]) * (pb[0] - pa[0]) / (pb[1] - pa[1])
          inside = !inside if point[0] < cross_x
        end
        inside
      end

    end
  end
end
