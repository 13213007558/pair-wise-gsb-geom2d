# -*- coding: utf-8 -*-
# frozen_string_literal: true

require 'test_helper'
require 'geom2d'
require 'geom2d/algorithms'

CDT = Geom2D::Algorithms::ConstrainedDelaunayTriangulation

describe CDT do
  before do
    @square = Geom2D::Polygon([0, 0], [10, 0], [10, 10], [0, 10])
    @square_hole = Geom2D::Polygon([3, 3], [3, 7], [7, 7], [7, 3]) # clockwise
  end

  def orient(a, b, c)
    (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
  end

  def proper_cross?(a, b, c, d)
    o1 = orient(a, b, c)
    o2 = orient(a, b, d)
    o3 = orient(c, d, a)
    o4 = orient(c, d, b)
    o1 * o2 < -1e-9 && o3 * o4 < -1e-9
  end

  def assert_valid_mesh(result)
    pts = result.points.map(&:to_a)
    tris = result.triangles

    assert(tris.size >= 1)
    tris.each do |tri|
      a = pts[tri.vertices[0]]
      b = pts[tri.vertices[1]]
      c = pts[tri.vertices[2]]
      assert(orient(a, b, c) > 1e-9, 'all triangles must be counterclockwise with positive area')
      assert_equal(3, tri.vertices.uniq.size)
    end

    tris.each_with_index do |tri, i|
      tri.neighbors.each_with_index do |neighbor, slot|
        next if neighbor.nil?
        e0 = tri.vertices[(slot + 1) % 3]
        e1 = tri.vertices[(slot + 2) % 3]
        other = tris[neighbor]
        assert(other.vertices.include?(e0) && other.vertices.include?(e1),
               'neighbor must share the expected edge')
        assert(other.neighbors.include?(i), 'adjacency must be symmetric')
      end
    end

    tris.each_with_index do |t1, i|
      p1 = t1.vertices.map {|v| pts[v] }
      tris.each_with_index do |t2, j|
        next if i >= j
        p2 = t2.vertices.map {|v| pts[v] }
        next if (t1.vertices & t2.vertices).size >= 2
        edges = [[0, 1], [1, 2], [2, 0]]
        crossing = edges.any? do |e1|
          edges.any? {|e2| proper_cross?(p1[e1[0]], p1[e1[1]], p2[e2[0]], p2[e2[1]]) }
        end
        refute(crossing, 'triangle interiors must not overlap')
      end
    end
  end

  def incircle(pa, pb, pc, pd)
    adx = pa[0] - pd[0]; ady = pa[1] - pd[1]
    bdx = pb[0] - pd[0]; bdy = pb[1] - pd[1]
    cdx = pc[0] - pd[0]; cdy = pc[1] - pd[1]
    a2 = adx * adx + ady * ady
    b2 = bdx * bdx + bdy * bdy
    c2 = cdx * cdx + cdy * cdy
    a2 * (bdx * cdy - cdx * bdy) -
      b2 * (adx * cdy - cdx * ady) +
      c2 * (adx * bdy - bdx * ady)
  end

  def assert_locally_delaunay(result)
    pts = result.points.map(&:to_a)
    result.triangles.each_with_index do |tri, i|
      tri.neighbors.each_with_index do |neighbor, slot|
        next if neighbor.nil? || neighbor <= i
        e0 = tri.vertices[(slot + 1) % 3]
        e1 = tri.vertices[(slot + 2) % 3]
        next if result.constraint_edge?(e0, e1)
        op1 = tri.vertices[slot]
        op2 = result.triangles[neighbor].vertices.find {|v| !tri.vertices.include?(v) }
        ccw = (orient(pts[e0], pts[e1], pts[op1]) > 0 ? [e0, e1, op1] : [e1, e0, op1])
        value = incircle(pts[ccw[0]], pts[ccw[1]], pts[ccw[2]], pts[op2])
        assert(value <= 1e-6, 'non-constraint edge must satisfy the empty circle condition')
      end
    end
  end

  it 'triangulates a single polygon' do
    result = CDT.run(@square)
    assert_valid_mesh(result)
    assert_locally_delaunay(result)
    assert_in_delta(100.0, result.area, 1e-9)
    assert_equal(2, result.triangles.size)
  end

  it 'accepts a single Polygon as well as a PolygonSet' do
    from_polygon = CDT.run(@square)
    from_set = CDT.run(Geom2D::PolygonSet(@square))
    assert_equal(from_polygon.triangles.size, from_set.triangles.size)
  end

  it 'keeps holes empty and the total area correct' do
    result = CDT.run(Geom2D::PolygonSet(@square, @square_hole))
    assert_valid_mesh(result)
    assert_locally_delaunay(result)
    assert_in_delta(84.0, result.area, 1e-9)
    assert(result.triangles.size > 2)
    assert(result.constraint_edges.size >= 8)
    result.constraint_edges.each do |a, b|
      assert(result.triangles.any? {|t| t.vertices.include?(a) && t.vertices.include?(b) })
    end
  end

  it 'triangulates multiple disjoint outer contours' do
    p1 = Geom2D::Polygon([0, 0], [4, 0], [4, 4], [0, 4])
    p2 = Geom2D::Polygon([20, 20], [24, 20], [24, 24], [20, 24])
    result = CDT.run(Geom2D::PolygonSet(p1, p2))
    assert_valid_mesh(result)
    assert_in_delta(32.0, result.area, 1e-9)
  end

  it 'keeps collinear boundary vertices and their constraint edges' do
    polygon = Geom2D::Polygon([0, 0], [5, 0], [10, 0], [10, 10], [0, 10])
    result = CDT.run(polygon)
    assert_valid_mesh(result)
    assert(result.points.map(&:to_a).include?([5, 0]))
    assert_in_delta(100.0, result.area, 1e-9)
  end

  it 'preserves a very thin, narrow hole' do
    outer = Geom2D::Polygon([0, 0], [100, 0], [100, 100], [0, 100])
    thin = Geom2D::Polygon([10, 49.999], [10, 50.001], [90, 50.001], [90, 49.999])
    result = CDT.run(Geom2D::PolygonSet(outer, thin))
    assert_valid_mesh(result)
    assert_in_delta(10_000 - 0.16, result.area, 1e-6)
  end

  it 'triangulates a non-convex polygon' do
    polygon = Geom2D::Polygon([0, 0], [10, 0], [10, 10], [7, 10], [7, 3], [3, 3], [3, 10], [0, 10])
    result = CDT.run(polygon)
    assert_valid_mesh(result)
    assert_locally_delaunay(result)
    assert_in_delta(72.0, result.area, 1e-9)
  end

  it 'handles nested contours (island inside a hole)' do
    outer = Geom2D::Polygon([0, 0], [20, 0], [20, 20], [0, 20])
    hole = Geom2D::Polygon([4, 4], [4, 16], [16, 16], [16, 4])
    island = Geom2D::Polygon([8, 8], [12, 8], [12, 12], [8, 12])
    result = CDT.run(Geom2D::PolygonSet(outer, hole, island))
    assert_valid_mesh(result)
    assert_in_delta(400.0 - 144.0 + 16.0, result.area, 1e-9)
  end

  it 'is deterministic for co-circular input' do
    first = CDT.run(@square).triangles.map(&:vertices)
    second = CDT.run(@square).triangles.map(&:vertices)
    assert_equal(first, second)
  end

  it 'rejects properly crossing constraint edges with conflict details' do
    t1 = Geom2D::Polygon([0, 0], [10, 10], [0, 10])
    t2 = Geom2D::Polygon([0, 10], [10, 0], [10, 10])
    error = assert_raises(CDT::ConstraintError) do
      CDT.run(Geom2D::PolygonSet(t1, t2))
    end
    assert_equal(1, error.conflicts.size)
    conflict = error.conflicts.first
    assert_equal(2, conflict[:edges].size)
    assert_kind_of(Geom2D::Segment, conflict[:edges].first)
    assert_kind_of(Geom2D::Point, conflict[:point])
    assert_in_delta(5.0, conflict[:point].x, 1e-9)
    assert_in_delta(5.0, conflict[:point].y, 1e-9)
  end

  it 'rejects dangling partial collinear overlaps' do
    p1 = Geom2D::Polygon([0, 0], [8, 0], [8, -4], [0, -4])
    p2 = Geom2D::Polygon([2, 0], [10, 0], [10, 4], [2, 4])
    assert_raises(CDT::ConstraintError) { CDT.run(Geom2D::PolygonSet(p1, p2)) }
  end

  it 'allows a T-junction where one edge ends inside another edge' do
    outer = Geom2D::Polygon([0, 0], [10, 0], [10, 10], [0, 10])
    attached = Geom2D::Polygon([2, 0], [8, 0], [8, -5], [2, -5])
    result = CDT.run(Geom2D::PolygonSet(outer, attached))
    assert_valid_mesh(result)
    assert_in_delta(130.0, result.area, 1e-9)
  end

  it 'integrates with the boolean polygon operations' do
    subject = Geom2D::PolygonSet(Geom2D::Polygon([0, 0], [10, 0], [10, 10], [0, 10]))
    clipping = Geom2D::PolygonSet(Geom2D::Polygon([5, 5], [15, 5], [15, 15], [5, 15]))
    difference = Geom2D::Algorithms::PolygonOperation.run(subject, clipping, :difference)
    result = CDT.run(difference)
    assert_valid_mesh(result)
    assert_in_delta(75.0, result.area, 1e-9)
  end

  it 'rejects a polygon with fewer than three distinct vertices' do
    assert_raises(ArgumentError) { CDT.run(Geom2D::Polygon([0, 0], [1, 0])) }
  end

  it 'returns a useful inspection string' do
    result = CDT.run(@square)
    assert_match(/triangles/, result.inspect)
    assert_match(/points/, result.inspect)
  end
end
