# -*- frozen_string_literal: true -*-

require 'test_helper'
require 'geom2d/algorithms/constrained_delaunay'

describe Geom2D::Algorithms::ConstrainedDelaunay do
  before do
    @klass = Geom2D::Algorithms::ConstrainedDelaunay
  end

  def triangulate(*polygons)
    @polygon_set = Geom2D::PolygonSet.new(polygons)
    @klass.run(@polygon_set)
  end

  def assert_valid(triangulation, polygon_set = @polygon_set)
    problems = triangulation.validate(polygon_set)
    assert(problems.empty?, "Triangulation should be valid but:\n#{problems.join("\n")}")
  end

  def assert_adjacency_consistent(triangulation)
    triangulation.neighbors.each_with_index do |neighbor_indices, ti|
      neighbor_indices.each_with_index do |ni, k|
        next if ni.nil?
        shared = triangulation.triangles[ti] & triangulation.triangles[ni]
        assert_equal(2, shared.size, "Neighbors #{ti} and #{ni} must share an edge")
        assert(triangulation.neighbors[ni].include?(ti),
               "Neighbor relation between #{ti} and #{ni} must be symmetric")
        edge = triangulation.triangles[ti] - shared
        assert_equal(1, edge.size)
        assert_equal(triangulation.triangles[ti][k], edge.first,
                     "Neighbor #{ni} must lie across the edge opposite vertex #{k}")
      end
    end
  end

  describe "simple shapes" do
    it "triangulates a square into two triangles with adjacency" do
      result = triangulate(Geom2D::Polygon([0, 0], [4, 0], [4, 4], [0, 4]))
      assert_equal(2, result.triangles.size)
      assert_in_delta(16.0, result.area)
      assert_equal(1, result.neighbors[0].compact.size)
      assert_equal(1, result.neighbors[1].compact.size)
      assert_adjacency_consistent(result)
      assert_valid(result)
    end

    it "returns an empty triangulation for an empty polygon set" do
      result = @klass.run(Geom2D::PolygonSet.new)
      assert_equal([], result.triangles)
      assert_equal([], result.neighbors)
      assert_equal(0.0, result.area)
      assert(result.valid?(Geom2D::PolygonSet.new))
    end

    it "supports multiple outer contours" do
      result = triangulate(Geom2D::Polygon([0, 0], [4, 0], [4, 4], [0, 4]),
                           Geom2D::Polygon([10, 0], [14, 0], [14, 4], [10, 4]))
      assert_equal(4, result.triangles.size)
      assert_in_delta(32.0, result.area)
      assert_adjacency_consistent(result)
      assert_valid(result)
    end

    it "supports contours touching at a shared vertex" do
      result = triangulate(Geom2D::Polygon([0, 0], [2, 0], [2, 2], [0, 2]),
                           Geom2D::Polygon([2, 2], [4, 2], [4, 4], [2, 4]))
      assert_in_delta(8.0, result.area)
      assert_valid(result)
    end

    it "handles contours with collinear consecutive vertices" do
      result = triangulate(Geom2D::Polygon([0, 0], [5, 0], [10, 0], [10, 10], [0, 10]))
      assert_in_delta(100.0, result.area)
      assert_valid(result)
    end

    it "raises an error for contours with less than three distinct vertices" do
      assert_raises(ArgumentError) do
        triangulate(Geom2D::Polygon([0, 0], [1, 1]))
      end
    end
  end

  describe "holes" do
    before do
      @outer = Geom2D::Polygon([0, 0], [10, 0], [10, 10], [0, 10])
    end

    it "keeps holes empty" do
      hole = Geom2D::Polygon([2, 2], [8, 2], [8, 8], [2, 8])
      result = triangulate(@outer, hole)
      assert_in_delta(64.0, result.area)
      result.each_triangle do |(a, b, c)|
        cx = (a.x + b.x + c.x) / 3.0
        cy = (a.y + b.y + c.y) / 3.0
        assert(cx < 2 || cx > 8 || cy < 2 || cy > 8,
               "No triangle centroid may lie inside the hole")
      end
      assert_adjacency_consistent(result)
      assert_valid(result)
    end

    it "does not depend on the hole orientation" do
      hole_cw = Geom2D::Polygon([2, 2], [2, 8], [8, 8], [8, 2])
      result = triangulate(@outer, hole_cw)
      assert_in_delta(64.0, result.area)
      assert_valid(result)
    end

    it "supports islands inside holes" do
      outer = Geom2D::Polygon([0, 0], [20, 0], [20, 20], [0, 20])
      hole = Geom2D::Polygon([4, 4], [16, 4], [16, 16], [4, 16])
      island = Geom2D::Polygon([8, 8], [12, 8], [12, 12], [8, 12])
      result = triangulate(outer, hole, island)
      assert_in_delta(400.0 - 144.0 + 16.0, result.area)
      assert_valid(result)
    end

    it "supports many holes" do
      holes = (0...5).map do |i|
        x = 5 + i * 18
        Geom2D::Polygon([x, 5], [x + 10, 5], [x + 10, 90], [x, 90])
      end
      result = triangulate(Geom2D::Polygon([0, 0], [100, 0], [100, 100], [0, 100]), *holes)
      assert_in_delta(10_000.0 - 5 * 850.0, result.area)
      assert_valid(result)
    end

    it "preserves narrow sliver holes" do
      outer = Geom2D::Polygon([0, 0], [100, 0], [100, 100], [0, 100])
      sliver = Geom2D::Polygon([10, 10], [90, 10], [90, 10.0001], [10, 10.0001])
      result = triangulate(outer, sliver)
      assert_in_delta(10_000.0 - 80 * 0.0001, result.area)
      assert_valid(result)
    end
  end

  describe "degenerate configurations" do
    it "handles cocircular points deterministically" do
      square = Geom2D::Polygon([0, 0], [1, 0], [1, 1], [0, 1])
      result1 = triangulate(square)
      result2 = @klass.run(@polygon_set)
      assert_equal(result1.triangles, result2.triangles)
      assert_valid(result1)
    end

    it "handles regular polygons whose vertices are all cocircular" do
      outer = Geom2D::Polygon(*(0...20).map do |i|
        angle = 2 * Math::PI * i / 20
        [50 * Math.cos(angle), 50 * Math.sin(angle)]
      end)
      hole = Geom2D::Polygon(*(0...10).map do |i|
        angle = 2 * Math::PI * i / 10
        [20 * Math.cos(angle), 20 * Math.sin(angle)]
      end)
      result = triangulate(outer, hole)
      assert_valid(result)
    end

    it "handles vertices collinear with constraint edges" do
      outer = Geom2D::Polygon([0, 0], [5, 0], [10, 0], [10, 5], [10, 10], [5, 10], [0, 10], [0, 5])
      hole = Geom2D::Polygon([2, 2], [4, 2], [6, 2], [8, 2], [8, 8], [2, 8])
      result = triangulate(outer, hole)
      assert_in_delta(100.0 - 36.0, result.area)
      assert_valid(result)
    end

    it "merges nearly identical points" do
      result = triangulate(Geom2D::Polygon([0, 0], [4, 0], [4, 4], [0, 4], [0, 1e-12]))
      assert_valid(result)
    end
  end

  describe "constraint conflicts" do
    it "rejects crossing constraint edges with a detailed error" do
      t1 = Geom2D::Polygon([0, 0], [4, 0], [2, 3])
      t2 = Geom2D::Polygon([0, 1], [4, 1], [2, -1])
      error = assert_raises(Geom2D::Algorithms::ConstrainedDelaunay::ConstraintConflictError) do
        triangulate(t1, t2)
      end
      assert_equal([0, 0], error.location1)
      assert_equal([1, 1], error.location2)
      assert_equal(Geom2D::Point(3, 0), error.intersection)
      assert_match(/contour 0/, error.message)
      assert_match(/contour 1/, error.message)
    end

    it "rejects self-intersecting contours" do
      bowtie = Geom2D::Polygon([0, 0], [4, 4], [4, 0], [0, 4])
      error = assert_raises(Geom2D::Algorithms::ConstrainedDelaunay::ConstraintConflictError) do
        triangulate(bowtie)
      end
      assert_equal(Geom2D::Point(2, 2), error.intersection)
    end

    it "rejects collinearly overlapping constraint edges" do
      c1 = Geom2D::Polygon([0, 0], [4, 0], [4, 4], [0, 4])
      c2 = Geom2D::Polygon([2, 0], [6, 0], [6, 2], [2, 2])
      error = assert_raises(Geom2D::Algorithms::ConstrainedDelaunay::ConstraintConflictError) do
        triangulate(c1, c2)
      end
      assert_kind_of(Geom2D::Segment, error.intersection)
    end

    it "rejects constraint vertices lying on other constraint edges" do
      outer = Geom2D::Polygon([0, 0], [10, 0], [10, 10], [0, 10])
      touching_hole = Geom2D::Polygon([5, 0], [7, 2], [3, 2])
      error = assert_raises(Geom2D::Algorithms::ConstrainedDelaunay::ConstraintConflictError) do
        triangulate(outer, touching_hole)
      end
      assert_equal(Geom2D::Point(5, 0), error.intersection)
    end
  end

  describe "validation" do
    it "detects area mismatches" do
      result = triangulate(Geom2D::Polygon([0, 0], [4, 0], [4, 4], [0, 4]))
      other = Geom2D::PolygonSet.new([Geom2D::Polygon([0, 0], [5, 0], [5, 5], [0, 5])])
      problems = result.validate(other)
      assert(problems.any? {|p| p.include?("area") })
    end

    it "detects triangles inside holes" do
      result = triangulate(Geom2D::Polygon([0, 0], [10, 0], [10, 10], [0, 10]))
      with_hole = Geom2D::PolygonSet.new([Geom2D::Polygon([0, 0], [10, 0], [10, 10], [0, 10]),
                                          Geom2D::Polygon([2, 2], [8, 2], [8, 8], [2, 8])])
      problems = result.validate(with_hole)
      assert(problems.any? {|p| p.include?("hole") })
    end

    it "detects overlapping triangles" do
      vertices = [Geom2D::Point(0, 0), Geom2D::Point(4, 0), Geom2D::Point(4, 4),
                  Geom2D::Point(0, 4), Geom2D::Point(2, 1)]
      triangulation = Geom2D::Algorithms::ConstrainedDelaunay::Triangulation.new(
        vertices, [[0, 1, 2], [0, 2, 3], [0, 1, 4]],
        [[nil, nil, nil], [nil, nil, nil], [nil, nil, nil]], []
      )
      problems = triangulation.validate(Geom2D::PolygonSet.new)
      assert(problems.any? {|p| p.include?("overlap") })
    end

    it "detects missing constraint edges" do
      result = triangulate(Geom2D::Polygon([0, 0], [4, 0], [4, 4], [0, 4]))
      other = Geom2D::PolygonSet.new([Geom2D::Polygon([0, 0], [4, 0], [4, 4], [0, 4]),
                                      Geom2D::Polygon([1, 1], [3, 1], [3, 3], [1, 3])])
      problems = result.validate(other)
      assert(problems.any? {|p| p.include?("not preserved") })
    end
  end

  describe "result structure" do
    it "marks constraint edges" do
      result = triangulate(Geom2D::Polygon([0, 0], [4, 0], [4, 4], [0, 4]))
      assert_equal(4, result.constraint_edges.size)
      result.constraint_edges.each do |(i, j)|
        assert(result.constrained_edge?(i, j))
        assert(result.constrained_edge?(j, i))
      end
    end

    it "yields counterclockwise triangles" do
      result = triangulate(Geom2D::Polygon([0, 0], [10, 0], [10, 10], [0, 10]),
                           Geom2D::Polygon([2, 2], [8, 2], [8, 8], [2, 8]))
      result.each_triangle do |(a, b, c)|
        cross = (b.x - a.x) * (c.y - a.y) - (c.x - a.x) * (b.y - a.y)
        assert(cross > 0, "Triangle #{[a, b, c]} must be counterclockwise")
      end
    end
  end
end
