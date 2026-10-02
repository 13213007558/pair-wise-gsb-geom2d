# -*- frozen_string_literal: true -*-

require 'test_helper'
require 'geom2d/algorithms'
require 'geom2d/point'

describe Geom2D::Algorithms do
  describe "convex_hull" do
    it "returns vertices in counter-clockwise order" do
      points = [Geom2D.Point(0, 1), Geom2D.Point(1, 1), Geom2D.Point(1, 0),
                Geom2D.Point(0, 0)]
      hull = Geom2D::Algorithms.convex_hull(points)
      assert_equal([[0, 0], [1, 0], [1, 1], [0, 1]], hull.map(&:to_a))
    end

    it "omits points inside hull edges" do
      points = [Geom2D.Point(0, 0), Geom2D.Point(1, 0), Geom2D.Point(2, 0),
                Geom2D.Point(1, 1)]
      hull = Geom2D::Algorithms.convex_hull(points)
      assert_equal([[0, 0], [2, 0], [1, 1]], hull.map(&:to_a))
    end

    it "represents coincident points only once" do
      points = [Geom2D.Point(0, 0), Geom2D.Point(0, 0), Geom2D.Point(1, 0)]
      hull = Geom2D::Algorithms.convex_hull(points)
      assert_equal([[0, 0], [1, 0]], hull.map(&:to_a))
    end
  end

  describe "ccw" do
    before do
      @p1 = Geom2D::Point(1, 2)
      @p2 = Geom2D::Point(3, 4)
    end

    it "returns 1 for counterclockwise turn" do
      assert_equal(1, Geom2D::Algorithms.ccw(@p1, @p2, Geom2D::Point(0, 3)))
    end

    it "returns -1 for clockwise turn" do
      assert_equal(-1, Geom2D::Algorithms.ccw(@p1, @p2, Geom2D::Point(2, 0)))
    end

    it "returns 0 for collinear points" do
      assert_equal(0, Geom2D::Algorithms.ccw(@p1, @p2, @p2 + (@p2 - @p1)))
    end
  end

  describe "moving_convex_hull" do
    MovingPoint = Geom2D::Algorithms::MovingConvexHull::MovingPoint

    def time_values(result)
      result[:events].map {|event| event[:time].to_f }
    end

    it "reports entry and exit events with interval hulls" do
      points = [
        MovingPoint.new([0, 0], [0, 0]),
        MovingPoint.new([2, 0], [0, 0]),
        MovingPoint.new([2, 2], [0, 0]),
        MovingPoint.new([0, 2], [0, 0]),
        MovingPoint.new([-1, 2], [10_000, 0], id: :moving)
      ]
      result = Geom2D::Algorithms.moving_convex_hull(points, 0, 2)

      assert_in_delta(0.0001, time_values(result)[0], 1e-12)
      assert_in_delta(0.0003, time_values(result)[1], 1e-12)
      assert_includes(result[:intervals][0][:hull_vertex_ids], :moving)
      refute_includes(result[:intervals][1][:hull_vertex_ids], :moving)
      assert_includes(result[:intervals][2][:hull_vertex_ids], :moving)
      assert_equal([true] * 3, result[:intervals].map {|interval| interval[:open] })
    end

    it "groups simultaneous events and shows all intervening hull states" do
      points = [
        MovingPoint.new([0, 0], [0, 0]),
        MovingPoint.new([2, 0], [0, 0]),
        MovingPoint.new([2, 2], [0, 0]),
        MovingPoint.new([0, 2], [0, 0]),
        MovingPoint.new([0, 1], [0, 1], id: :left),
        MovingPoint.new([2, 1], [0, -1], id: :right)
      ]
      result = Geom2D::Algorithms.moving_convex_hull(points, 0, 2)

      assert_equal([1], time_values(result))
      assert_equal(1, result[:events].length)
      assert_includes(result[:events][0][:triples], [0, 3, 4])
      assert_includes(result[:events][0][:triples], [1, 2, 5])
      assert_includes(result[:intervals][1][:hull_vertex_ids], :left)
      assert_includes(result[:intervals][1][:hull_vertex_ids], :right)
      assert_equal([0, 1, 2, 3], result[:intervals][0][:hull_vertex_indices])
      assert_equal([0, 5, 2, 4], result[:intervals][1][:hull_vertex_indices])
    end

    it "keeps a hull vertex that exists only in a very short interval" do
      points = [
        MovingPoint.new([-5, 0], [-20, 1]),
        MovingPoint.new([2, 3], [-10, 1], id: :brief),
        MovingPoint.new([-2, -8], [0, 9]),
        MovingPoint.new([-1, 6], [6, -3]),
        MovingPoint.new([-6, -6], [-19, 6]),
        MovingPoint.new([-3, 0], [13, 9])
      ]
      result = Geom2D::Algorithms.moving_convex_hull(points, 0, 1)
      brief = result[:intervals].min_by do |interval|
        (interval[:t1] - interval[:t0]).to_f
      end

      refute_nil(brief)
      assert_operator(brief[:t1] - brief[:t0], :<, 0.01)
      assert_includes(brief[:hull_vertex_ids], :brief)
    end

    it "groups points coincident throughout the interval" do
      points = [
        MovingPoint.new([0, 0], [0, 0], id: :a),
        MovingPoint.new([0, 0], [0, 0], id: :b),
        MovingPoint.new([1, 0], [0, 0], id: :c)
      ]
      result = Geom2D::Algorithms.moving_convex_hull(points, 0, 1)

      assert_empty(result[:events])
      group_ids = result[:coincident_groups].map do |group|
        group.map {|index| result[:point_ids][index] }
      end
      assert_equal([%i[a b]], group_ids)
      assert_equal([:a, :c], result[:intervals][0][:hull_vertex_ids])
      assert_equal([%i[a b], [:c]], result[:intervals][0][:hull_vertex_groups])
    end

    it "reports a two-point encounter at the exact interior time" do
      points = [
        MovingPoint.new([0, 0], [1, 0], id: :a),
        MovingPoint.new([2, 0], [-1, 0], id: :b)
      ]
      result = Geom2D::Algorithms.moving_convex_hull(points, 0, 2)

      assert_equal([1], time_values(result))
      assert_equal([[:a, :b], [:a, :b]],
                   result[:intervals].map {|interval| interval[:hull_vertex_ids] })
      assert_equal([[0, 1]], result[:events][0][:pairs])
    end

    it "reports an instant coincidence at an existing hull vertex" do
      points = [
        MovingPoint.new([0, 0], [0, 0], id: :corner),
        MovingPoint.new([1, 0], [0, 0]),
        MovingPoint.new([1, 1], [0, 0]),
        MovingPoint.new([0, 1], [0, 0]),
        MovingPoint.new([-1, 0], [1, 0], id: :touch)
      ]
      result = Geom2D::Algorithms.moving_convex_hull(points, 0, 2)

      assert_equal([1], time_values(result))
      assert_equal([[0, 4]], result[:events][0][:pairs])
      assert_equal([1, 2, 3, 4], result[:intervals][0][:hull_vertex_indices])
      assert_equal([0, 1, 2, 3], result[:intervals][1][:hull_vertex_indices])
    end

    it "matches a static hull for arbitrary times inside each interval" do
      rng = Random.new(42)
      20.times do
        points = 8.times.map do |index|
          MovingPoint.new([rng.rand(-4.0..4.0), rng.rand(-4.0..4.0)],
                          [rng.rand(-4.0..4.0), rng.rand(-4.0..4.0)], id: index)
        end
        result = Geom2D::Algorithms.moving_convex_hull(points, 0, 1)

        result[:intervals].each do |interval|
          t = interval[:t0] + (interval[:t1] - interval[:t0]) * 0.37
          static_hull = Geom2D::Algorithms::MovingConvexHull
                       .hull_vertex_indices_at(points, t)
          assert_equal(static_hull, interval[:hull_vertex_indices])
        end
      end
    end
  end
end
