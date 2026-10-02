# -*- frozen_string_literal: true -*-

require 'test_helper'
require 'geom2d/algorithms/kinetic_convex_hull'
require 'geom2d/algorithms/convex_hull'
require 'geom2d/point'

describe Geom2D::Algorithms::KineticConvexHull do
  def analyze(points, velocities, t_start, t_end)
    Geom2D::Algorithms::KineticConvexHull.run(points, velocities, t_start, t_end)
  end

  def static_hull(points, velocities, t)
    moved = points.each_with_index.map do |point, i|
      Geom2D::Point(point.x + velocities[i].x * t, point.y + velocities[i].y * t)
    end
    Geom2D::Algorithms::ConvexHull.hull_indices(moved)
  end

  it "raises if the number of points and velocities differs" do
    assert_raises(ArgumentError) do
      Geom2D::Algorithms::KineticConvexHull.new([Geom2D::Point(0, 0)], [])
    end
  end

  it "raises if the time interval is invalid" do
    assert_raises(ArgumentError) { analyze([Geom2D::Point(0, 0)], [Geom2D::Point(0, 0)], 2, 1) }
  end

  it "reports no events for non-moving points" do
    points = [Geom2D::Point(0, 0), Geom2D::Point(4, 0), Geom2D::Point(2, 3)]
    velocities = [Geom2D::Point(0, 0), Geom2D::Point(0, 0), Geom2D::Point(0, 0)]
    result = analyze(points, velocities, 0, 5)
    assert_equal([], result.events)
    assert_equal(1, result.intervals.length)
    assert_equal(0, result.intervals[0].from)
    assert_equal(5, result.intervals[0].to)
    assert_equal([0, 1, 2], result.intervals[0].vertices)
  end

  it "detects a point entering the hull" do
    points = [Geom2D::Point(0, 0), Geom2D::Point(4, 0), Geom2D::Point(2, 3), Geom2D::Point(2, 1)]
    velocities = [Geom2D::Point(0, 0), Geom2D::Point(0, 0), Geom2D::Point(0, 0), Geom2D::Point(0, -1)]
    result = analyze(points, velocities, 0, 3)
    assert_equal(1, result.events.length)
    assert_in_delta(1.0, result.events[0], 1e-9)
    assert_equal([0, 1, 2], result.intervals[0].vertices)
    assert_equal([0, 3, 1, 2], result.intervals[1].vertices)
  end

  it "does not miss a vertex that is on the hull only briefly" do
    # D moves straight up through the triangle: it is a hull vertex below the triangle,
    # interior between t=1 and t=4 (where C is a vertex) and replaces C afterwards.
    points = [Geom2D::Point(0, 0), Geom2D::Point(4, 0), Geom2D::Point(2, 3), Geom2D::Point(2, -1)]
    velocities = [Geom2D::Point(0, 0), Geom2D::Point(0, 0), Geom2D::Point(0, 0), Geom2D::Point(0, 1)]
    result = analyze(points, velocities, 0, 5)
    assert_equal(2, result.events.length)
    assert_in_delta(1.0, result.events[0], 1e-9)
    assert_in_delta(4.0, result.events[1], 1e-9)
    assert_equal([0, 3, 1, 2], result.intervals[0].vertices)
    assert_equal([0, 1, 2], result.intervals[1].vertices)
    assert_equal([0, 1, 3], result.intervals[2].vertices)
  end

  it "merges simultaneous events into a single event time" do
    points = [
      Geom2D::Point(0, 0), Geom2D::Point(4, 0), Geom2D::Point(2, 3),
      Geom2D::Point(1, -1), Geom2D::Point(3, -1),
    ]
    velocities = [
      Geom2D::Point(0, 0), Geom2D::Point(0, 0), Geom2D::Point(0, 0),
      Geom2D::Point(0, 1), Geom2D::Point(0, 1),
    ]
    result = analyze(points, velocities, 0, 3)
    assert_equal(1, result.events.length)
    assert_in_delta(1.0, result.events[0], 1e-9)
    assert_equal([0, 3, 4, 1, 2], result.intervals[0].vertices)
    assert_equal([0, 1, 2], result.intervals[1].vertices)
  end

  it "handles permanently coincident points" do
    points = [Geom2D::Point(0, 0), Geom2D::Point(4, 0), Geom2D::Point(2, 3), Geom2D::Point(0, 0)]
    velocities = [Geom2D::Point(0, 0), Geom2D::Point(0, 0), Geom2D::Point(0, 0), Geom2D::Point(0, 0)]
    result = analyze(points, velocities, 0, 2)
    assert_equal([], result.events)
    assert_equal([0, 1, 2], result.intervals[0].vertices)
  end

  it "ignores collinearities that do not change the hull" do
    # E moves through the square's interior and is collinear with the diagonal A-C at t=1
    # without ever becoming a hull vertex.
    points = [
      Geom2D::Point(0, 0), Geom2D::Point(4, 0), Geom2D::Point(4, 4),
      Geom2D::Point(0, 4), Geom2D::Point(2, 1),
    ]
    velocities = [
      Geom2D::Point(0, 0), Geom2D::Point(0, 0), Geom2D::Point(0, 0),
      Geom2D::Point(0, 0), Geom2D::Point(0, 1),
    ]
    result = analyze(points, velocities, 0, 2)
    assert_equal([], result.events)
    assert_equal(1, result.intervals.length)
    assert_equal([0, 1, 2, 3], result.intervals[0].vertices)
  end

  it "accepts arrays as velocities" do
    points = [Geom2D::Point(0, 0), Geom2D::Point(4, 0), Geom2D::Point(2, 3), Geom2D::Point(2, 1)]
    result = analyze(points, [[0, 0], [0, 0], [0, 0], [0, -1]], 0, 3)
    assert_equal(1, result.events.length)
    assert_in_delta(1.0, result.events[0], 1e-9)
  end

  it "produces intervals whose hulls match the static hull at any sampled time" do
    random = Random.new(12345)
    5.times do
      count = 6
      points = Array.new(count) { Geom2D::Point(random.rand(-5.0..5.0), random.rand(-5.0..5.0)) }
      velocities = Array.new(count) { Geom2D::Point(random.rand(-1.0..1.0), random.rand(-1.0..1.0)) }
      result = analyze(points, velocities, 0, 4)

      # The intervals cover [0, 4] contiguously and events match the interval boundaries
      assert_equal(0, result.intervals.first.from)
      assert_equal(4, result.intervals.last.to)
      result.intervals.each_cons(2) {|a, b| assert_equal(a.to, b.from) }
      assert_equal(result.intervals[1..].map(&:from), result.events)

      # Any time sampled inside an interval yields the interval's hull
      result.intervals.each do |interval|
        [0.05, 0.37, 0.5, 0.83, 0.95].each do |fraction|
          t = interval.from + (interval.to - interval.from) * fraction
          assert_equal(interval.vertices, static_hull(points, velocities, t),
                       "hull mismatch at t=#{t} in interval (#{interval.from}, #{interval.to})")
        end
      end
    end
  end
end
