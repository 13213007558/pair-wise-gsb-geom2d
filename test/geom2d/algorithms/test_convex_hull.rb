# -*- frozen_string_literal: true -*-

require 'test_helper'
require 'geom2d/algorithms/convex_hull'
require 'geom2d/point'

describe Geom2D::Algorithms::ConvexHull do
  def hull_indices(points)
    Geom2D::Algorithms::ConvexHull.hull_indices(points)
  end

  it "returns an empty hull for no points" do
    assert_equal([], hull_indices([]))
  end

  it "returns the point itself for a single point" do
    assert_equal([0], hull_indices([Geom2D::Point(1, 2)]))
  end

  it "returns both points for two points" do
    assert_equal([0, 1], hull_indices([Geom2D::Point(1, 2), Geom2D::Point(3, 4)]))
  end

  it "returns the vertices counterclockwise, ignoring interior points" do
    points = [
      Geom2D::Point(0, 0), Geom2D::Point(4, 0), Geom2D::Point(4, 4),
      Geom2D::Point(0, 4), Geom2D::Point(2, 2),
    ]
    assert_equal([0, 1, 2, 3], hull_indices(points))
  end

  it "excludes collinear points on hull edges" do
    points = [
      Geom2D::Point(0, 0), Geom2D::Point(2, 0), Geom2D::Point(4, 0),
      Geom2D::Point(4, 4), Geom2D::Point(0, 4),
    ]
    assert_equal([0, 2, 3, 4], hull_indices(points))
  end

  it "returns the two extreme points if all points are collinear" do
    points = [Geom2D::Point(0, 0), Geom2D::Point(1, 1), Geom2D::Point(2, 2)]
    assert_equal([0, 2], hull_indices(points))
  end

  it "collapses coincident points to the lowest index" do
    points = [
      Geom2D::Point(0, 0), Geom2D::Point(4, 0), Geom2D::Point(4, 4),
      Geom2D::Point(0, 0), Geom2D::Point(0, 4),
    ]
    assert_equal([0, 1, 2, 4], hull_indices(points))
  end

  it "returns a single index if all points coincide" do
    points = [Geom2D::Point(1, 1), Geom2D::Point(1, 1), Geom2D::Point(1, 1)]
    assert_equal([0], hull_indices(points))
  end

  it "returns points with ::run" do
    points = [Geom2D::Point(0, 0), Geom2D::Point(4, 0), Geom2D::Point(0, 4)]
    assert_equal(points, Geom2D::Algorithms::ConvexHull.run(points))
  end
end
