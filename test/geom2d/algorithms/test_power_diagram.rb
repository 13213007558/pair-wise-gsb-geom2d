# -*- frozen_string_literal: true -*-

require 'test_helper'
require 'geom2d/algorithms/power_diagram'

describe Geom2D::Algorithms::PowerDiagram do
  before do
    @domain = Geom2D::Polygon([0, 0], [1, 0], [1, 1], [0, 1])
    @sites = [[0.25, 0.25], [0.75, 0.25], [0.75, 0.75], [0.25, 0.75]]
  end

  it 'restores equal areas for a symmetric layout' do
    result = Geom2D::Algorithms::PowerDiagram.solve(@domain, @sites, [0.25] * 4)

    assert_equal(:optimal, result.status)
    assert_equal(0, result.iterations)
    assert_in_delta(0.0, result.weights.sum, 1e-12)
    result.areas.each {|area| assert_in_delta(0.25, area, 1e-10) }
    assert_in_delta(0.0, result.coverage_residual, 1e-12)
  end

  it 'moves boundaries for unequal capacity targets' do
    result = Geom2D::Algorithms::PowerDiagram.solve(@domain, @sites, [0.7, 0.1, 0.1, 0.1])

    assert_equal(:optimal, result.status)
    assert(result.weights.first > result.weights.last)
    assert_in_delta(0.7, result.areas.first, 5e-7)
    result.areas[1..].each {|area| assert_in_delta(0.1, area, 5e-7) }
    assert_in_delta(0.0, result.area_residual_norm, 5e-7)
    assert_in_delta(1.0, result.areas.sum, 1e-12)
  end

  it 'covers the original polygon after clipping all cells' do
    result = Geom2D::Algorithms::PowerDiagram.solve(@domain, @sites, [0.7, 0.1, 0.1, 0.1])
    union = result.cells.reduce(Geom2D::PolygonSet.new) do |combined, cell|
      Geom2D::Algorithms::PolygonOperation.run(combined, cell, :union)
    end

    assert_in_delta(1.0, result.areas.sum, 1e-10)
    assert_in_delta(1.0, union.polygons.sum {|polygon| polygon_area(polygon) }, 1e-10)
  end

  it 'supports intentionally empty cells and reports their residual' do
    result = Geom2D::Algorithms::PowerDiagram.solve(@domain, @sites, [0.5, 0.5, 0.0, 0.0])

    assert_equal(:optimal, result.status)
    assert_in_delta(0.0, result.areas[2], 1e-10)
    assert_in_delta(0.0, result.areas[3], 1e-10)
    assert_in_delta(0.0, result.area_residuals[2], 1e-10)
    assert_in_delta(1.0, result.areas.sum, 1e-10)
  end

  it 'rejects duplicate sites with positive target areas' do
    sites = [[0.25, 0.25], [0.75, 0.25], [0.75, 0.25], [0.25, 0.75]]
    error = assert_raises(Geom2D::Algorithms::PowerDiagram::DuplicateSitesError) do
      Geom2D::Algorithms::PowerDiagram.solve(@domain, sites, [0.4, 0.3, 0.3, 0.0])
    end
    assert_match(/coincide/, error.message)
  end

  it 'allows coincident sites when every duplicate except one has zero target' do
    sites = [[0.25, 0.25], [0.75, 0.25], [0.75, 0.25], [0.25, 0.75]]
    result = Geom2D::Algorithms::PowerDiagram.solve(@domain, sites, [0.7, 0.3, 0.0, 0.0])

    assert_equal(:optimal, result.status)
    assert_in_delta(0.0, result.areas[2], 1e-10)
  end

  it 'rejects target area totals that do not match the domain area' do
    assert_raises(Geom2D::Algorithms::PowerDiagram::InfeasibleTargetError) do
      Geom2D::Algorithms::PowerDiagram.solve(@domain, @sites, [0.4] * 4)
    end
  end

  it 'reports iteration exhaustion with residuals instead of raising' do
    result = Geom2D::Algorithms::PowerDiagram.solve(
      @domain, @sites, [0.7, 0.1, 0.1, 0.1], max_iterations: 0,
      tolerance: 1e-15
    )

    assert_equal(:iteration_exhausted, result.status)
    assert(result.area_residual_norm > 1e-15)
    assert_raises(Geom2D::Algorithms::PowerDiagram::IterationExhausted) do
      Geom2D::Algorithms::PowerDiagram.solve!(
        @domain, @sites, [0.7, 0.1, 0.1, 0.1], max_iterations: 0,
        tolerance: 1e-15
      )
    end
  end

  it 'works on a concave bounded polygon' do
    domain = Geom2D::Polygon([0, 0], [2, 0], [2, 1], [1, 1], [1, 2], [0, 2])
    sites = [[0.4, 0.4], [1.5, 0.5], [0.5, 1.5]]
    result = Geom2D::Algorithms::PowerDiagram.solve(domain, sites, [1.5, 0.75, 0.75])

    assert_equal(:optimal, result.status)
    assert_in_delta(1.5, result.areas.first, 1e-6)
    assert_in_delta(0.75, result.areas[1], 1e-6)
    assert_in_delta(0.75, result.areas[2], 1e-6)
    assert_in_delta(3.0, result.areas.sum, 1e-9)
  end

  private

  def polygon_area(polygon)
    vertices = polygon.to_a
    signed_area = vertices.each_with_index.sum do |point, i|
      other = polygon.to_a[(i + 1) % polygon.nr_of_vertices]
      (other.x - point.x) * (other.y + point.y)
    end / 2.0
    signed_area.abs
  end
end
