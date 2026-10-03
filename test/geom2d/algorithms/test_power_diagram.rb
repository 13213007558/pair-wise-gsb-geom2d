# -*- frozen_string_literal: true -*-

require 'test_helper'
require 'geom2d/algorithms/power_diagram'

describe Geom2D::Algorithms::PowerDiagram do
  before do
    @pd = Geom2D::Algorithms::PowerDiagram
    @square = Geom2D::Polygon([0, 0], [10, 0], [10, 10], [0, 10])
    @sites = [[2.5, 2.5], [7.5, 2.5], [7.5, 7.5], [2.5, 7.5]]
  end

  def cell_area(cell)
    vertices = cell.to_a
    return 0.0 if vertices.size < 3
    area = vertices[-1].wedge(vertices[0])
    0.upto(vertices.size - 2) {|i| area += vertices[i].wedge(vertices[i + 1]) }
    area.abs * 0.5
  end

  describe '.cells' do
    it 'partitions the polygon into the ordinary Voronoi diagram without weights' do
      cells = @pd.cells(@square, @sites)
      assert_equal(4, cells.size)
      cells.each {|cell| assert_in_delta(25.0, cell_area(cell), 1e-8) }
      assert_in_delta(100.0, cells.sum {|cell| cell_area(cell) }, 1e-8)
    end

    it 'moves the cell boundaries when weights are given' do
      cells = @pd.cells(@square, @sites, [10, 0, 0, 0])
      assert_operator(cell_area(cells[0]), :>, 25.0)
      assert_in_delta(100.0, cells.sum {|cell| cell_area(cell) }, 1e-8)
    end

    it 'returns an empty polygon for an empty cell' do
      cells = @pd.cells(@square, [[5, 5], [50, 50]])
      assert_equal(0, cells[1].nr_of_vertices)
    end
  end

  describe '.fit with a symmetric layout' do
    it 'recovers the equal partition for equal target areas' do
      result = @pd.fit(@square, @sites, [25, 25, 25, 25])

      assert(result.converged?)
      assert_equal(:converged, result.status)
      result.areas.each {|area| assert_in_delta(25.0, area, 1e-6) }
      result.residuals.each {|residual| assert_in_delta(0.0, residual, 1e-6) }
    end

    it 'fixes the additive degree of freedom by centering the weights' do
      result = @pd.fit(@square, @sites, [25, 25, 25, 25])

      assert_in_delta(0.0, result.weights.sum, 1e-8)
      result.weights.each {|weight| assert_in_delta(result.weights[0], weight, 1e-6) }
    end
  end

  describe '.fit with unequal capacities' do
    before do
      @targets = [50, 20, 20, 10]
      @result = @pd.fit(@square, @sites, @targets)
    end

    it 'converges to the target areas' do
      assert(@result.converged?)
      @targets.each_with_index do |target, i|
        assert_in_delta(target, @result.areas[i], 1e-6)
        assert_in_delta(0.0, @result.residuals[i], 1e-6)
      end
    end

    it 'actually moves the boundaries away from the equal partition' do
      equal_areas = @pd.cells(@square, @sites).map {|cell| cell_area(cell) }
      assert_in_delta(25.0, equal_areas[0], 1e-8)
      refute_in_delta(equal_areas[0], @result.areas[0], 5.0)
      assert(@result.weights.uniq.size > 1)
    end

    it 'keeps the cells covering the whole polygon without overlaps' do
      cells = @pd.cells(@square, @sites, @result.weights)
      assert_in_delta(100.0, cells.sum {|cell| cell_area(cell) }, 1e-6)

      # No cell crosses the radical line of another cell, so the cells are disjoint and
      # their union covers the polygon.
      points = @sites.map {|site| Geom2D::Point(site) }
      cells.each_with_index do |cell, i|
        cell.each_vertex do |vertex|
          own = (vertex - points[i]).dot(vertex - points[i]) - @result.weights[i]
          points.each_with_index do |site, j|
            other = (vertex - site).dot(vertex - site) - @result.weights[j]
            assert_operator(own, :<=, other + 1e-8)
          end
        end
      end
    end
  end

  describe '.fit input validation' do
    it 'rejects target areas that do not sum up to the polygon area' do
      error = assert_raises(ArgumentError) do
        @pd.fit(@square, @sites, [25, 25, 25, 20])
      end
      assert_match(/sum up/, error.message)
    end

    it 'rejects negative target areas' do
      assert_raises(ArgumentError) { @pd.fit(@square, @sites, [125, -25, 0, 0]) }
    end

    it 'rejects a mismatched number of target areas' do
      assert_raises(ArgumentError) { @pd.fit(@square, @sites, [100]) }
    end

    it 'rejects an empty list of sites' do
      assert_raises(ArgumentError) { @pd.fit(@square, [], []) }
    end
  end

  describe '.fit failure modes' do
    it 'detects unreachable targets caused by coincident sites' do
      result = @pd.fit(@square, [[5, 5], [5, 5]], [50, 50])

      refute(result.converged?)
      assert_equal(:coincident_sites, result.status)
      assert_equal([50.0, 50.0], result.residuals)
    end

    it 'detects unreachable targets caused by empty cells' do
      result = @pd.fit(@square, [[5, 5], [50, 50]], [50, 50])

      refute(result.converged?)
      assert_equal(:empty_cell, result.status)
      assert_in_delta(0.0, result.areas[1], 1e-10)
    end

    it 'reports exhausted iterations and still provides the area residuals' do
      result = @pd.fit(@square, @sites, [50, 20, 20, 10], max_iterations: 0)

      refute(result.converged?)
      assert_equal(:max_iterations, result.status)
      assert_equal(4, result.residuals.size)
      assert_in_delta(100.0, result.areas.sum, 1e-8)
      result.residuals.each_with_index do |residual, i|
        assert_in_delta(result.areas[i] - [50, 20, 20, 10][i], residual, 1e-10)
      end
    end
  end

  it 'handles a single site covering the whole polygon' do
    result = @pd.fit(@square, [[5, 5]], [100])

    assert(result.converged?)
    assert_equal([0.0], result.weights)
    assert_in_delta(100.0, result.areas[0], 1e-8)
  end
end
