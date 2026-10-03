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
require 'geom2d/polygon'

module Geom2D
  module Algorithms

    # Computes power diagrams (also known as Laguerre-Voronoi diagrams) of a set of sites
    # clipped to a bounded polygon and fits the power weights so that the clipped cells match
    # given target areas.
    #
    # The power distance of a point +p+ to a site +s+ with weight +w+ is
    # |p - s|^2 - w. The cell of a site consists of all points of the polygon for which its
    # power distance is minimal. Adding a constant to all weights does not change the diagram,
    # so PowerDiagram.fit fixes this degree of freedom by keeping the weights centered on zero.
    class PowerDiagram

      include Utils

      # The result of a PowerDiagram.fit call.
      class Result

        # The fitted power weights, one per site, centered so that they sum to zero.
        attr_reader :weights

        # The areas of the clipped cells for the fitted weights.
        attr_reader :areas

        # The signed area residuals, i.e. areas[i] - target_areas[i] for each site.
        attr_reader :residuals

        # The number of iterations that were performed.
        attr_reader :iterations

        # The status of the fitting process, one of:
        #
        # [:converged]        All area residuals are within the requested tolerance.
        # [:max_iterations]   The iteration limit was exhausted before convergence.
        # [:empty_cell]       A cell became empty although its target area is positive, so the
        #                     target is unreachable.
        # [:coincident_sites] Two or more sites with positive target areas coincide, making
        #                     the targets unreachable.
        attr_reader :status

        def initialize(weights:, areas:, residuals:, iterations:, status:) # :nodoc:
          @weights = weights
          @areas = areas
          @residuals = residuals
          @iterations = iterations
          @status = status
        end

        # Returns +true+ if the fitting converged within the requested tolerance.
        def converged?
          @status == :converged
        end

      end

      class << self

        include Utils

        # Returns the cells of the power diagram of the given +sites+ (point-like objects)
        # with the given +weights+, clipped to the bounded +polygon+.
        #
        # If +weights+ is +nil+, all weights default to zero and the result is the ordinary
        # Voronoi diagram clipped to the polygon.
        #
        # The returned array contains one Geom2D::Polygon per site, in the same order. A
        # polygon has no vertices if the corresponding cell is empty.
        def cells(polygon, sites, weights = nil)
          sites = sites.map {|site| Geom2D::Point(site) }
          weights ||= Array.new(sites.size, 0.0)
          if weights.size != sites.size
            raise ArgumentError, 'Number of weights has to match number of sites'
          end

          sites.each_index.map do |i|
            cell = polygon.to_a
            sites.each_index do |j|
              next if i == j
              cell = clip_half_plane(cell, sites[i], weights[i], sites[j], weights[j])
              break if cell.empty?
            end
            Geom2D::Polygon.new(cell)
          end
        end

        # Fits power weights for the given +sites+ so that their cells clipped to the bounded
        # +polygon+ approximate the given +target_areas+.
        #
        # The target areas have to be non-negative and sum up to the area of the polygon,
        # otherwise an ArgumentError is raised.
        #
        # The optional arguments are:
        #
        # [+max_iterations+] The maximum number of Newton iterations (default: 100).
        # [+tolerance+]      The fitting stops once every absolute area residual is at most
        #                    this value (default: 1e-8).
        #
        # Returns a Result object containing the fitted weights, the achieved areas, the
        # signed area residuals, the number of performed iterations and a status. The weights
        # are centered on zero to fix the additive constant degree of freedom.
        def fit(polygon, sites, target_areas, max_iterations: 100, tolerance: 1e-8)
          sites = sites.map {|site| Geom2D::Point(site) }
          targets = target_areas.map {|area| Float(area) }
          total_area = polygon_area(polygon.to_a)

          raise ArgumentError, 'At least one site is needed' if sites.empty?
          if targets.size != sites.size
            raise ArgumentError, 'Number of target areas has to match number of sites'
          end
          if targets.any? {|area| area < 0 }
            raise ArgumentError, 'Target areas have to be non-negative'
          end
          raise ArgumentError, 'Polygon area has to be positive' if total_area <= 0
          if (targets.sum - total_area).abs > 1e-6 * total_area
            raise ArgumentError, 'Target areas have to sum up to the polygon area'
          end

          n = sites.size
          weights = Array.new(n, 0.0)

          if unreachable_coincident_sites?(sites, targets, tolerance)
            areas = cells(polygon, sites, weights).map {|cell| polygon_area(cell.to_a) }
            residuals = areas.each_index.map {|i| areas[i] - targets[i] }
            return Result.new(weights: weights, areas: areas, residuals: residuals,
                              iterations: 0, status: :coincident_sites)
          end

          if n == 1
            return Result.new(weights: weights, areas: [total_area],
                              residuals: [total_area - targets[0]], iterations: 0,
                              status: :converged)
          end

          status = :max_iterations
          iterations = 0

          max_iterations.times do |iteration|
            iterations = iteration
            areas = cells(polygon, sites, weights).map {|cell| polygon_area(cell.to_a) }
            residuals = areas.each_index.map {|i| areas[i] - targets[i] }

            empty = sites.each_index.any? do |i|
              areas[i] <= Utils.precision && targets[i] > tolerance
            end
            if empty
              status = :empty_cell
              break
            end
            if residuals.all? {|residual| residual.abs <= tolerance }
              status = :converged
              break
            end

            step = newton_step(polygon, sites, weights, residuals)
            break unless step

            accepted = false
            damping = 1.0
            current_error = residuals.map(&:abs).max
            25.times do
              candidate = weights.each_index.map {|i| weights[i] + damping * step[i] }
              center!(candidate)
              candidate_areas = cells(polygon, sites, candidate).map do |cell|
                polygon_area(cell.to_a)
              end
              error = candidate_areas.each_index.map do |i|
                (candidate_areas[i] - targets[i]).abs
              end.max
              if error < current_error
                weights = candidate
                accepted = true
                break
              end
              damping *= 0.5
            end
            break unless accepted

            iterations = iteration + 1
          end

          areas = cells(polygon, sites, weights).map {|cell| polygon_area(cell.to_a) }
          residuals = areas.each_index.map {|i| areas[i] - targets[i] }
          if status == :max_iterations &&
              residuals.all? {|residual| residual.abs <= tolerance }
            status = :converged
          end

          Result.new(weights: weights, areas: areas, residuals: residuals,
                     iterations: iterations, status: status)
        end

        private

        # Returns the area of the polygon defined by the given vertices using the shoelace
        # formula.
        def polygon_area(vertices)
          return 0.0 if vertices.size < 3
          area = vertices[-1].wedge(vertices[0])
          0.upto(vertices.size - 2) {|i| area += vertices[i].wedge(vertices[i + 1]) }
          area.abs * 0.5
        end

        # Returns +true+ if two or more sites with positive target areas coincide. Their
        # targets cannot be reached because coincident sites always share the same cell.
        def unreachable_coincident_sites?(sites, targets, tolerance)
          positive = sites.each_index.select {|i| targets[i] > tolerance }
          positive.combination(2).any? do |i, j|
            sites[i].distance(sites[j]) < Utils.precision
          end
        end

        # Centers the weights on zero to fix the additive constant degree of freedom.
        def center!(weights)
          mean = weights.sum / weights.size
          weights.map! {|weight| weight - mean }
        end

        # Clips the polygon given by +vertices+ to the half-plane of points whose power
        # distance to site +si+ with weight +wi+ is at most the power distance to site +sj+
        # with weight +wj+ (Sutherland-Hodgman style) and returns the resulting vertices.
        def clip_half_plane(vertices, si, wi, sj, wj)
          return [] if vertices.empty?

          power_difference = lambda do |point|
            (point.x - si.x)**2 + (point.y - si.y)**2 - wi -
              (point.x - sj.x)**2 - (point.y - sj.y)**2 + wj
          end

          output = []
          previous = vertices.last
          previous_value = power_difference.call(previous)
          vertices.each do |current|
            current_value = power_difference.call(current)
            if current_value <= 0
              if previous_value > 0
                output << bisector_intersection(previous, previous_value,
                                                current, current_value)
              end
              output << current
            elsif previous_value <= 0
              output << bisector_intersection(previous, previous_value,
                                              current, current_value)
            end
            previous = current
            previous_value = current_value
          end
          output
        end

        # Returns the intersection point of the segment from +a+ to +b+ with the radical
        # line of two sites, given the signed power distance differences +value_a+ and
        # +value_b+ at the endpoints.
        def bisector_intersection(a, value_a, b, value_b)
          t = value_a / (value_a - value_b)
          Geom2D::Point(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y))
        end

        # Computes a Newton step for the area residuals using the analytic Jacobian of the
        # cell areas with respect to the weights: dA_i/dw_j is minus the length of the
        # shared edge of cells i and j divided by twice the distance of the sites (for
        # i != j), and dA_i/dw_i is the sum over all neighbors. The step sums to zero so
        # that the weights stay centered. Returns +nil+ if the system cannot be solved.
        def newton_step(polygon, sites, weights, residuals)
          n = sites.size
          jacobian = Array.new(n) { Array.new(n, 0.0) }
          scale = [Math.hypot(polygon.bbox.width, polygon.bbox.height), 1.0].max
          edge_tolerance = 1e-8 * scale * scale

          cells(polygon, sites, weights).each_with_index do |cell, i|
            cell.each_segment do |segment|
              sites.each_index do |j|
                next if i == j
                distance = sites[i].distance(sites[j])
                next if distance < Utils.precision
                on_bisector = [segment.start_point, segment.end_point].all? do |point|
                  ((point.x - sites[i].x)**2 + (point.y - sites[i].y)**2 - weights[i] -
                    (point.x - sites[j].x)**2 - (point.y - sites[j].y)**2 +
                    weights[j]).abs < edge_tolerance
                end
                next unless on_bisector
                derivative = segment.length / (2.0 * distance)
                jacobian[i][j] -= derivative
                jacobian[i][i] += derivative
                break
              end
            end
          end

          # Only sites that can still influence the partition take part in the solve;
          # sites with an empty cell and a zero target area get a zero step.
          active = sites.each_index.select do |i|
            jacobian[i].any? {|value| value.abs > 0 } || residuals[i].abs > 0
          end
          return Array.new(n, 0.0) if active.empty?

          # Gauge condition: the step sums to zero. Eliminate the last active unknown and
          # solve the reduced system for the remaining ones.
          reduced = active[0..-2]
          matrix = reduced.map do |i|
            reduced.map {|j| jacobian[i][j] - jacobian[i][active[-1]] }
          end
          rhs = reduced.map {|i| -residuals[i] }
          solution = solve_linear_system(matrix, rhs)
          return nil unless solution

          step = Array.new(n, 0.0)
          reduced.each_with_index {|row, index| step[row] = solution[index] }
          step[active[-1]] = -solution.sum
          step
        end

        # Solves a dense linear system with Gaussian elimination and partial pivoting.
        # Returns +nil+ if the matrix is singular.
        def solve_linear_system(matrix, rhs)
          n = rhs.size
          augmented = matrix.each_with_index.map {|row, i| row.dup << rhs[i] }

          n.times do |column|
            pivot = (column...n).max_by {|row| augmented[row][column].abs }
            return nil if augmented[pivot][column].abs < 1e-14
            augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
            (column + 1...n).each do |row|
              factor = augmented[row][column] / augmented[column][column]
              (column..n).each do |k|
                augmented[row][k] -= factor * augmented[column][k]
              end
            end
          end

          result = Array.new(n, 0.0)
          (n - 1).downto(0) do |row|
            sum = augmented[row][n]
            (row + 1...n).each {|k| sum -= augmented[row][k] * result[k] }
            result[row] = sum / augmented[row][row]
          end
          result
        end

      end

    end

  end
end
