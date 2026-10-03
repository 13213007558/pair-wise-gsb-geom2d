# -*- frozen_string: true -*-
#
#--
# geom2d - 2D Geometric Objects and Algorithms
# Copyright (C) 2018-2023 Thomas Leitner
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
#++

require 'geom2d/algorithms'
require 'geom2d/algorithms/polygon_operation'
require 'geom2d/polygon_set'

module Geom2D
  module Algorithms

    # Computes weights for a clipped power diagram with prescribed cell areas.
    class PowerDiagram

      Error = Class.new(StandardError)
      DuplicateSitesError = Class.new(Error)
      InfeasibleTargetError = Class.new(Error)

      class ResultError < Error
        attr_reader :result

        def initialize(message, result)
          @result = result
          super(message)
        end
      end

      EmptyCellError = Class.new(ResultError)
      IterationExhausted = Class.new(ResultError)

      Result = Struct.new(:weights, :areas, :area_residuals, :area_residual_norm,
                          :total_area_residual, :coverage_residual, :cells,
                          :iterations, :status, keyword_init: true) do
        def converged?
          status == :optimal
        end
      end

      def self.solve(domain, sites, target_areas, **options)
        new(domain, sites, target_areas, **options).solve
      end

      def self.solve!(domain, sites, target_areas, **options)
        result = solve(domain, sites, target_areas, **options)
        case result.status
        when :empty_cell
          raise EmptyCellError.new('A site with positive target has an empty clipped cell', result)
        when :iteration_exhausted
          raise IterationExhausted.new('The area iteration did not converge', result)
        else
          result
        end
      end

      def initialize(domain, sites, target_areas, max_iterations: 100,
                     tolerance: nil, area_tolerance: nil, site_tolerance: 1e-10)
        @domain = domain
        @sites = sites.map {|site| Geom2D::Point(site) }
        @target_areas = target_areas.map(&:to_f)
        @max_iterations = max_iterations
        @tolerance = tolerance
        @area_tolerance = area_tolerance
        @site_tolerance = site_tolerance
      end

      def solve
        validate!

        active = @target_areas.each_index.select {|i| @target_areas[i].positive? }
        zero = @sites.each_index.to_a - active
        weights = Array.new(@sites.size, 0.0)
        warm_start!(weights, active)
        areas = Array.new(@sites.size, 0.0)
        cells = Array.new(@sites.size) { PolygonSet.new }
        jacobian = Array.new(active.size) { Array.new(active.size, 0.0) }
        iterations = 0

        if active.one?
          active_areas, active_cells, jacobian = evaluate(weights, active)
          assign_evaluation(areas, cells, active, active_areas, active_cells)
        elsif active.size > 1
          iterations = iterate(weights, active, areas, cells, jacobian)
        end

        zero.each {|i| weights[i] = weights.values_at(*active).min - 2.0 * span_squared } if active.any?

        build_result(weights, areas, cells, iterations, final_status(areas))
      end

      private

      def iterate(weights, active, areas, cells, jacobian)
        active_areas, active_cells, jacobian = evaluate(weights, active)
        iteration = 0

        while iteration < @max_iterations
          assign_evaluation(areas, cells, active, active_areas, active_cells)
          residual = active.map {|i| areas[i] - @target_areas[i] }
          return iteration if converged?(residual)

          direction = newton_direction(residual, jacobian)
          accepted = step(weights, active, direction, active_areas, active_cells, jacobian)

          unless accepted
            residual = active.map {|i| active_areas[i] - @target_areas[i] }
            gradient = active.map.with_index do |_site_index, k|
              active.map.with_index {|_other, j| residual[j] * jacobian[j][k] }.sum
            end
            mean = gradient.sum / gradient.size
            direction = gradient.map {|value| -(value - mean) / [span_squared, 1.0].max }
            accepted = step(weights, active, direction, active_areas, active_cells, jacobian)
          end

          break unless accepted

          active_areas, active_cells, jacobian = evaluate(weights, active)
          iteration += 1
        end

        assign_evaluation(areas, cells, active, active_areas, active_cells)
        iteration
      end

      def step(weights, active, direction, old_areas, old_cells, old_jacobian)
        old_loss = loss(active.map {|i| old_areas[i] - @target_areas[i] })
        scale = 1.0
        old_weights = active.map {|i| weights[i] }

        12.times do
          active.each_with_index {|i, k| weights[i] = old_weights[k] + scale * direction[k] }
          center_active_weights!(weights, active)
          candidate = candidate_evaluation(weights, active)
          unless candidate
            active.each_with_index {|i, k| weights[i] = old_weights[k] }
            scale *= 0.5
            next
          end

          new_areas, new_cells, new_jacobian = candidate
          residual = active.map.with_index {|_, k| new_areas[k] - @target_areas[active[k]] }

          if loss(residual) <= old_loss + 1e-14 * @domain_area**2
            old_areas.replace(new_areas)
            old_cells.replace(new_cells)
            old_jacobian.replace(new_jacobian)
            return true
          end

          scale *= 0.5
        end

        active.each_with_index {|i, k| weights[i] = old_weights[k] }
        false
      end

      def newton_direction(residual, jacobian)
        size = residual.size
        reduced = Array.new(size - 1) do |i|
          Array.new(size - 1) do |j|
            jacobian[i][j] - jacobian[i][-1] - jacobian[-1][j] + jacobian[-1][-1]
          end
        end
        rhs = Array.new(size - 1) {|i| -(residual[i] - residual[-1]) }
        delta = solve_linear(reduced, rhs)
        delta + [-delta.sum]
      end

      def candidate_evaluation(weights, active)
        evaluate(weights, active)
      rescue StandardError
        nil
      end

      def solve_linear(matrix, rhs)
        ridge = 1e-12 * span_squared
        8.times do
          begin
            return gaussian_solve(matrix, rhs, ridge)
          rescue ArgumentError
            ridge = [ridge * 10.0, 1e-8 * span_squared].max
          end
        end
        Array.new(rhs.size, 0.0)
      end

      def gaussian_solve(matrix, rhs, ridge)
        size = rhs.size
        a = matrix.map.with_index {|row, i| row.map.with_index {|value, j| value + (i == j ? ridge : 0.0) } }
        b = rhs.dup

        size.times do |pivot|
          row = (pivot...size).max_by {|i| a[i][pivot].abs }
          if a[row][pivot].abs <= 1e-18 * [span_squared, 1.0].max
            raise ArgumentError, 'Singular power diagram Jacobian'
          end

          a[pivot], a[row] = a[row], a[pivot]
          b[pivot], b[row] = b[row], b[pivot]

          (pivot + 1...size).each do |row_index|
            factor = a[row_index][pivot] / a[pivot][pivot]
            (pivot...size).each {|column| a[row_index][column] -= factor * a[pivot][column] }
            b[row_index] -= factor * b[pivot]
          end
        end

        solution = Array.new(size, 0.0)
        (size - 1).downto(0) do |i|
          sum = (i + 1...size).sum {|j| a[i][j] * solution[j] }
          solution[i] = (b[i] - sum) / a[i][i]
        end
        solution
      end

      def validate!
        raise ArgumentError, 'At least one site is required' if @sites.empty?
        unless @sites.size == @target_areas.size
          raise ArgumentError, 'Sites and target areas must have the same size'
        end
        raise ArgumentError, 'Target areas must be finite' unless @target_areas.all?(&:finite?)
        raise InfeasibleTargetError, 'Target areas must be non-negative' if @target_areas.any?(&:negative?)
        validate_sites!

        @domain_area = polygon_set_area(@domain)
        raise ArgumentError, 'The bounded polygon must have positive area' unless @domain_area.positive?

        tolerance = @area_tolerance || [@domain_area * 1e-9, 1e-9].max
        total = @target_areas.sum
        return if (total - @domain_area).abs <= tolerance

        raise InfeasibleTargetError,
              "Target area sum #{total} does not match domain area #{@domain_area}"
      end

      def validate_sites!
        positive_sites = @target_areas.each_index.select {|i| @target_areas[i].positive? }
        positive_sites.combination(2).each do |i, j|
          next unless @sites[i].distance(@sites[j]) <= @site_tolerance

          raise DuplicateSitesError,
                "Sites at #{i} and #{j} coincide but both have positive target area"
        end
      end


      def warm_start!(weights, active)
        return if active.size < 2

        areas, _cells, jacobian = evaluate(weights, active)
        if active.each_with_index.any? {|_, k| areas[k].zero? }
          expand_active_cells!(weights, active, areas)
          areas, _cells, jacobian = evaluate(weights, active)
          active.size.times do
            break unless active.each_with_index.any? {|_, k| areas[k].zero? }

            expand_active_cells!(weights, active, areas)
            areas, _cells, jacobian = evaluate(weights, active)
          end
        end
        30.times do
          residual = active.map.with_index {|_, k| areas[k] - @target_areas[active[k]] }
          return if residual.all? {|value| value.abs <= convergence_tolerance }

          direction = newton_direction(residual, jacobian)
          direction.map! {|value| [[value, -span_squared].max, span_squared].min }
          old_weights = active.map {|i| weights[i] }
          old_loss = loss(residual)
          scale = 0.25

          accepted = 8.times.any? do
            active.each_with_index {|i, k| weights[i] = old_weights[k] + scale * direction[k] }
            center_active_weights!(weights, active)
            candidate = candidate_evaluation(weights, active)
            unless candidate
              active.each_with_index {|i, k| weights[i] = old_weights[k] }
              scale *= 0.5
              next
            end

            new_areas, _new_cells, new_jacobian = candidate
            new_residual = active.map.with_index {|_, k| new_areas[k] - @target_areas[active[k]] }
            if new_areas.all?(&:positive?) && loss(new_residual) < old_loss
              areas = new_areas
              jacobian = new_jacobian
              true
            else
              scale *= 0.5
              false
            end
          end

          unless accepted
            active.each_with_index {|i, k| weights[i] = old_weights[k] }
            return
          end
        end
      end

      def convergence_tolerance
        @tolerance || [@domain_area * 5e-8, 5e-8].max
      end

      def expand_active_cells!(weights, active, areas)
        empty = active.each_with_index.select {|_i, k| areas[k].zero? }.map(&:first)
        return if empty.empty?

        points = interior_sample_points.shuffle(random: Random.new(0))
        empty.each do |site_index|
          point = points.max_by {|sample| -distance_squared(sample, @sites[site_index]) }
          break unless point

          points.delete(point)
          minimum = active.map {|other_index| power_value(point, other_index, weights) }.min
          weights[site_index] = distance_squared(point, @sites[site_index]) - minimum
        end
        empty.each {|site_index| weights[site_index] += 1e-6 * span_squared }
        center_active_weights!(weights, active)
      end

      def interior_sample_points
        return @interior_sample_points if defined?(@interior_sample_points)

        box = bounding_box
        @interior_sample_points = nil
        2.upto(40) do |resolution|
          points = []
          resolution.times do |ix|
            resolution.times do |iy|
              point = Geom2D::Point(box.min_x + (ix + 0.5) * box.width / resolution,
                                   box.min_y + (iy + 0.5) * box.height / resolution)
              points << point if point_in_polygon_set?(point, @domain)
            end
          end
          if points.size >= @sites.size
            @interior_sample_points = points
            break
          end
        end
        @interior_sample_points ||= []
      end

      def center_active_weights!(weights, active)
        offset = weights.values_at(*active).sum / active.size
        active.each {|i| weights[i] -= offset }
      end

      def evaluate(weights, active)
        polygons = active.map {|i| power_polygon(i, weights, active) }
        areas = Array.new(active.size, 0.0)
        cells = Array.new(active.size) { PolygonSet.new }

        polygons.each_with_index do |polygon, k|
          next if polygon.nil?

          cells[k] = intersect_domain(polygon)
          areas[k] = polygon_set_area(cells[k])
        end

        jacobian = Array.new(active.size) { Array.new(active.size, 0.0) }
        build_jacobian!(jacobian, polygons, active, weights)
        [areas, cells, jacobian]
      end

      def power_polygon(site_index, weights, active)
        box = bounding_box
        weight_scale = Math.sqrt(active.map {|i| weights[i].abs }.max + 1.0)
        margin = [box.width, box.height, 1.0].max * 100.0 + weight_scale * 1.0e3
        vertices = [[box.min_x - margin, box.min_y - margin, nil],
                    [box.max_x + margin, box.min_y - margin, nil],
                    [box.max_x + margin, box.max_y + margin, nil],
                    [box.min_x - margin, box.max_y + margin, nil]]

        active.each do |other_index|
          next if other_index == site_index
          vertices = clip_half_plane(vertices, site_index, other_index, weights)
          return nil if vertices.empty?
        end
        vertices = clean_vertices(vertices)
        return nil if vertices.size < 3

        Polygon.new(vertices.map {|x, y, _| [x, y] })
      end

      def clean_vertices(vertices)
        result = []
        vertices.each do |vertex|
          if result.empty? || !same_vertex?(result.last, vertex)
            result << vertex
          else
            result[-1] = [result.last[0], result.last[1], result.last[2] || vertex[2]]
          end
        end
        result.pop if result.size > 1 && same_vertex?(result.first, result.last)
        result.size.times { result = remove_collinear_vertex(result) }
        result
      end

      def same_vertex?(a, b)
        (a[0] - b[0]).abs <= 1e-10 && (a[1] - b[1]).abs <= 1e-10
      end

      def remove_collinear_vertex(vertices)
        vertices.each_with_index do |vertex, index|
          previous = vertices[(index - 1) % vertices.size]
          next_vertex = vertices[(index + 1) % vertices.size]
          cross = (vertex[0] - previous[0]) * (next_vertex[1] - previous[1]) -
                  (vertex[1] - previous[1]) * (next_vertex[0] - previous[0])
          next unless cross.abs <= 1e-9 * [span_squared, 1.0].max

          return vertices[0...index] + vertices[(index + 1)..]
        end
        vertices
      end

      def clip_half_plane(vertices, site_index, other_index, weights)
        result = []
        site = @sites[site_index]
        other = @sites[other_index]
        a = 2.0 * (other.x - site.x)
        b = 2.0 * (other.y - site.y)
        c = other.x**2 + other.y**2 - site.x**2 - site.y**2 -
            weights[other_index] + weights[site_index]

        previous = vertices.last
        vertices.each do |current|
          current_value = a * current[0] + b * current[1] - c
          previous_value = a * previous[0] + b * previous[1] - c

          if current_value <= @site_tolerance
            result << intersect_clip_edge(previous, current, previous_value,
                                          current_value, a, b, c, other_index)
            result << current
          elsif previous_value <= @site_tolerance
            result << intersect_clip_edge(previous, current, previous_value,
                                          current_value, a, b, c, other_index)
          end
          previous = current
        end
        result.compact
      end

      def intersect_clip_edge(previous, current, previous_value, current_value, a, b, c, label)
        return previous if previous_value == current_value

        t = previous_value / (previous_value - current_value)
        x = previous[0] + (current[0] - previous[0]) * t
        y = previous[1] + (current[1] - previous[1]) * t
        boundary = (a * x + b * y - c).abs <= 1e-9 * [span_squared, 1.0].max
        [x, y, boundary ? label : nil]
      end

      def build_jacobian!(jacobian, polygons, active, weights)
        polygons.each_with_index do |polygon, k|
          next if polygon.nil?

          polygon.each_segment do |segment|
            label = segment_label(polygon, segment, k, active, weights)
            next unless label

            length = clipped_segment_length(segment)
            next unless length.positive?

            j = active.index(label)
            jacobian[k][j] -= length / 2.0
            jacobian[k][k] += length / 2.0
          end
        end
      end

      def segment_label(polygon, segment, k, active, weights)
        midpoint = Geom2D::Point((segment.start_point.x + segment.end_point.x) / 2.0,
                                 (segment.start_point.y + segment.end_point.y) / 2.0)
        tolerance = 1e-8 * [span_squared, 1.0].max
        site_index = active[k]
        value = power_value(midpoint, site_index, weights)
        active.each do |other_index|
          next if other_index == site_index
          if (power_value(midpoint, other_index, weights) - value).abs <= tolerance
            return other_index
          end
        end
        nil
      end

      def clipped_segment_length(segment)
        if as_polygon_set(@domain).polygons.size == 1
          domain = as_polygon_set(@domain).polygons.first
          return convex_domain_segment_length(segment, domain) if convex_polygon?(domain)
        end

        intersections = [0.0, 1.0]
        @domain.each_segment do |boundary|
          intersection = segment.intersect(boundary)
          add_segment_root!(intersections, segment, intersection)
        end
        intersections.sort!
        intersections.uniq {|value| (value * 1e10).round }

        length = 0.0
        intersections.each_cons(2) do |start_t, end_t|
          next unless end_t - start_t > 1e-10

          midpoint = point_on_segment(segment, (start_t + end_t) / 2.0)
          length += segment.length * (end_t - start_t) if point_in_polygon_set?(midpoint, @domain)
        end
        length
      end

      def convex_domain_segment_length(segment, domain)
        a = segment.start_point
        b = segment.end_point
        direction_x = b.x - a.x
        direction_y = b.y - a.y
        lower = 0.0
        upper = 1.0
        ccw = domain.ccw?

        domain.each_segment do |boundary|
          edge_x = boundary.end_point.x - boundary.start_point.x
          edge_y = boundary.end_point.y - boundary.start_point.y
          value = edge_x * (a.y - boundary.start_point.y) -
                  edge_y * (a.x - boundary.start_point.x)
          derivative = edge_x * direction_y - edge_y * direction_x
          value = -value unless ccw
          derivative = -derivative unless ccw

          if derivative.abs < 1e-15
            return 0.0 if value < 0
          else
            parameter = -value / derivative
            if derivative.positive?
              lower = [lower, parameter].max
            else
              upper = [upper, parameter].min
            end
          end
        end

        return 0.0 if upper <= lower

        segment.length * (upper - lower)
      end

      def add_segment_root!(roots, segment, intersection)
        case intersection
        when NilClass
          nil
        when Geom2D::Point
          root = segment_parameter(segment, intersection)
          roots << root if root
        when Geom2D::Segment
          root = segment_parameter(segment, intersection.start_point)
          roots << root if root
          root = segment_parameter(segment, intersection.end_point)
          roots << root if root
        end
      end

      def segment_parameter(segment, point)
        dx = segment.end_point.x - segment.start_point.x
        dy = segment.end_point.y - segment.start_point.y
        t = if dx.abs > dy.abs
              (point.x - segment.start_point.x) / dx
            else
              (point.y - segment.start_point.y) / dy
            end
        return nil unless t.finite? && t >= -1e-9 && t <= 1.0 + 1e-9

        [[t, 0.0].max, 1.0].min
      end

      def point_on_segment(segment, parameter)
        Geom2D::Point(segment.start_point.x + (segment.end_point.x - segment.start_point.x) * parameter,
                      segment.start_point.y + (segment.end_point.y - segment.start_point.y) * parameter)
      end

      def point_in_polygon_set?(point, polygon_set)
        as_polygon_set(polygon_set).polygons.any? {|polygon| point_in_polygon?(point, polygon) }
      end

      def point_in_polygon?(point, polygon)
        return false if polygon.each_segment.any? {|segment| point_on_segment_boundary?(point, segment) }

        inside = false
        polygon.each_segment do |segment|
          a = segment.start_point
          b = segment.end_point
          next unless (a.y > point.y) != (b.y > point.y)

          cross_x = a.x + (point.y - a.y) * (b.x - a.x) / (b.y - a.y)
          inside = !inside if point.x < cross_x
        end
        inside
      end

      def point_on_segment_boundary?(point, segment)
        cross = ((point.x - segment.start_point.x) * (segment.end_point.y - segment.start_point.y) -
                 (point.y - segment.start_point.y) * (segment.end_point.x - segment.start_point.x)).abs
        return false if cross > 1e-9 * [span_squared, 1.0].max

        dot = (point - segment.start_point) * (point - segment.end_point)
        dot <= 1e-9 * [span_squared, 1.0].max
      end

      def intersect_domain(polygon)
        domain_set = as_polygon_set(@domain)
        if domain_set.polygons.size == 1 && convex_polygon?(domain_set.polygons.first)
          clipped = clip_to_convex_domain(polygon, domain_set.polygons.first)
          return clipped if clipped
        end

        Algorithms::PolygonOperation.run(PolygonSet.new([polygon]), domain_set, :intersection)
      end

      def convex_polygon?(polygon)
        vertices = polygon.to_a
        return false if vertices.size < 3

        signs = vertices.each_with_index.map do |point, i|
          other = vertices[(i + 1) % vertices.size]
          next_vertex = vertices[(i + 2) % vertices.size]
          Algorithms.ccw(point, other, next_vertex)
        end
        signs.all?(&:positive?) || signs.all?(&:negative?)
      end

      def clip_to_convex_domain(polygon, domain)
        return nil if domain.each_segment.any?(&:degenerate?)

        ccw = domain.ccw?
        vertices = polygon.to_a
        domain.each_segment do |segment|
          a = segment.start_point
          b = segment.end_point
          vertices = clip_to_domain_edge(vertices, a, b, ccw)
          return PolygonSet.new if vertices.size < 3
        end
        PolygonSet.new([Polygon.new(vertices)])
      rescue StandardError
        nil
      end

      def clip_to_domain_edge(vertices, a, b, ccw)
        result = []
        previous = vertices.last
        previous_value = edge_side(previous, a, b)
        previous_inside = ccw ? previous_value >= 0 : previous_value <= 0

        vertices.each do |current|
          current_value = edge_side(current, a, b)
          current_inside = ccw ? current_value >= 0 : current_value <= 0
          if current_inside != previous_inside
            denominator = previous_value - current_value
            t = previous_value / denominator
            result << Geom2D::Point(previous.x + (current.x - previous.x) * t,
                                   previous.y + (current.y - previous.y) * t)
          end
          result << current if current_inside
          previous = current
          previous_value = current_value
          previous_inside = current_inside
        end
        result
      end

      def edge_side(point, a, b)
        (b.x - a.x) * (point.y - a.y) - (b.y - a.y) * (point.x - a.x)
      end

      def as_polygon_set(value)
        value.is_a?(PolygonSet) ? value : PolygonSet.new([value])
      end

      def polygon_set_area(value)
        as_polygon_set(value).polygons.sum {|polygon| polygon_area(polygon).abs }
      end

      def polygon_area(polygon)
        vertices = polygon.to_a
        return 0.0 if vertices.size < 3

        vertices.each_with_index.sum do |point, i|
          other = vertices[(i + 1) % vertices.size]
          (other.x - point.x) * (other.y + point.y)
        end.abs / 2.0
      end

      def bounding_box
        @bounding_box ||= begin
          box = as_polygon_set(@domain).bbox
          @sites.each {|site| box.add!(site) }
          box
        end
      end

      def span_squared
        box = bounding_box
        @span_squared ||= box.width**2 + box.height**2
      end

      def power_value(point, site_index, weights)
        distance_squared(point, @sites[site_index]) - weights[site_index]
      end

      def distance_squared(point, site)
        (point.x - site.x)**2 + (point.y - site.y)**2
      end

      def loss(residual)
        residual.sum {|value| value * value }
      end

      def converged?(residual)
        tolerance = convergence_tolerance
        residual.all? {|value| value.abs <= tolerance }
      end

      def assign_evaluation(areas, cells, active, active_areas, active_cells)
        active.each_with_index do |site_index, k|
          areas[site_index] = active_areas[k]
          cells[site_index] = active_cells[k]
        end
      end

      def final_status(areas)
        active = @target_areas.each_index.select {|i| @target_areas[i].positive? }
        positive_empty = active.any? {|i| areas[i] <= 1e-12 * [@domain_area, 1.0].max }
        return :empty_cell if positive_empty

        residual = active.map {|i| areas[i] - @target_areas[i] }
        converged?(residual) ? :optimal : :iteration_exhausted
      end

      def build_result(weights, areas, cells, iterations, status)
        center_weights!(weights)
        residuals = @target_areas.each_index.map {|i| areas[i] - @target_areas[i] }
        coverage = areas.sum - @domain_area
        Result.new(
          weights: weights,
          areas: areas,
          area_residuals: residuals,
          area_residual_norm: Math.sqrt(residuals.sum {|value| value * value }),
          total_area_residual: coverage,
          coverage_residual: coverage.abs,
          cells: cells,
          iterations: iterations,
          status: status
        )
      end

      def center_weights!(weights)
        offset = weights.sum / weights.size
        weights.map! {|weight| weight - offset }
      end

    end
  end
end
