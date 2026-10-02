# -*- frozen_string_literal: true -*-
#
#--
# geom2d - 2D Geometric Objects and Algorithms
# Copyright (C) 2018-2023 Thomas Leitner
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
#++

require 'bigdecimal'
require 'bigdecimal/util'
require 'bigdecimal/math'
require 'geom2d/point'

module Geom2D
  module Algorithms

    # Performs an event-driven analysis of the convex hull of linearly moving points.
    class MovingConvexHull

      NUMERIC_PRECISION = 200
      private_constant(:NUMERIC_PRECISION)

      # A point at +initial_point+ moving with constant +velocity+.
      class MovingPoint

        # The position at the beginning of the analyzed interval.
        attr_reader :initial_point

        # The constant velocity vector.
        attr_reader :velocity

        # An optional identifier used in the analysis result.
        attr_reader :id

        # Creates a moving point.
        def initialize(initial_point, velocity = [0, 0], id: nil)
          @initial_point = MovingConvexHull.coerce_point(initial_point)
          @velocity = MovingConvexHull.coerce_point(velocity)
          @id = id
        end

        # Returns the point's position at time +t+.
        def position(t)
          Point.new(initial_point.x + velocity.x * t,
                    initial_point.y + velocity.y * t)
        end

      end

      # Analyzes +moving_points+ on the finite closed interval [+t0+, +t1+].
      #
      # Candidate events are pair-coincidence times and triple-collinearity times strictly inside
      # the interval. Simultaneous candidates are evaluated as one event, and adjacent open
      # intervals are never removed merely because they are short. Event records are retained when
      # the cyclic hull order changes or when distinct points coincide on a hull vertex at the
      # exact event time.
      #
      # Points whose initial positions and velocities agree within +epsilon+ are treated as one
      # trajectory for the complete interval; the smallest index is the representative and all IDs
      # in the group remain available in the result.
      #
      # +epsilon+ is an absolute geometric tolerance for cross products, coincident coordinates,
      # and always-identical trajectories. Candidate roots are solved with high-precision exact
      # sign tests, so +epsilon+ does not suppress short intervals. +time_tolerance+ defaults to
      # zero and only clusters adjacent numerical event times when explicitly provided.
      #
      # The returned Hash contains +point_ids+, +coincident_groups+, retained +events+, and
      # +intervals+. Event pair and triple indices identify representative trajectories. Each
      # interval contains ordered representative indices and IDs, the full ID groups at those
      # vertices, moving points, and its two open endpoints.
      def self.analyze(moving_points, t0, t1,
                       epsilon: Utils.precision, time_tolerance: 0.0)
        precision = epsilon.to_d
        start_time = t0.to_d
        end_time = t1.to_d
        raise ArgumentError, 't1 must be greater than t0' unless end_time > start_time

        points = moving_points.map.with_index {|point, index| normalize(point, index) }
        validate_ids!(points)

        representatives, groups = coincident_representatives(points, precision)
        event_times = candidate_event_times(points, representatives, start_time, end_time)
        event_times = group_times(event_times, time_tolerance.to_d)
        changes = hull_changes(points, representatives, event_times, start_time, end_time,
                               precision)
        intervals = build_intervals(points, groups, representatives, changes,
                                    start_time, end_time, precision)

        {
          t0: t0,
          t1: t1,
          point_ids: points.map(&:id),
          coincident_groups: groups.select {|group| group.length > 1 },
          events: changes.map {|change| event_result(change, intervals) },
          intervals: intervals
        }
      end

      # Returns the ordered representative indices at +t+.
      def self.hull_vertex_indices_at(moving_points, t, epsilon: Utils.precision)
        points = moving_points.map.with_index {|point, index| normalize(point, index) }
        representatives, = coincident_representatives(points, epsilon.to_d)
        hull_at(points, representatives, t.to_d, epsilon.to_d)
      end

      def self.coerce_point(value)
        return value if value.is_a?(Point)
        return Point.new(value.x, value.y) if value.respond_to?(:x) && value.respond_to?(:y)

        Point.new(*value)
      end

      def self.normalize(point, index)
        case point
        when MovingPoint
          point.id.nil? ? MovingPoint.new(point.initial_point, point.velocity, id: index) : point
        when Hash
          initial = point.fetch(:initial_point, point[:initial])
          MovingPoint.new(initial, point.fetch(:velocity, [0, 0]),
                          id: point.fetch(:id, index))
        else
          MovingPoint.new(point.initial_point, point.velocity,
                          id: point.respond_to?(:id) ? point.id : index)
        end
      end

      private_class_method :normalize

      def self.validate_ids!(points)
        ids = points.map(&:id).reject(&:nil?)
        raise ArgumentError, 'Moving point IDs must be unique' unless ids.uniq.length == ids.length
      end

      private_class_method :validate_ids!

      def self.coincident_representatives(points, epsilon)
        groups = []
        representatives = []

        points.each_with_index do |point, index|
          group = groups.find {|existing| same_trajectory?(point, points[existing[0]], epsilon) }
          if group
            group << index
          else
            groups << [index]
            representatives << index
          end
        end

        [representatives, groups]
      end

      private_class_method :coincident_representatives

      def self.same_trajectory?(point, other, epsilon)
        [:x, :y].all? do |coordinate|
          (point.initial_point.public_send(coordinate).to_d -
            other.initial_point.public_send(coordinate).to_d).abs < epsilon &&
            (point.velocity.public_send(coordinate).to_d -
              other.velocity.public_send(coordinate).to_d).abs < epsilon
        end
      end

      private_class_method :same_trajectory?

      def self.candidate_event_times(points, indices, start_time, end_time)
        times = []

        indices.combination(3) do |i, j, k|
          roots(*triple_coefficients(points, i, j, k)).each do |time|
            times << time if time > start_time && time < end_time
          end
        end

        indices.combination(2) do |i, j|
          time = pair_meeting_time(points, i, j)
          times << time if time && time > start_time && time < end_time
        end

        times.sort.uniq
      end

      private_class_method :candidate_event_times

      def self.triple_coefficients(points, i, j, k)
        p0 = point_value(points, i, :initial_point)
        p1 = point_value(points, j, :initial_point)
        p2 = point_value(points, k, :initial_point)
        v0 = point_value(points, i, :velocity)
        v1 = point_value(points, j, :velocity)
        v2 = point_value(points, k, :velocity)
        initial01 = p1 - p0
        velocity01 = v1 - v0
        initial02 = p2 - p0
        velocity02 = v2 - v0

        [
          velocity01.wedge(velocity02),
          velocity01.wedge(initial02) + initial01.wedge(velocity02),
          initial01.wedge(initial02)
        ]
      end

      private_class_method :triple_coefficients

      def self.point_value(points, index, attribute)
        point = points[index].public_send(attribute)
        Point.new(point.x.to_d, point.y.to_d)
      end

      private_class_method :point_value

      def self.pair_meeting_time(points, i, j)
        dx = points[j].initial_point.x.to_d - points[i].initial_point.x.to_d
        dy = points[j].initial_point.y.to_d - points[i].initial_point.y.to_d
        vx = points[j].velocity.x.to_d - points[i].velocity.x.to_d
        vy = points[j].velocity.y.to_d - points[i].velocity.y.to_d

        if vx.abs >= vy.abs
          return nil if vx.zero?

          time = -dx / vx
        else
          return nil if vy.zero?

          time = -dy / vy
        end

        (dx + vx * time).zero? && (dy + vy * time).zero? ? time : nil
      end

      private_class_method :pair_meeting_time

      def self.roots(a, b, c)
        if a.zero?
          return [] if b.zero?

          return [-c / b]
        end

        discriminant = b * b - 4 * a * c
        return [] if discriminant.negative?
        return [-b / (2 * a)] if discriminant.zero?

        root = BigMath.sqrt(discriminant, NUMERIC_PRECISION)
        [(-b - root) / (2 * a), (-b + root) / (2 * a)].uniq.sort
      end

      private_class_method :roots

      def self.group_times(times, tolerance)
        return times if tolerance.zero?

        groups = []
        times.each do |time|
          if groups.empty? || (time - groups.last).abs > tolerance
            groups << time
          end
        end
        groups
      end

      private_class_method :group_times

      def self.hull_changes(points, indices, event_times, start_time, end_time, epsilon)
        boundaries = [start_time]

        event_times.filter_map do |time|
          index = event_times.index(time)
          right_time = index + 1 < event_times.length ? event_times[index + 1] : end_time
          previous = hull_at(points, indices, midpoint(boundaries.last, time), epsilon)
          following = hull_at(points, indices, midpoint(time, right_time), epsilon)
          event_hull = hull_at(points, indices, time, epsilon)
          hull_pairs = indices.combination(2).select do |i, j|
            left = points[i].position(time)
            right = points[j].position(time)
            (right.x - left.x).abs < epsilon &&
              (right.y - left.y).abs < epsilon &&
              ([i, j] & event_hull).any?
          end
          boundaries << time
          next if previous == following && hull_pairs.empty?

          {
            time: time,
            triples: triples_at(points, indices, time, epsilon),
            pairs: hull_pairs,
            previous_hull: previous,
            next_hull: following
          }
        end
      end

      private_class_method :hull_changes

      def self.midpoint(left, right)
        left + (right - left) * 0.5.to_d
      end

      private_class_method :midpoint

      def self.triples_at(points, indices, time, epsilon)
        indices.combination(3).select do |i, j, k|
          a, b, c = triple_coefficients(points, i, j, k)
          (a * time * time + b * time + c).abs < epsilon
        end
      end

      private_class_method :triples_at

      def self.build_intervals(points, groups, indices, changes, start_time, end_time, epsilon)
        boundaries = [start_time, *changes.map {|change| change[:time] }, end_time]
        boundaries.each_cons(2).with_index.map do |(left, right), index|
          hull = if changes.empty?
                   hull_at(points, indices, midpoint(left, right), epsilon)
                 elsif index.zero?
                   changes.first[:previous_hull]
                 else
                   changes[index - 1][:next_hull]
                 end
          interval_result(left, right, hull, points, groups)
        end
      end

      private_class_method :build_intervals

      def self.hull_at(points, indices, time, epsilon)
        entries = indices.map do |index|
          point = points[index].position(time)
          [point.x, point.y, index]
        end
        hull = Algorithms.send(:indexed_convex_hull, entries, epsilon)
        canonicalize(hull.map {|_, _, index| index })
      end

      private_class_method :hull_at

      def self.canonicalize(indices)
        return indices.sort if indices.length <= 2

        offset = indices.each_with_index.min[1]
        indices.rotate(offset)
      end

      private_class_method :canonicalize

      def self.interval_result(left, right, hull, points, groups)
        group_by_vertex = groups.each_with_object({}) do |group, result|
          result[group[0]] = group.map {|index| points[index].id }
        end

        {
          t0: left,
          t1: right,
          open: true,
          hull_vertex_indices: hull,
          hull_vertex_ids: hull.map {|index| points[index].id },
          hull_vertex_groups: hull.map {|index| group_by_vertex[index] },
          hull_moving_points: hull.map {|index| points[index] }
        }
      end

      private_class_method :interval_result

      def self.event_result(change, intervals)
        {
          time: change[:time],
          previous_interval_index: intervals.index {|interval| interval[:t1] == change[:time] },
          next_interval_index: intervals.index {|interval| interval[:t0] == change[:time] },
          triples: change[:triples],
          pairs: change[:pairs]
        }
      end

      private_class_method :event_result

    end

  end
end
