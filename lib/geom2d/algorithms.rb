# -*- frozen_string_literal: true -*-
#
#--
# geom2d - 2D Geometric Objects and Algorithms
# Copyright (C) 2018-2023 Thomas Leitner <t_leitner@gmx.at>
#
# This software may be modified and distributed under the terms
# of the MIT license.  See the LICENSE file for details.
#++

require 'geom2d/utils'
require 'geom2d/point'

module Geom2D

  # This module contains helper functions as well as classes implementing algorithms.
  module Algorithms

    autoload(:PolygonOperation, 'geom2d/algorithms/polygon_operation')
    autoload(:MovingConvexHull, 'geom2d/algorithms/moving_convex_hull')

    extend Utils

    # Determines whether the three points form a counterclockwise turn.
    #
    # Returns
    #
    # * +1 if the points a -> b -> c form a counterclockwise angle,
    # * -1 if the points a -> b -> c from a clockwise angle, and
    # *  0 if the points are collinear.
    def self.ccw(a, b, c)
      float_compare((b.x - a.x) * (c.y - a.y), (c.x - a.x) * (b.y - a.y))
    end

    # Returns the strict, ordered convex hull of the given points.
    #
    # The vertices are returned in counter-clockwise order. Points that coincide are represented
    # only once, and points in the interior of a hull edge, including collinear points at an event
    # time, are omitted.
    #
    # The optional +epsilon+ defines the absolute numerical tolerance used for cross products and
    # coincident points; it does not set a minimum distance between events. It defaults to
    # Utils.precision.
    def self.convex_hull(points, epsilon: Utils.precision)
      entries = points.each_with_index.map {|point, index| [point.x, point.y, index] }
      indexed_convex_hull(entries, epsilon).map {|x, y, _| Point.new(x, y) }
    end

    # Analyzes the convex hull of points moving with constant velocity on the finite closed
    # interval [+t0+, +t1+].
    #
    # Each moving point may be a MovingConvexHull::MovingPoint or an object responding to
    # +initial_point+ and +velocity+. See MovingConvexHull.analyze for the returned event
    # structure and all options.
    def self.moving_convex_hull(moving_points, t0 = 0.0, t1 = 1.0, **options)
      MovingConvexHull.analyze(moving_points, t0, t1, **options)
    end

    # :nodoc:
    def self.indexed_convex_hull(entries, epsilon)
      return [] if entries.empty?

      ordered = entries.sort_by {|x, y, id| [x, y, id] }
      unique = []
      ordered.each do |entry|
        x, y = entry
        next if unique.any? {|other| (x - other[0]).abs < epsilon && (y - other[1]).abs < epsilon }

        unique << entry
      end
      return unique.dup if unique.length == 1

      lower = []
      unique.each do |entry|
        while lower.length >= 2
          ax, ay = lower[-2]
          bx, by = lower[-1]
          cx, cy = entry
          cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
          break unless cross <= epsilon

          lower.pop
        end
        lower << entry
      end

      upper = []
      unique.reverse_each do |entry|
        while upper.length >= 2
          ax, ay = upper[-2]
          bx, by = upper[-1]
          cx, cy = entry
          cross = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
          break unless cross <= epsilon

          upper.pop
        end
        upper << entry
      end

      lower[0...-1] + upper[0...-1]
    end

    private_class_method :indexed_convex_hull

  end

end
