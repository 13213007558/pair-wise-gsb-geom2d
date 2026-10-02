# -*- frozen_string_literal: true -*-

require 'geom2d/algorithms'
require 'geom2d/algorithms/convex_hull'
require 'geom2d/point'
require 'geom2d/utils'

module Geom2D
  module Algorithms

    # Analyzes the convex hull of points moving with constant velocities over a finite time
    # interval (a "kinetic" or "moving" convex hull).
    #
    # Each point i moves as p_i(t) = position_i + t * velocity_i. Instead of recomputing the
    # hull at fixed time steps, this class determines the exact event times at which the hull
    # changes its vertex set, together with the ordered hull vertices for each open interval
    # between two adjacent events.
    #
    # == How events are found
    #
    # The convex hull (as an ordered vertex list) can only change when three points become
    # collinear. For linearly moving points the orientation of any triple is a quadratic
    # polynomial in t, so all candidate event times are found by solving one quadratic per
    # triple - direct enumeration that is exact and therefore cannot miss vertices that stay
    # on the hull only for a very short time. Two points coinciding is covered as well: if
    # p_i(t) == p_j(t), the triple (i, j, k) is degenerate (collinear) for every k.
    #
    # The hull for each open interval between adjacent candidate events is computed at the
    # interval midpoint with ConvexHull. Adjacent intervals with identical vertex lists are
    # merged, so only events at which the hull actually changes members are reported (e.g. a
    # collinearity of interior points or a tangential touch produces no reported event).
    #
    # == Numeric tolerances
    #
    # All tolerances derive from Utils.precision (default 1e-10):
    #
    # * Quadratic and linear coefficients of the orientation polynomial with an absolute value
    #   below Utils.precision are treated as zero. A triple whose orientation polynomial is
    #   identically zero (permanently collinear) generates no events.
    # * A negative discriminant with absolute value below Utils.precision is treated as zero,
    #   i.e. as a tangential (double) root.
    # * Event times closer together than Utils.precision are merged and treated as simultaneous;
    #   a single event time is reported for them.
    # * Events within Utils.precision of the interval boundaries are attributed to the boundary
    #   and are not reported as interior events.
    #
    # Consequently, hull configurations living on intervals shorter than Utils.precision may be
    # merged into their neighbors; anything longer is guaranteed to be detected.
    #
    # == Coincident points
    #
    # Points that coincide (permanently or only at some instants) are handled by the static
    # ConvexHull: nearly coincident points are collapsed to a single representative (the lowest
    # input index), so hull vertex lists never contain duplicates and permanently coincident
    # points never cause spurious events.
    #
    # == Usage
    #
    #   result = Geom2D::Algorithms::KineticConvexHull.run(points, velocities, 0.0, 10.0)
    #   result.events    # => [t1, t2, ...] event times inside (0.0, 10.0)
    #   result.intervals # => [#<Interval from:, to:, vertices:>, ...] covering [0.0, 10.0]
    #
    # +vertices+ are arrays of indices into the input arrays, ordered counterclockwise. Sampling
    # any time strictly inside an interval and computing the static convex hull there yields
    # exactly the interval's vertex list.
    class KineticConvexHull

      extend Utils

      # One open time interval (from, to) together with the ordered convex hull vertex indices
      # valid throughout it.
      Interval = Struct.new(:from, :to, :vertices)

      # The result of a kinetic convex hull analysis: the event times at which the hull changes
      # members and the Interval objects covering the analyzed time range.
      Result = Struct.new(:events, :intervals)

      # Creates a new kinetic convex hull analyzer.
      #
      # +points+ are the initial positions (Geom2D::Point instances or [x, y] arrays) and
      # +velocities+ the constant velocities in the same format. Both arrays must have the same
      # length.
      def initialize(points, velocities)
        if points.length != velocities.length
          raise ArgumentError, "Need as many velocities as points"
        end
        @px = points.map {|p| Geom2D::Point(p).x }
        @py = points.map {|p| Geom2D::Point(p).y }
        @vx = velocities.map {|v| Geom2D::Point(v).x }
        @vy = velocities.map {|v| Geom2D::Point(v).y }
      end

      # Runs the analysis for the time interval [t_start, t_end] and returns a Result.
      def analyze(t_start, t_end)
        raise ArgumentError, "t_start must be smaller than t_end" unless t_start < t_end

        bounds = [t_start, *event_times(t_start, t_end), t_end]
        intervals = []
        bounds.each_cons(2) do |from, to|
          vertices = hull_at((from + to) / 2.0)
          if intervals.empty? || intervals[-1].vertices != vertices
            intervals << Interval.new(from, to, vertices)
          else
            intervals[-1].to = to
          end
        end
        Result.new(intervals[1..].map(&:from), intervals)
      end

      # Runs the analysis for the time interval [t_start, t_end] and returns a Result.
      #
      # See #initialize for the format of +points+ and +velocities+.
      def self.run(points, velocities, t_start, t_end)
        new(points, velocities).analyze(t_start, t_end)
      end

      private

      # Returns the sorted, merged event times strictly inside (t_start, t_end).
      def event_times(t_start, t_end)
        eps = Utils.precision
        roots = []
        (0...@px.length).each do |i|
          ((i + 1)...@px.length).each do |j|
            ((j + 1)...@px.length).each do |k|
              orientation_roots(i, j, k).each do |t|
                roots << t if t > t_start + eps && t < t_end - eps
              end
            end
          end
        end
        roots.sort!

        merged = []
        roots.each do |t|
          merged << t if merged.empty? || t - merged[-1] >= eps
        end
        merged
      end

      # Returns the times at which the points i, j, k are collinear.
      #
      # The orientation of the triple is a quadratic a*t^2 + b*t + c in t.
      def orientation_roots(i, j, k)
        eps = Utils.precision
        dxj = @px[j] - @px[i]
        dyj = @py[j] - @py[i]
        dxk = @px[k] - @px[i]
        dyk = @py[k] - @py[i]
        wxj = @vx[j] - @vx[i]
        wyj = @vy[j] - @vy[i]
        wxk = @vx[k] - @vx[i]
        wyk = @vy[k] - @vy[i]

        a = wxj * wyk - wyj * wxk
        b = dxj * wyk + wxj * dyk - dyj * wxk - wyj * dxk
        c = dxj * dyk - dyj * dxk

        if a.abs < eps
          return [] if b.abs < eps # constantly collinear (or constant non-zero): no event
          return [-c / b]
        end

        discriminant = b * b - 4 * a * c
        return [] if discriminant < 0 && discriminant.abs >= eps
        return [-b / (2 * a)] if discriminant.abs < eps

        sq = Math.sqrt(discriminant)
        [(-b - sq) / (2 * a), (-b + sq) / (2 * a)]
      end

      # Returns the convex hull vertex indices of the points at time t.
      def hull_at(t)
        points = (0...@px.length).map do |i|
          Geom2D::Point.new(@px[i] + @vx[i] * t, @py[i] + @vy[i] * t)
        end
        ConvexHull.hull_indices(points)
      end

    end

  end
end
