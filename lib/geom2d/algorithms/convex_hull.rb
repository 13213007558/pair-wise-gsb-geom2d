# -*- frozen_string_literal: true -*-

require 'geom2d/algorithms'
require 'geom2d/utils'

module Geom2D
  module Algorithms

    # Computes the convex hull of a set of points using Andrew's monotone chain algorithm.
    #
    # The hull is returned as a list of indices into the input array, ordered counterclockwise
    # and canonically rotated to start with the lowest index, which makes it easy to track and
    # compare point identities (needed, for example, by the kinetic convex hull analysis in
    # KineticConvexHull).
    #
    # Details of the computation:
    #
    # * Points that are (nearly) equal with respect to Utils.precision are treated as a single
    #   point; the lowest index among them is used as representative. This means that permanently
    #   coincident points never produce duplicate or degenerate hull vertices.
    #
    # * Points lying on a hull edge (i.e. collinear with their neighbors, within Utils.precision)
    #   are not reported as vertices; only the extreme vertices of each edge are returned.
    #
    # * If all points are collinear, the hull consists of the two extreme points. If all points
    #   coincide, it consists of a single index.
    class ConvexHull

      extend Utils

      # Returns the indices of the convex hull vertices of the given points, ordered
      # counterclockwise and rotated to start with the lowest index.
      #
      # The +points+ argument needs to be an array of Geom2D::Point instances (or objects
      # responding to #x and #y).
      def self.hull_indices(points)
        sorted = points.each_with_index.sort_by {|point, _index| [point.x, point.y] }

        # Collapse (nearly) coincident points, keeping the lowest original index
        unique = []
        sorted.each do |point, index|
          if unique.empty? ||
              (point.x - unique[-1][0].x).abs >= Utils.precision ||
              (point.y - unique[-1][0].y).abs >= Utils.precision
            unique << [point, index]
          else
            unique[-1][1] = index if index < unique[-1][1]
          end
        end

        return canonical(unique.map {|_, index| index }) if unique.length <= 2

        cross = lambda do |o, a, b|
          (a[0].x - o[0].x) * (b[0].y - o[0].y) - (a[0].y - o[0].y) * (b[0].x - o[0].x)
        end

        build_half = lambda do |list|
          half = []
          list.each do |entry|
            while half.length >= 2 && float_compare(cross.call(half[-2], half[-1], entry), 0) <= 0
              half.pop
            end
            half << entry
          end
          half
        end

        lower = build_half.call(unique)
        upper = build_half.call(unique.reverse)

        canonical((lower[0...-1] + upper[0...-1]).map {|_, index| index })
      end

      # Returns the convex hull vertices of the given points as Geom2D::Point instances, ordered
      # counterclockwise.
      #
      # See hull_indices for details on how coincident and collinear points are handled.
      def self.run(points)
        hull_indices(points).map {|index| points[index] }
      end

      # Rotates the hull so that it starts with the lowest index, making the representation of a
      # given vertex set unique.
      def self.canonical(hull)
        return hull if hull.empty?
        hull.rotate(hull.index(hull.min))
      end
      private_class_method :canonical

    end

  end
end
