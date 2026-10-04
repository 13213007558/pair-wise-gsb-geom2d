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

module Geom2D
  module Algorithms

    # Exact solver for the guillotine cutting problem on a single rectangular stock sheet.
    #
    # A set of order items (rectangles with a value, a quantity limit and an individual rotation
    # allowance) has to be cut from one rectangular stock sheet such that the total value of the
    # placed items is maximized. Every cut is a guillotine cut, i.e. it runs straight across the
    # full width or height of the sub-sheet it is applied to, and each cut consumes a saw kerf of
    # integer width. Leftover waste is allowed.
    #
    # The solver uses dynamic programming over (width, height, remaining quantities) states and
    # is therefore only suited for small problem sizes. The total amount of explored
    # state/split combinations can be limited via a budget; when the budget is exhausted the
    # best feasible solution found so far is returned but the result is not marked as proven
    # optimal.
    #
    # The result is a binary cutting tree (see CutNode, ItemNode and WasteNode) that can be
    # executed cut by cut from the root downwards. GuillotineCutting.verify checks a tree for
    # geometric and order-related consistency.
    #
    # All dimensions (sheet, items, kerf) have to be integers.
    class GuillotineCutting

      # The default budget (number of explored state/split combinations) used by #solve.
      DEFAULT_BUDGET = 1_000_000

      # Raised when a cutting tree fails verification.
      class VerificationError < StandardError; end

      # Describes a single order item: a rectangle of a certain size and value of which at most
      # +limit+ pieces may be cut from the stock sheet.
      class Item

        # The width of the item.
        attr_reader :width

        # The height of the item.
        attr_reader :height

        # The value of a single placed piece.
        attr_reader :value

        # The maximum number of pieces that may be placed.
        attr_reader :limit

        def initialize(width, height, value: 1, limit: 1, rotatable: false)
          unless width.kind_of?(Integer) && height.kind_of?(Integer) && width > 0 && height > 0
            raise ArgumentError, 'width and height have to be positive integers'
          end
          unless limit.kind_of?(Integer) && limit >= 0
            raise ArgumentError, 'limit has to be a non-negative integer'
          end
          @width = width
          @height = height
          @value = value
          @limit = limit
          @rotatable = !!rotatable
        end

        # Returns +true+ if the item may be rotated by 90 degrees.
        def rotatable?
          @rotatable
        end

        # Returns the allowed orientations as list of [width, height, rotated] triples.
        def orientations
          @orientations ||= begin
            list = [[@width, @height, false]]
            list << [@height, @width, true] if @rotatable && @width != @height
            list
          end
        end

      end

      # Base class for the nodes of a cutting tree. Each node covers the rectangle with the
      # bottom-left corner (x, y) and the given width and height.
      class Node

        # The x-coordinate of the bottom-left corner of the covered rectangle.
        attr_reader :x

        # The y-coordinate of the bottom-left corner of the covered rectangle.
        attr_reader :y

        # The width of the covered rectangle.
        attr_reader :width

        # The height of the covered rectangle.
        attr_reader :height

        def initialize(x, y, width, height)
          @x = x
          @y = y
          @width = width
          @height = height
        end

        # Returns +true+ if this node is a leaf of the cutting tree.
        def leaf?
          true
        end

      end

      # A leaf of the cutting tree representing waste material.
      class WasteNode < Node

        # The value of waste is always zero.
        def value
          0
        end

      end

      # A leaf of the cutting tree representing a placed item. The covered rectangle may be
      # larger than the item itself; the surrounding area is waste.
      class ItemNode < Node

        # The placed Item.
        attr_reader :item

        def initialize(x, y, width, height, item, rotated)
          super(x, y, width, height)
          @item = item
          @rotated = rotated
        end

        # Returns +true+ if the item is placed rotated by 90 degrees.
        def rotated?
          @rotated
        end

        # The value of the placed item.
        def value
          @item.value
        end

        # The width of the placed item, taking the rotation into account.
        def item_width
          @rotated ? @item.height : @item.width
        end

        # The height of the placed item, taking the rotation into account.
        def item_height
          @rotated ? @item.width : @item.height
        end

      end

      # An inner node of the cutting tree representing a single guillotine cut. The cut runs
      # across the full height (direction :vertical) or width (direction :horizontal) of the
      # covered rectangle at the given offset from the node's origin and splits the node into
      # the +first+ and +second+ child, consuming +kerf+ units of material between them.
      class CutNode < Node

        # The direction of the cut, either :vertical or :horizontal.
        attr_reader :direction

        # The position of the cut as offset from the node's origin along the cut axis.
        attr_reader :offset

        # The child before the cut (left or bottom part).
        attr_reader :first

        # The child after the cut (right or top part).
        attr_reader :second

        def initialize(x, y, width, height, direction, offset, first, second)
          super(x, y, width, height)
          @direction = direction
          @offset = offset
          @first = first
          @second = second
        end

        # Returns +false+ since a cut node is an inner node.
        def leaf?
          false
        end

        # The combined value of both children.
        def value
          @first.value + @second.value
        end

      end

      # The result of a call to GuillotineCutting#solve.
      class Result

        # The root node of the binary cutting tree.
        attr_reader :tree

        # The number of states that were fully explored during solving.
        attr_reader :states

        def initialize(tree, optimal, states)
          @tree = tree
          @optimal = optimal
          @states = states
        end

        # The total value of all placed items.
        def value
          @tree.value
        end

        # Returns +true+ if the solution is proven to be optimal, i.e. the search was not
        # cut short by the budget.
        def optimal?
          @optimal
        end

      end

      # The width of the stock sheet.
      attr_reader :width

      # The height of the stock sheet.
      attr_reader :height

      # The width of the saw kerf consumed by every cut.
      attr_reader :kerf

      # The order items added via #add_item.
      attr_reader :items

      def initialize(width, height, kerf: 0)
        unless width.kind_of?(Integer) && height.kind_of?(Integer) && width > 0 && height > 0
          raise ArgumentError, 'width and height have to be positive integers'
        end
        unless kerf.kind_of?(Integer) && kerf >= 0
          raise ArgumentError, 'kerf has to be a non-negative integer'
        end
        @width = width
        @height = height
        @kerf = kerf
        @items = []
      end

      # Adds an order item of the given size and returns it.
      #
      # At most +limit+ pieces of the item are placed, each contributing +value+ to the total.
      # If +rotatable+ is true, the item may additionally be placed rotated by 90 degrees.
      def add_item(width, height, value: 1, limit: 1, rotatable: false)
        item = Item.new(width, height, value: value, limit: limit, rotatable: rotatable)
        @items << item
        item
      end

      # Solves the cutting problem and returns a Result with the cutting tree.
      #
      # At most +budget+ state/split combinations are explored. If the budget is exhausted
      # before the search is complete, the returned result still holds a feasible cutting tree
      # but Result#optimal? is false.
      def solve(budget: DEFAULT_BUDGET)
        unless budget.kind_of?(Integer) && budget >= 0
          raise ArgumentError, 'budget has to be a non-negative integer'
        end
        @budget = budget
        @work = 0
        @states = 0
        @proven = true
        @memo = {}
        @split_cache = {}
        @fits_cache = {}
        remaining = @items.map(&:limit)
        best(@width, @height, remaining)
        tree = build_node(0, 0, @width, @height, remaining)
        Result.new(tree, @proven, @states)
      end

      # Verifies that the cutting tree rooted at +tree+ is geometrically executable and
      # consistent with the given order +items+ and +kerf+.
      #
      # Checks performed along the tree:
      #
      # * every cut spans the full width or height of the sub-sheet it is applied to and the
      #   children exactly tile the parent minus the saw kerf,
      # * every placed item fits into its leaf rectangle in its declared orientation,
      # * the quantity limits of the order items are respected and
      # * the values of the leaves sum up to the value of the root node.
      #
      # Returns +true+ if the tree is valid and raises a VerificationError otherwise.
      def self.verify(tree, items, kerf)
        counts = Hash.new(0)
        value = verify_node(tree, items, kerf, counts)
        counts.each do |index, count|
          if count > items[index].limit
            raise VerificationError,
              "quantity limit exceeded for item #{index}: #{count} > #{items[index].limit}"
          end
        end
        if value != tree.value
          raise VerificationError, "leaf values sum to #{value} but root value is #{tree.value}"
        end
        true
      end

      # Verifies a single node and returns the sum of the values of its leaves.
      def self.verify_node(node, items, kerf, counts) # :nodoc:
        unless node.width.kind_of?(Integer) && node.height.kind_of?(Integer) &&
            node.width > 0 && node.height > 0
          raise VerificationError, 'node covers an invalid rectangle'
        end
        case node
        when WasteNode
          0
        when ItemNode
          index = items.index(node.item)
          raise VerificationError, 'placed item is not part of the order' unless index
          unless node.item_width <= node.width && node.item_height <= node.height
            raise VerificationError,
              "item does not fit into its leaf rectangle at " \
              "(#{node.x}, #{node.y}, #{node.width}, #{node.height})"
          end
          counts[index] += 1
          node.value
        when CutNode
          verify_cut(node, kerf)
          verify_node(node.first, items, kerf, counts) +
            verify_node(node.second, items, kerf, counts)
        else
          raise VerificationError, "unknown node type #{node.class}"
        end
      end
      private_class_method :verify_node

      # Verifies the geometry of a single guillotine cut.
      def self.verify_cut(node, kerf) # :nodoc:
        first, second = node.first, node.second
        case node.direction
        when :vertical
          expected = [[node.x, node.y, node.offset, node.height],
                      [node.x + node.offset + kerf, node.y,
                       node.width - node.offset - kerf, node.height]]
        when :horizontal
          expected = [[node.x, node.y, node.width, node.offset],
                      [node.x, node.y + node.offset + kerf,
                       node.width, node.height - node.offset - kerf]]
        else
          raise VerificationError, "unknown cut direction #{node.direction.inspect}"
        end
        if node.offset < 1 || expected[1][2] < 1 || expected[1][3] < 1
          raise VerificationError, 'cut leaves a non-positive sub-sheet'
        end
        actual = [[first.x, first.y, first.width, first.height],
                  [second.x, second.y, second.width, second.height]]
        if actual != expected
          raise VerificationError,
            "children #{actual.inspect} do not tile parent " \
            "(#{node.x}, #{node.y}, #{node.width}, #{node.height}) " \
            "cut #{node.direction} at #{node.offset} with kerf #{kerf}"
        end
      end
      private_class_method :verify_cut

      private

      # Returns the [value, decision] pair for the best solution of the sub-sheet (width,
      # height) with the given remaining item quantities.
      #
      # A decision is one of:
      #
      # * [:waste] - the whole sub-sheet is waste
      # * [:item, index, rotated] - a single item is placed, the rest is waste
      # * [:cut, direction, offset, first_remaining] - a guillotine cut, the remaining
      #   quantities of the second child are the complement of first_remaining
      def best(width, height, remaining)
        remaining = cap_remaining(width, height, remaining)
        key = [width, height, remaining]
        entry = @memo[key]
        return entry if entry

        value = 0
        decision = [:waste]
        @items.each_with_index do |item, index|
          next if remaining[index] == 0 || item.value <= value
          item.orientations.each do |item_width, item_height, rotated|
            next unless item_width <= width && item_height <= height
            value = item.value
            decision = [:item, index, rotated]
            break
          end
        end

        if @work < @budget
          @states += 1
          value, decision = explore_cuts(width, height, remaining, value, decision)
        else
          @proven = false
        end

        @memo[key] = [value, decision]
      end

      # Explores all guillotine cuts of the sub-sheet (width, height) and returns the best
      # [value, decision] pair found, considering the given incumbent.
      #
      # Stops early, marking the result as not proven optimal, once the work budget is
      # exhausted.
      def explore_cuts(width, height, remaining, value, decision)
        [[:vertical, width, height], [:horizontal, height, width]].each do |direction, size, across|
          1.upto(size - @kerf - 1) do |offset|
            if direction == :vertical
              first_width, first_height = offset, across
              second_width, second_height = size - offset - @kerf, across
            else
              first_width, first_height = across, offset
              second_width, second_height = across, size - offset - @kerf
            end
            first_fits = fits_any?(first_width, first_height, remaining)
            second_fits = fits_any?(second_width, second_height, remaining)
            next unless first_fits || second_fits

            if first_fits && second_fits
              splits(remaining).each do |first_remaining|
                return [value, decision] if budget_exhausted?
                second_remaining = remaining.zip(first_remaining).map {|a, b| a - b }
                total = best(first_width, first_height, first_remaining)[0] +
                  best(second_width, second_height, second_remaining)[0]
                if total > value
                  value = total
                  decision = [:cut, direction, offset, first_remaining]
                end
              end
            else
              return [value, decision] if budget_exhausted?
              first_remaining = (first_fits ? remaining : remaining.map { 0 })
              second_remaining = (first_fits ? remaining.map { 0 } : remaining)
              total = best(first_width, first_height, first_remaining)[0] +
                best(second_width, second_height, second_remaining)[0]
              if total > value
                value = total
                decision = [:cut, direction, offset, first_remaining]
              end
            end
          end
        end
        [value, decision]
      end

      # Accounts one unit of work and returns +true+ once the budget is exhausted, marking
      # the result as not proven optimal in that case.
      def budget_exhausted?
        @work += 1
        @proven = false if @work >= @budget
        @work >= @budget
      end

      # Returns all vectors first_remaining with 0 <= first_remaining[i] <= remaining[i].
      def splits(remaining)
        @split_cache[remaining] ||= begin
          ranges = remaining.map {|count| (0..count).to_a }
          ranges[0].product(*ranges[1..-1])
        end
      end

      # Returns +true+ if any remaining item fits into a sub-sheet of the given size.
      def fits_any?(width, height, remaining)
        @fits_cache[[width, height, remaining]] ||= @items.each_with_index.any? do |item, index|
          remaining[index] > 0 &&
            item.orientations.any? {|item_width, item_height, _|
              item_width <= width && item_height <= height
            }
        end
      end

      # Caps the remaining quantities at the number of pieces that can geometrically fit into a
      # sub-sheet of the given size. This increases the effectiveness of the memoization
      # without changing the result.
      def cap_remaining(width, height, remaining)
        @items.each_with_index.map do |item, index|
          next 0 if remaining[index] == 0
          max_fit = item.orientations.map do |item_width, item_height, _|
            (width / item_width) * (height / item_height)
          end.max
          [remaining[index], max_fit].min
        end
      end

      # Builds the cutting tree node for the sub-sheet with the bottom-left corner (x, y) using
      # the decisions stored by #best.
      def build_node(x, y, width, height, remaining)
        remaining = cap_remaining(width, height, remaining)
        _, decision = @memo[[width, height, remaining]]
        case decision[0]
        when :waste
          WasteNode.new(x, y, width, height)
        when :item
          ItemNode.new(x, y, width, height, @items[decision[1]], decision[2])
        when :cut
          _, direction, offset, first_remaining = decision
          second_remaining = remaining.zip(first_remaining).map {|a, b| a - b }
          if direction == :vertical
            first = build_node(x, y, offset, height, first_remaining)
            second = build_node(x + offset + @kerf, y, width - offset - @kerf, height,
                                second_remaining)
          else
            first = build_node(x, y, width, offset, first_remaining)
            second = build_node(x, y + offset + @kerf, width, height - offset - @kerf,
                                second_remaining)
          end
          CutNode.new(x, y, width, height, direction, offset, first, second)
        end
      end

    end

  end
end
