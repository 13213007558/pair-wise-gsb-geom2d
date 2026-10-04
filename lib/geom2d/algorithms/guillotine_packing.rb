# -*- frozen_string_literal: true -*-

module Geom2D
  module Algorithms

    # Finds an optimal guillotine cutting pattern for a single rectangular board.
    class GuillotinePacking

      # One ordered rectangle type.
      class Order
        attr_reader :id, :width, :height, :quantity, :value

        def initialize(id, width, height, quantity, value, rotatable)
          @id = id
          @width = width
          @height = height
          @quantity = quantity
          @value = value
          @rotatable = rotatable
        end

        def rotatable?
          @rotatable == true
        end
      end

      # Common state for a panel in a cutting tree.
      class Node
        attr_reader :x, :y, :width, :height, :value

        def initialize(x, y, width, height, value)
          @x = x
          @y = y
          @width = width
          @height = height
          @value = value
        end

        def leaf?
          !cut?
        end

        def order?
          false
        end

        def waste?
          false
        end

        def cut?
          false
        end
      end

      # A terminal panel containing either one ordered piece or waste.
      class Leaf < Node
        attr_reader :order_id, :rotated

        def initialize(x, y, width, height, value: 0, order_id: nil, rotated: false)
          super(x, y, width, height, value)
          @order_id = order_id
          @rotated = rotated
        end

        def order?
          !order_id.nil?
        end

        def waste?
          order_id.nil?
        end

        def rotated?
          @rotated == true
        end
      end

      # A guillotine cut and the two panels produced by it.
      class Cut < Node
        attr_reader :direction, :position, :kerf_width, :first_child, :second_child

        def initialize(direction, x, y, width, height, position, kerf_width,
                       first_child, second_child)
          super(x, y, width, height, first_child.value + second_child.value)
          @direction = direction
          @position = position
          @kerf_width = kerf_width
          @first_child = first_child
          @second_child = second_child
        end

        def children
          [first_child, second_child]
        end

        def vertical?
          direction == :vertical
        end

        def horizontal?
          direction == :horizontal
        end

        def cut?
          true
        end
      end

      # A complete cutting plan and optimization metadata.
      class Solution
        attr_reader :root, :total_value, :used_quantities

        def initialize(root, orders, proven_optimal, kerf_width)
          @root = root
          @orders = orders
          @proven_optimal = proven_optimal
          @kerf_width = kerf_width
          @total_value = root.value
          @used_quantities = orders.each_with_object({}) {|order, counts| counts[order.id] = 0 }
          each_leaf {|leaf| @used_quantities[leaf.order_id] += 1 if leaf.order? }
        end

        def optimal?
          @proven_optimal == true
        end
        alias proven_optimal? optimal?

        def leaves
          result = []
          each_leaf {|leaf| result << leaf }
          result
        end

        def cuts
          result = []
          walk(root) {|node| result << node if node.cut? }
          result
        end

        # Raises RuntimeError if the binary tree is geometrically or economically inconsistent.
        def validate!
          validate_node(root, 0, 0, root.width, root.height)
          counts = @orders.each_with_object({}) {|order, counts| counts[order.id] = 0 }
          value = 0

          each_leaf do |leaf|
            next if leaf.waste?
            order = @orders.find {|candidate| candidate.id == leaf.order_id }
            raise "unknown order #{leaf.order_id.inspect}" unless order

            counts[order.id] += 1
            raise "quantity exceeded for order #{order.id.inspect}" if counts[order.id] > order.quantity
            raise "wrong leaf value for order #{order.id.inspect}" unless leaf.value == order.value
            value += leaf.value
          end

          raise 'leaf quantities do not match solution metadata' unless counts == @used_quantities
          raise 'leaf values do not match solution total' unless value == @total_value

          self
        end

        def valid?
          validate!
          true
        rescue RuntimeError
          false
        end

        def each_leaf(&block)
          return enum_for(:each_leaf) unless block_given?
          walk(root) {|node| yield(node) if node.leaf? }
        end

        private

        def walk(node, &block)
          yield(node)
          node.children.each {|child| walk(child, &block) } if node.cut?
        end

        def validate_node(node, x, y, width, height)
          unless [node.x, node.y, node.width, node.height] == [x, y, width, height]
            raise "panel geometry mismatch at (#{x}, #{y}, #{width}, #{height})"
          end

          return validate_leaf(node) if node.leaf?
          raise 'all cuts must consume the configured integer kerf width' unless node.kerf_width == @kerf_width

          first = node.first_child
          second = node.second_child
          if node.vertical?
            unless width - node.kerf_width - node.position == second.width
              raise 'vertical cut children do not span the panel with the kerf'
            end
            validate_node(first, x, y, node.position, height)
            validate_node(second, x + node.position + @kerf_width, y, second.width, height)
          else
            unless height - node.kerf_width - node.position == second.height
              raise 'horizontal cut children do not span the panel with the kerf'
            end
            validate_node(first, x, y, width, node.position)
            validate_node(second, x, y + node.position + @kerf_width, width, second.height)
          end

          raise 'cut value is not the sum of its child values' unless node.value == first.value + second.value
        end

        def validate_leaf(node)
          return if node.waste? && node.value.zero?

          order = @orders.find {|candidate| candidate.id == node.order_id }
          raise "unknown order #{node.order_id.inspect}" unless order

          if node.rotated?
            unless order.rotatable? && [order.height, order.width] == [node.width, node.height]
              raise "order #{order.id.inspect} is not allowed in this orientation"
            end
          else
            unless [order.width, order.height] == [node.width, node.height]
              raise "ordered leaf #{order.id.inspect} has wrong dimensions"
            end
          end
        end
      end

      Candidate = Struct.new(:node, :used, :value)

      def self.solve(board_width, board_height, orders,
                     kerf_width: 0, time_limit: nil, node_limit: nil)
        new(board_width, board_height, orders,
            kerf_width: kerf_width, time_limit: time_limit, node_limit: node_limit).solve
      end

      def self.optimize(*arguments, **options)
        solve(*arguments, **options)
      end

      def initialize(board_width, board_height, orders,
                     kerf_width:, time_limit:, node_limit:)
        @board_width = board_width
        @board_height = board_height
        @orders = normalize_orders(orders)
        @kerf_width = kerf_width
        @time_limit = time_limit
        @node_limit = node_limit
        @start_time = Time.now
        @memo = {}
        @exhausted = false
        validate_input
      end

      def solve
        candidates = subproblem(@board_width, @board_height, @orders.map(&:quantity))
        best = candidates.max_by {|candidate| [candidate.value, candidate.used.sum] }
        Solution.new(best.node, @orders, !@exhausted, @kerf_width)
      end

      private

      def normalize_orders(orders)
        raise ArgumentError, 'orders must be an array' unless orders.is_a?(Array)

        orders.each_with_index.map do |attributes, index|
          raise ArgumentError, "order at index #{index} must be a hash" unless attributes.is_a?(Hash)

          data = attributes.each_with_object({}) {|(key, value), result| result[key.to_sym] = value }
          quantity = data.fetch(:quantity) { data.fetch(:max_quantity) { data.fetch(:limit) } }
          Order.new(data.fetch(:id, index), data.fetch(:width), data.fetch(:height),
                    quantity, data.fetch(:value), data.fetch(:rotatable, false))
        end
      end

      def validate_input
        unless @board_width.is_a?(Integer) && @board_width.positive?
          raise ArgumentError, 'board width must be a positive integer'
        end
        unless @board_height.is_a?(Integer) && @board_height.positive?
          raise ArgumentError, 'board height must be a positive integer'
        end
        unless @kerf_width.is_a?(Integer) && !@kerf_width.negative?
          raise ArgumentError, 'kerf width must be a non-negative integer'
        end
        unless @time_limit.nil? || (@time_limit.is_a?(Numeric) && @time_limit.positive?)
          raise ArgumentError, 'time limit must be nil or a positive number'
        end
        unless @node_limit.nil? || (@node_limit.is_a?(Integer) && @node_limit >= 0)
          raise ArgumentError, 'node limit must be nil or a non-negative integer'
        end

        if @orders.map(&:id).uniq.size != @orders.size
          raise ArgumentError, 'order ids must be unique'
        end

        @orders.each_with_index do |order, index|
          unless order.width.is_a?(Integer) && order.width.positive?
            raise ArgumentError, "order at index #{index} has invalid width"
          end
          unless order.height.is_a?(Integer) && order.height.positive?
            raise ArgumentError, "order at index #{index} has invalid height"
          end
          unless order.quantity.is_a?(Integer) && order.quantity.positive?
            raise ArgumentError, "order at index #{index} has invalid quantity"
          end
          unless order.value.is_a?(Numeric) && order.value.positive?
            raise ArgumentError, "order at index #{index} has invalid value"
          end
        end
      end

      def subproblem(width, height, available)
        key = [width, height, available]
        return @memo[key] if @memo.key?(key)

        candidates = [waste_candidate(width, height)]
        add_leaf_candidates(candidates, width, height, available)

        if budget_exhausted?
          @exhausted = true
          result = candidates
          @memo[key] = result
          return result
        end

        add_cut_candidates(candidates, :vertical, width, height, available)
        add_cut_candidates(candidates, :horizontal, width, height, available)
        result = nondominated(candidates)
        @memo[key] = result
        result
      end

      def budget_exhausted?
        return true if @exhausted
        return true if @node_limit && @memo.size >= @node_limit

        @time_limit && Time.now - @start_time >= @time_limit
      end

      def waste_candidate(width, height)
        Candidate.new(Leaf.new(0, 0, width, height), zeros, 0)
      end

      def add_leaf_candidates(candidates, width, height, available)
        @orders.each_with_index do |order, index|
          next if available[index].zero?

          if order.width == width && order.height == height
            used = zeros.dup
            used[index] = 1
            candidates << Candidate.new(ordered_leaf(order, width, height, false), used, order.value)
          end

          next unless order.rotatable? && order.height == width && order.width == height

          used = zeros.dup
          used[index] = 1
          candidates << Candidate.new(ordered_leaf(order, width, height, true), used, order.value)
        end
      end

      def ordered_leaf(order, width, height, rotated)
        Leaf.new(0, 0, width, height,
                 value: order.value, order_id: order.id, rotated: rotated)
      end

      def add_cut_candidates(candidates, direction, width, height, available)
        cut_dimension = direction == :vertical ? width : height
        cross_dimension = direction == :vertical ? height : width
        max_position = cut_dimension - @kerf_width - 1
        return if max_position < 1

        (1..max_position).each do |position|
          break if @exhausted
          second_size = cut_dimension - @kerf_width - position

          if direction == :vertical
            first_candidates = subproblem(position, cross_dimension, available)
          else
            first_candidates = subproblem(cross_dimension, position, available)
          end

          first_candidates.each do |first|
            break if @exhausted
            remaining = available.each_with_index.map {|count, index| count - first.used[index] }

            if direction == :vertical
              second_candidates = subproblem(second_size, cross_dimension, remaining)
            else
              second_candidates = subproblem(cross_dimension, second_size, remaining)
            end

            second_candidates.each do |second|
              used = first.used.each_with_index.map {|count, index| count + second.used[index] }
              node = build_cut(direction, width, height, position, first.node, second.node)
              candidates << Candidate.new(node, used, first.value + second.value)
            end
          end
        end
      end

      def build_cut(direction, width, height, position, first_node, second_node)
        if direction == :vertical
          first = first_node
          second = translate_node(second_node, position + @kerf_width, 0)
        else
          first = first_node
          second = translate_node(second_node, 0, position + @kerf_width)
        end

        Cut.new(direction, 0, 0, width, height, position, @kerf_width, first, second)
      end

      def translate_node(node, delta_x, delta_y)
        if node.cut?
          first = translate_node(node.first_child, delta_x, delta_y)
          second = translate_node(node.second_child, delta_x, delta_y)
          Cut.new(node.direction, node.x + delta_x, node.y + delta_y,
                  node.width, node.height, node.position, node.kerf_width, first, second)
        else
          Leaf.new(node.x + delta_x, node.y + delta_y, node.width, node.height,
                   value: node.value, order_id: node.order_id, rotated: node.rotated?)
        end
      end

      def nondominated(candidates)
        unique = candidates.uniq {|candidate| [candidate.used, candidate.value] }
        unique.reject do |candidate|
          unique.any? do |other|
            next false if other.equal?(candidate)
            next false unless other.value >= candidate.value
            next false unless other.used.each_with_index.all? {|count, index| count >= candidate.used[index] }

            other.value > candidate.value || other.used != candidate.used
          end
        end
      end

      def zeros
        @zeros ||= Array.new(@orders.size, 0)
      end
    end
  end
end
