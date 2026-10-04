# -*- frozen_string_literal: true -*-

require 'test_helper'
require 'geom2d/algorithms'

Packing = Geom2D::Algorithms::GuillotinePacking
describe Packing do
  def brute_best(width, height, available, orders, kerf)
    @brute_memo ||= {}
    key = [width, height, available, kerf]
    return @brute_memo[key] if @brute_memo.key?(key)

    candidates = [0]
    orders.each_with_index do |order, index|
      next if available[index].zero?

      if order[:width] == width && order[:height] == height
        used = available.dup
        used[index] -= 1
        candidates << order[:value]
      end

      if order.fetch(:rotatable, false) && order[:height] == width && order[:width] == height
        candidates << order[:value]
      end
    end

    each_used(available) do |first_used|
      remaining = available.each_with_index.map {|count, i| count - first_used[i] }

      (1..(width - kerf - 1)).each do |position|
        second_width = width - kerf - position
        candidates << brute_best(position, height, first_used, orders, kerf) +
                      brute_best(second_width, height, remaining, orders, kerf)
      end

      (1..(height - kerf - 1)).each do |position|
        second_height = height - kerf - position
        candidates << brute_best(width, position, first_used, orders, kerf) +
                      brute_best(width, second_height, remaining, orders, kerf)
      end
    end

    @brute_memo[key] = candidates.max
  end

  def each_used(limits)
    return yield(limits) if limits.empty?

    each_used(limits[1..-1]) do |tail|
      (0..limits.first).each {|count| yield([count] + tail) }
    end
  end

  before do
    @brute_memo = {}
  end

  it 'cuts pieces with an integer vertical saw kerf and validates the tree' do
    orders = [
      { width: 2, height: 2, quantity: 2, value: 10, rotatable: true }
    ]
    solution = Packing.solve(5, 4, orders, kerf_width: 1)
    solution.validate!

    assert_equal(20, solution.total_value)
    assert(solution.optimal?)
    assert_equal({ 0 => 2 }, solution.used_quantities)
    assert_equal([10, 10], solution.leaves.select(&:order?).map(&:value))
    assert(solution.cuts.any? {|cut| cut.vertical? && cut.position == 2 && cut.kerf_width == 1 })
    assert_equal(brute_best(5, 4, [2], orders, 1), solution.total_value)
  end

  it 'respects per-order rotation restrictions' do
    orders = [
      { id: :wide, width: 3, height: 2, quantity: 1, value: 10 },
      { id: :square, width: 2, height: 2, quantity: 1, value: 1 }
    ]
    solution = Packing.solve(2, 3, orders, kerf_width: 0)
    solution.validate!

    wide_leaf = solution.leaves.find {|leaf| leaf.order_id == :wide }
    assert_nil(wide_leaf)
    assert_equal(1, solution.total_value)
    assert_equal(brute_best(2, 3, [1, 1], orders, 0), solution.total_value)
  end

  it 'uses rotation when allowed' do
    orders = [
      { id: :wide, width: 3, height: 2, quantity: 1, value: 10, rotatable: true }
    ]
    solution = Packing.solve(2, 3, orders)
    solution.validate!

    leaf = solution.leaves.find(&:order?)
    assert_equal(:wide, leaf.order_id)
    assert_equal([2, 3], [leaf.width, leaf.height])
    assert(leaf.rotated?)
    assert_equal(10, solution.total_value)
    assert_equal(brute_best(2, 3, [1], orders, 0), solution.total_value)
  end

  it 'does not exceed order quantities on mixed small boards' do
    orders = [
      { id: :a, width: 2, height: 2, quantity: 2, value: 10 },
      { id: :b, width: 1, height: 3, quantity: 1, value: 7, rotatable: true }
    ]
    solution = Packing.solve(4, 3, orders, kerf_width: 1)
    solution.validate!

    solution.used_quantities.each do |id, count|
      order = orders.find {|item| item[:id] == id }
      assert(count <= order[:quantity])
    end
    assert_equal(brute_best(4, 3, [2, 1], orders, 1), solution.total_value)
    assert(solution.optimal?)
  end

  it 'leaves waste when no full order fits' do
    orders = [{ width: 6, height: 6, quantity: 1, value: 10 }]
    solution = Packing.solve(5, 5, orders)
    solution.validate!

    assert_equal(0, solution.total_value)
    assert(solution.leaves.one?(&:waste?))
    assert_equal({ 0 => 0 }, solution.used_quantities)
  end

  it 'returns a feasible non-optimal solution when the node budget is exhausted' do
    orders = [
      { width: 1, height: 1, quantity: 100, value: 1 }
    ]
    solution = Packing.solve(1, 1, orders, node_limit: 0)
    solution.validate!

    assert_equal(1, solution.total_value)
    refute(solution.optimal?)
    assert_equal({ 0 => 1 }, solution.used_quantities)
  end

  it 'compares exhaustive enumeration across several small boards' do
    orders = [
      { id: :a, width: 2, height: 1, quantity: 2, value: 5, rotatable: true },
      { id: :b, width: 1, height: 2, quantity: 1, value: 7 },
      { id: :c, width: 1, height: 1, quantity: 3, value: 2 }
    ]

    [3, 4].each do |width|
      [3, 4].each do |height|
        [0, 1].each do |kerf|
          next if width <= kerf + 1 && height <= kerf + 1

          solution = Packing.solve(width, height, orders, kerf_width: kerf)
          solution.validate!
          expected = brute_best(width, height, [2, 1, 3], orders, kerf)
          assert_equal(expected, solution.total_value,
                       "#{width}x#{height} kerf=#{kerf}")
        end
      end
    end
  end
end
