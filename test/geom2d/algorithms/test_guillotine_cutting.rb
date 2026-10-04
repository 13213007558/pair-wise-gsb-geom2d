# -*- frozen_string_literal: true -*-

require 'test_helper'
require 'geom2d/algorithms/guillotine_cutting'

describe Geom2D::Algorithms::GuillotineCutting do
  GuillotineCutting = Geom2D::Algorithms::GuillotineCutting

  # Independent complete enumeration of all guillotine cutting patterns, used as the
  # reference for the solver on small sheets.
  def brute_force_value(width, height, kerf, items, remaining, memo = {})
    key = [width, height, remaining]
    return memo[key] if memo.key?(key)

    best = 0
    items.each_with_index do |item, index|
      next if remaining[index] == 0
      item.orientations.each do |item_width, item_height, _|
        best = item.value if item_width <= width && item_height <= height && item.value > best
      end
    end
    [[width, height], [height, width]].each_with_index do |(size, _), axis|
      1.upto(size - kerf - 1) do |offset|
        child_sizes = if axis == 0
                        [[offset, height], [size - offset - kerf, height]]
                      else
                        [[width, offset], [width, size - offset - kerf]]
                      end
        splits(remaining).each do |first_remaining|
          second_remaining = remaining.zip(first_remaining).map {|a, b| a - b }
          total = brute_force_value(*child_sizes[0], kerf, items, first_remaining, memo) +
            brute_force_value(*child_sizes[1], kerf, items, second_remaining, memo)
          best = total if total > best
        end
      end
    end
    memo[key] = best
  end

  def splits(remaining, &block)
    ranges = remaining.map {|count| (0..count).to_a }
    ranges[0].product(*ranges[1..-1]).each(&block)
  end

  def solve(width, height, kerf, item_specs, budget: 100_000)
    cutter = GuillotineCutting.new(width, height, kerf: kerf)
    item_specs.each do |spec|
      cutter.add_item(spec[0], spec[1], value: spec[2], limit: spec[3], rotatable: spec[4])
    end
    result = cutter.solve(budget: budget)
    assert(GuillotineCutting.verify(result.tree, cutter.items, kerf))
    [result, cutter]
  end

  def collect_leaves(node, leaves)
    if node.leaf?
      leaves << node
    else
      collect_leaves(node.first, leaves)
      collect_leaves(node.second, leaves)
    end
  end

  it 'places a single item and leaves the rest as waste' do
    result, = solve(10, 10, 0, [[4, 4, 7, 1, false]])
    assert_equal(7, result.value)
    assert(result.optimal?)
    assert_kind_of(GuillotineCutting::ItemNode, result.tree)
  end

  it 'packs multiple items up to their quantity limits' do
    result, = solve(10, 10, 0, [[4, 4, 5, 2, false]])
    assert_equal(10, result.value)
    assert(result.optimal?)
    leaves = []
    collect_leaves(result.tree, leaves)
    assert_equal(2, leaves.count {|leaf| leaf.kind_of?(GuillotineCutting::ItemNode) })
  end

  it 'accounts for the saw kerf between pieces' do
    # two 4x4 pieces need 4 + kerf + 4 <= 9 next to each other
    result, = solve(9, 4, 1, [[4, 4, 1, 2, false]])
    assert_equal(2, result.value)
    result, = solve(9, 4, 2, [[4, 4, 1, 2, false]])
    assert_equal(1, result.value)
  end

  it 'honours the per-item rotation allowance' do
    result, = solve(2, 4, 0, [[4, 2, 5, 1, false]])
    assert_equal(0, result.value)
    result, = solve(2, 4, 0, [[4, 2, 5, 1, true]])
    assert_equal(5, result.value)
    assert(result.tree.rotated?)
  end

  it 'allows leaving waste when nothing fits exactly' do
    # the item does not fill the sheet, the uncovered area is waste
    result, = solve(5, 5, 0, [[3, 3, 7, 1, false]])
    assert_equal(7, result.value)
    assert_kind_of(GuillotineCutting::ItemNode, result.tree)
    assert_equal(5, result.tree.width)

    # trees with explicit waste leaves are accepted by the verification
    cutter = GuillotineCutting.new(8, 4, kerf: 0)
    item = cutter.add_item(4, 4, value: 7, limit: 1)
    first = GuillotineCutting::ItemNode.new(0, 0, 4, 4, item, false)
    second = GuillotineCutting::WasteNode.new(4, 0, 4, 4)
    tree = GuillotineCutting::CutNode.new(0, 0, 8, 4, :vertical, 4, first, second)
    assert(GuillotineCutting.verify(tree, cutter.items, 0))
    assert_equal(7, tree.value)
  end

  it 'returns a feasible but unproven solution when the budget is exhausted' do
    result, cutter = solve(10, 10, 0, [[4, 4, 5, 4, false]], budget: 0)
    refute(result.optimal?)
    assert_equal(5, result.value) # only the leaf options of the root were evaluated
    assert(GuillotineCutting.verify(result.tree, cutter.items, 0))
  end

  it 'marks the solution as optimal when the budget suffices' do
    result, = solve(10, 10, 0, [[4, 4, 5, 4, false]])
    assert(result.optimal?)
    assert_equal(20, result.value)
  end

  it 'verifies tampered trees as invalid' do
    result, cutter = solve(9, 4, 1, [[4, 4, 1, 2, false]])
    cut = result.tree
    assert_kind_of(GuillotineCutting::CutNode, cut)
    # rebuild the tree without accounting for the kerf
    first = cut.first
    second = GuillotineCutting::ItemNode.new(cut.x + cut.offset, cut.y,
                                             cut.width - cut.offset, cut.height,
                                             cutter.items[0], false)
    tampered = GuillotineCutting::CutNode.new(cut.x, cut.y, cut.width, cut.height,
                                              cut.direction, cut.offset, first, second)
    assert_raises(GuillotineCutting::VerificationError) do
      GuillotineCutting.verify(tampered, cutter.items, 1)
    end
  end

  it 'rejects trees exceeding the quantity limits' do
    cutter = GuillotineCutting.new(8, 4, kerf: 0)
    item = cutter.add_item(4, 4, value: 1, limit: 1)
    first = GuillotineCutting::ItemNode.new(0, 0, 4, 4, item, false)
    second = GuillotineCutting::ItemNode.new(4, 0, 4, 4, item, false)
    tree = GuillotineCutting::CutNode.new(0, 0, 8, 4, :vertical, 4, first, second)
    assert_raises(GuillotineCutting::VerificationError) do
      GuillotineCutting.verify(tree, cutter.items, 0)
    end
  end

  it 'matches the complete enumeration on small sheets' do
    cases = [
      [6, 4, 1, [[2, 2, 3, 2, false], [3, 2, 4, 1, false]]],
      [6, 4, 0, [[2, 3, 5, 2, true], [2, 2, 3, 1, false]]],
      [5, 5, 1, [[3, 2, 4, 2, true], [2, 2, 2, 3, false]]],
      [7, 3, 0, [[3, 3, 6, 1, false], [2, 2, 3, 2, true]]],
      [4, 4, 2, [[2, 2, 5, 2, false], [1, 3, 3, 2, true]]],
    ]
    cases.each do |width, height, kerf, item_specs|
      result, cutter = solve(width, height, kerf, item_specs)
      expected = brute_force_value(width, height, kerf, cutter.items,
                                   cutter.items.map(&:limit))
      assert(result.optimal?, "expected a proven optimal result for #{item_specs.inspect}")
      assert_equal(expected, result.value,
                   "solver disagrees with complete enumeration for #{item_specs.inspect}")
    end
  end
end
