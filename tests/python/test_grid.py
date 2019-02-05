import unittest

import sycomore

class TestGrid(unittest.TestCase):
    def test_empty_constructor(self):
        grid = sycomore.Grid[complex]()
        self.assertEqual(grid.dimension(), 0)
        self.assertTrue(grid.origin().empty())
        self.assertTrue(grid.shape().empty())
        self.assertTrue(grid.stride().empty())

    def test_uninitialized_constructor(self):
        grid = sycomore.Grid[complex](
            sycomore.Index(-2, -3, -5), sycomore.Shape(5, 7, 11))
        self.assertEqual(grid.dimension(), 3)
        self.assertSequenceEqual(grid.origin(), [-2, -3, -5])
        self.assertSequenceEqual(grid.shape(), [5, 7, 11])
        self.assertSequenceEqual(grid.stride(), [1, 5, 35, 385])

    def test_initialized_constructor(self):
        grid = sycomore.Grid[complex](
            sycomore.Index(-5, -3), sycomore.Shape(11, 7), 1+2j)
        self.assertEqual(grid.dimension(), 2)
        self.assertSequenceEqual(grid.origin(), [-5, -3])
        self.assertSequenceEqual(grid.shape(), [11, 7])
        self.assertSequenceEqual(grid.stride(), [1, 11, 77])
        for x in grid:
            self.assertEqual(x, 1+2j)

    def test_accessor(self):
        grid = sycomore.Grid[complex](
            sycomore.Index(-5, -3), sycomore.Shape(11, 7), 0)
        grid[sycomore.Index(-2, -1)] = 1+2j;

        for index, offset in sycomore.GridScanner(grid.origin(), grid.shape()):
            self.assertEqual(
                grid[index], 1+2j if index == sycomore.Index(-2, -1) else 0)

    def test_scan_order(self):
        grid = sycomore.Grid[complex](
            sycomore.Index(-5, -3), sycomore.Shape(11, 7), 0)

        i = 0
        for index, _ in sycomore.GridScanner(grid.origin(), grid.shape()):
            grid[index] = i
            i += 1

        i = 0
        for x in grid:
            self.assertEqual(x, i)
            i += 1

if __name__ == "__main__":
    unittest.main()
