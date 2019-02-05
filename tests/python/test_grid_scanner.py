import unittest

import sycomore

class TestGridScanner(unittest.TestCase):
    def test_without_region(self):
        scanner = sycomore.GridScanner(
            sycomore.Index(-3, -5, -7), sycomore.Shape(7, 11, 15))
        data = range(7*11*15)

        iterator = iter(scanner)

        for z in range(-7, 8):
            for y in range(-5, 6):
                for x in range(-3, 4):
                    index, offset = next(iterator)
                    self.assertSequenceEqual(index, [x,y,z])
                    self.assertEqual(offset, data[offset])
        with self.assertRaises(StopIteration):
            next(iterator)

    def test_with_region_3d(self):
        scanner = sycomore.GridScanner(
            sycomore.Index(-3, -5, -7), sycomore.Shape(7, 11, 15),
            sycomore.Index(-1, -2, -3), sycomore.Shape(5, 7, 11))
        data = range(7*11*15)

        iterator = iter(scanner)

        for z in range(-3, 8):
            for y in range(-2, 5):
                for x in range(-1, 4):
                    index, offset = next(iterator)
                    self.assertSequenceEqual(index, [x,y,z])
                    self.assertEqual(offset, data[offset])
        with self.assertRaises(StopIteration):
            next(iterator)

    def test_with_region_2d(self):
        scanner = sycomore.GridScanner(
            sycomore.Index(-3, -5, -7), sycomore.Shape(7, 11, 15),
            sycomore.Index(-1, -2, -3), sycomore.Shape(5, 1, 11))
        data = range(7*11*15)

        iterator = iter(scanner)

        for z in range(-3, 8):
            for y in range(-2, -1):
                for x in range(-1, 4):
                    index, offset = next(iterator)
                    self.assertSequenceEqual(index, [x,y,z])
                    self.assertEqual(offset, data[offset])
        with self.assertRaises(StopIteration):
            next(iterator)


if __name__ == "__main__":
    unittest.main()
