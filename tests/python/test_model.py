import math
import unittest

import sycomore
from sycomore.units import *

class TestModel(unittest.TestCase):
    def test_pulse(self):
        model = sycomore.como.Model(
            sycomore.Species(1*s, 0.1*s),
            sycomore.Magnetization(0, 0, 1),
            [["dummy", sycomore.TimeInterval(0*s)]])

        model.apply_pulse(sycomore.Pulse(41*deg, 27*deg))

        grid = model.magnetization()
        for index, _ in sycomore.GridScanner(grid.origin(), grid.shape()):
            if index == sycomore.Index(0):
                self.assertAlmostEqual(
                    grid[index].p , 0.210607912662250-0.413341301933443j)
                self.assertAlmostEqual(grid[index].z, 0.754709580222772)
                self.assertAlmostEqual(
                    grid[index].m, 0.210607912662250+0.413341301933443j)
            else:
                self.assertEqual(grid[index].p, 0)
                self.assertAlmostEqual(grid[index].z, 0)
                self.assertAlmostEqual(grid[index].m, 0)

    def test_time_interval(self):
        model = sycomore.como.Model(
            sycomore.Species(math.log(2)*Hz, math.log(2)*Hz),
            sycomore.Magnetization(0, 0, 1), [
                ["foo", sycomore.TimeInterval(1*s)],
                ["bar", sycomore.TimeInterval(1*s)]])

        model.apply_pulse(sycomore.Pulse(45*deg, 90*deg))

        model.apply_time_interval("foo")

        grid = model.magnetization()
        for index, _ in sycomore.GridScanner(grid.origin(), grid.shape()):
            if index == sycomore.Index(-1, 0):
                self.assertEqual(grid[index].p, 0)
                self.assertEqual(grid[index].z, 0)
                self.assertAlmostEqual(grid[index].m, 0.25)
            elif index == sycomore.Index(0, 0):
                self.assertEqual(grid[index].p, 0)
                self.assertEqual(grid[index].z, 0.5*(1+math.sqrt(2)/2))
                self.assertEqual(grid[index].m, 0)
            elif index == sycomore.Index(1, 0):
                self.assertAlmostEqual(grid[index].p, 0.25)
                self.assertEqual(grid[index].z, 0)
                self.assertEqual(grid[index].m, 0)
            else:
                self.assertEqual(grid[index].p , 0)
                self.assertAlmostEqual(grid[index].z, 0)
                self.assertAlmostEqual(grid[index].m, 0)

        model.apply_time_interval("bar")

        grid = model.magnetization()
        for index, _ in sycomore.GridScanner(grid.origin(), grid.shape()):
            if index == sycomore.Index(-1, -1):
                self.assertEqual(grid[index].p, 0)
                self.assertEqual(grid[index].z, 0)
                self.assertAlmostEqual(grid[index].m, 0.125)
            elif index == sycomore.Index(0, 0):
                self.assertEqual(grid[index].p, 0)
                self.assertEqual(grid[index].z, 0.5+0.25*(1+math.sqrt(2)/2))
                self.assertEqual(grid[index].m, 0)
            elif index == sycomore.Index(1, 1):
                self.assertAlmostEqual(grid[index].p, 0.125)
                self.assertEqual(grid[index].z, 0)
                self.assertEqual(grid[index].m, 0)
            else:
                self.assertEqual(grid[index].p , 0)
                self.assertAlmostEqual(grid[index].z, 0)
                self.assertAlmostEqual(grid[index].m, 0)

        isochromat = model.isochromat()
        self.assertAlmostEqual(isochromat[0], 0.125*math.sqrt(2))
        self.assertAlmostEqual(isochromat[1], 0)
        self.assertAlmostEqual(isochromat[2], 0.5+0.25*(1+math.sqrt(2)/2))

        isochromat = model.isochromat(
            {sycomore.Index(0,0), sycomore.Index(-1, -1)})
        self.assertAlmostEqual(isochromat[0], 0.125*math.sqrt(2)/2)
        self.assertAlmostEqual(isochromat[1], 0)
        self.assertAlmostEqual(isochromat[2], 0.5+0.25*(1+math.sqrt(2)/2))

    def test_diffusion(self):
        model = sycomore.como.Model(
            sycomore.Species(0*Hz, 0*Hz, 1*um*um/ms),
            sycomore.Magnetization(0, 0, 1), [
                ["foo", sycomore.TimeInterval(500*ms, 0.1*rad/um)]])

        model.apply_pulse(sycomore.Pulse(40*deg, 0*deg))
        model.apply_time_interval("foo")

        grid = model.magnetization()
        for index, _ in sycomore.GridScanner(grid.origin(), grid.shape()):
            if index == sycomore.Index(-1):
                self.assertEqual(grid[index].p, 0)
                self.assertEqual(grid[index].z, 0)
                self.assertAlmostEqual(grid[index].m, 0+0.003062528150606j)
            elif index == sycomore.Index(0):
                self.assertEqual(grid[index].p, 0)
                self.assertAlmostEqual(grid[index].z, 0.766044443118978)
                self.assertEqual(grid[index].m, 0)
            elif index == sycomore.Index(1):
                self.assertAlmostEqual(grid[index].p, 0-0.003062528150606j)
                self.assertEqual(grid[index].z, 0)
                self.assertEqual(grid[index].m, 0)
            else:
                self.assertEqual(grid[index].p , 0)
                self.assertAlmostEqual(grid[index].z, 0)
                self.assertAlmostEqual(grid[index].m, 0)


if __name__ == "__main__":
    unittest.main()
